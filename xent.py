#!/usr/bin/env python3
"""Compare token-wise HyenaDNA cross-entropy before and after head fine-tuning."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn.functional as F
from matplotlib.patches import Patch
from transformers import AutoModelForCausalLM, AutoTokenizer


DEFAULT_MODEL = "LongSafari/hyenadna-tiny-1k-seqlen-hf"
DEFAULT_FASTA = Path("data/arabiopsis-thaliana.CCA1.fasta")
DEFAULT_PROJECTION = Path("projection_head_finetuned.pt")
DEFAULT_OUTPUT = Path("data/hyena_base_finetuned_xent_comparison.png")
DEFAULT_VALIDATION_FRACTION = 0.1
NUCLEOTIDES = "ACGT"

# This text is deliberately kept in one place so it is easy to revise.
FIGURE_EXPLANATION = (
    "Each letter is coloured by the negative log-probability (cross-entropy, in nats) "
    "assigned to that nucleotide given the preceding sequence. Lower values mean the "
    "nucleotide was more predictable. The first nucleotide is grey because it has no "
    "preceding context. Both heads are evaluated over the same A/C/G/T outcome space; "
    "colours above the displayed limit are clipped."
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot per-nucleotide cross-entropy for base and fine-tuned HyenaDNA heads."
    )
    parser.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    parser.add_argument("--projection-head", type=Path, default=DEFAULT_PROJECTION)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--model", default=DEFAULT_MODEL)
    parser.add_argument(
        "--max-tokens",
        type=int,
        default=None,
        metavar="N",
        help="Limit analysis to the first N nucleotides of the selected region (default: all).",
    )
    parser.add_argument(
        "--evaluation-region",
        choices=("validation", "all"),
        default="validation",
        help="Evaluate the held-out validation tail or the entire sequence (default: validation).",
    )
    parser.add_argument(
        "--validation-fraction",
        type=float,
        default=DEFAULT_VALIDATION_FRACTION,
        help="Must match training when evaluating the validation region (default: 0.1).",
    )
    parser.add_argument(
        "--device",
        choices=("cpu", "cuda", "mps", "auto"),
        default="cpu",
        help="Inference device (default: cpu; 'auto' prefers CUDA, then MPS).",
    )
    parser.add_argument(
        "--context-size",
        type=int,
        default=512,
        help="Maximum bases in each model window (default: 512, matching training).",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=256,
        help="New bases scored per overlapping window after the first (default: 256).",
    )
    parser.add_argument("--line-width", type=int, default=70)
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument(
        "--colour-percentile",
        type=float,
        default=99.0,
        help="Percentile used as the shared colour maximum (default: 99).",
    )
    return parser.parse_args()


def read_fasta(path: Path) -> str:
    """Read one or more FASTA records as a single uppercase A/C/G/T sequence."""
    raw = path.read_text(encoding="utf-8")
    sequence = "".join(
        line.strip() for line in raw.splitlines() if line.strip() and not line.startswith(">")
    ).upper()
    invalid = sorted(set(re.sub(r"\s", "", sequence)) - set(NUCLEOTIDES))
    if invalid:
        raise ValueError(
            f"{path} contains unsupported sequence symbols: {', '.join(invalid)}. "
            "This comparison expects only A, C, G and T."
        )
    if len(sequence) < 2:
        raise ValueError(f"{path} must contain at least two nucleotides.")
    return sequence


def select_evaluation_region(
    sequence: str,
    region: str,
    context_size: int,
    validation_fraction: float,
) -> tuple[str, str]:
    """Select the same contiguous validation tail reserved by the trainer."""
    if region == "all":
        return sequence, "complete sequence"
    if not 0.0 < validation_fraction < 1.0:
        raise ValueError("--validation-fraction must be between 0 and 1.")
    validation_size = max(context_size, math.ceil(len(sequence) * validation_fraction))
    training_size = len(sequence) - validation_size
    if training_size < context_size:
        raise ValueError(
            f"Sequence length {len(sequence)} is too short for non-overlapping training "
            f"and validation regions with context size {context_size}."
        )
    return sequence[training_size:], f"held-out validation bases {training_size + 1}-{len(sequence)}"


def choose_device(requested: str) -> torch.device:
    if requested != "auto":
        device = torch.device(requested)
        if requested == "cuda" and not torch.cuda.is_available():
            raise RuntimeError("CUDA was requested but is not available.")
        if requested == "mps" and not torch.backends.mps.is_available():
            raise RuntimeError("MPS was requested but is not available.")
        return device
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def load_projection(path: Path, hidden_size: int, device: torch.device) -> torch.nn.Linear:
    projection = torch.nn.Linear(hidden_size, len(NUCLEOTIDES), bias=False)
    try:
        state = torch.load(path, map_location="cpu", weights_only=True)
    except TypeError:  # Compatibility with PyTorch versions predating weights_only.
        state = torch.load(path, map_location="cpu")
    projection.load_state_dict(state)
    projection.to(device).eval()
    return projection


def encode_sequence(tokenizer, sequence: str) -> tuple[torch.Tensor, list[int]]:
    nucleotide_ids = [
        tokenizer(base, add_special_tokens=False)["input_ids"][0] for base in NUCLEOTIDES
    ]
    # Converting character tokens directly avoids asking the tokenizer to accept a
    # sequence longer than the model limit; inference itself is windowed below.
    input_ids = torch.tensor(tokenizer.convert_tokens_to_ids(list(sequence)), dtype=torch.long)
    expected_ids = torch.tensor([nucleotide_ids[NUCLEOTIDES.index(base)] for base in sequence])
    if not torch.equal(input_ids.cpu(), expected_ids):
        raise RuntimeError("Tokenizer output is not one token per nucleotide as expected.")
    return input_ids, nucleotide_ids


def compute_cross_entropies(
    sequence: str,
    input_ids: torch.Tensor,
    nucleotide_ids: list[int],
    causal_model,
    projection: torch.nn.Linear,
    device: torch.device,
    context_size: int,
    stride: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Score every target once, retaining overlap as left context between windows."""
    if context_size < 2:
        raise ValueError("--context-size must be at least 2.")
    if not 1 <= stride < context_size:
        raise ValueError("--stride must be at least 1 and smaller than --context-size.")

    configured_max = int(getattr(causal_model.config, "max_seq_len", context_size))
    if context_size > configured_max:
        raise ValueError(
            f"--context-size {context_size} exceeds the model maximum {configured_max}."
        )

    base_scores = np.full(len(sequence), np.nan, dtype=np.float32)
    tuned_scores = np.full(len(sequence), np.nan, dtype=np.float32)
    output_head = causal_model.get_output_embeddings()
    backbone = causal_model.base_model
    input_ids = input_ids.to(device)
    nucleotide_index = torch.tensor(nucleotide_ids, device=device)
    class_by_base = torch.tensor(
        [NUCLEOTIDES.index(base) for base in sequence], dtype=torch.long, device=device
    )

    cursor = 1  # Position zero has no left context and is intentionally left as NaN.
    window_number = 0
    with torch.inference_mode():
        while cursor < len(sequence):
            if cursor == 1:
                end = min(len(sequence), context_size)
            else:
                end = min(len(sequence), cursor + stride)
            start = max(0, end - context_size)
            window = input_ids[start:end].unsqueeze(0)

            hidden = backbone(window, return_dict=True).last_hidden_state
            base_logits = output_head(hidden).index_select(-1, nucleotide_index)
            tuned_logits = projection(hidden)

            targets = torch.arange(cursor, end, device=device)
            prediction_positions = targets - start - 1
            labels = class_by_base.index_select(0, targets)
            base_loss = F.cross_entropy(
                base_logits[0].index_select(0, prediction_positions), labels, reduction="none"
            )
            tuned_loss = F.cross_entropy(
                tuned_logits[0].index_select(0, prediction_positions), labels, reduction="none"
            )
            base_scores[cursor:end] = base_loss.float().cpu().numpy()
            tuned_scores[cursor:end] = tuned_loss.float().cpu().numpy()

            window_number += 1
            print(
                f"  window {window_number}: scored bases {cursor + 1:,}-{end:,} "
                f"using input bases {start + 1:,}-{end:,}"
            )
            cursor = end

    return base_scores, tuned_scores


def summary(scores: np.ndarray) -> tuple[float, float]:
    mean = float(np.nanmean(scores))
    return mean, math.exp(mean)


def draw_sequence_panel(
    ax: plt.Axes,
    sequence: str,
    scores: np.ndarray,
    title: str,
    line_width: int,
    font_size: float,
    cmap: mpl.colors.Colormap,
    norm: mpl.colors.Normalize,
) -> None:
    rows = math.ceil(len(sequence) / line_width)
    for index, (base, score) in enumerate(zip(sequence, scores, strict=True)):
        colour = "#a8a8a8" if np.isnan(score) else cmap(norm(float(score)))
        ax.text(
            index % line_width,
            index // line_width,
            base,
            color=colour,
            family="monospace",
            fontsize=font_size,
            fontweight="bold",
            ha="center",
            va="center",
        )
    ax.set_xlim(-1, line_width)
    ax.set_ylim(rows - 0.25, -1.0)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    ax.set_title(title, fontsize=12, pad=12)


def make_figure(
    sequence: str,
    base_scores: np.ndarray,
    tuned_scores: np.ndarray,
    output: Path,
    line_width: int,
    dpi: int,
    colour_percentile: float,
    region_label: str,
) -> None:
    if line_width < 10:
        raise ValueError("--line-width must be at least 10.")
    if not 50.0 <= colour_percentile <= 100.0:
        raise ValueError("--colour-percentile must be between 50 and 100.")

    # For short selections, narrower wrapping lets each nucleotide occupy more
    # space. The cap preserves the existing dense layout for long sequences.
    display_line_width = min(
        line_width,
        max(10, math.ceil(math.sqrt(2.0 * len(sequence)))),
    )
    font_size = min(18.0, 6.6 * line_width / display_line_width)

    all_scores = np.concatenate((base_scores[1:], tuned_scores[1:]))
    colour_max = float(np.percentile(all_scores, colour_percentile))
    norm = mpl.colors.Normalize(vmin=0.0, vmax=max(colour_max, 1e-6), clip=True)
    cmap = mpl.colormaps["coolwarm"]
    rows = math.ceil(len(sequence) / display_line_width)
    figure_height = max(6.0, rows * 0.13 + 3.0)
    fig, axes = plt.subplots(1, 2, figsize=(20, figure_height))

    base_mean, base_ppl = summary(base_scores)
    tuned_mean, tuned_ppl = summary(tuned_scores)
    draw_sequence_panel(
        axes[0], sequence, base_scores,
        f"Pretrained HyenaDNA head\nmean CE {base_mean:.3f} nats · perplexity {base_ppl:.2f}",
        display_line_width, font_size, cmap, norm,
    )
    draw_sequence_panel(
        axes[1], sequence, tuned_scores,
        f"Fine-tuned A/C/G/T projection head\nmean CE {tuned_mean:.3f} nats · perplexity {tuned_ppl:.2f}",
        display_line_width, font_size, cmap, norm,
    )

    fig.suptitle(
        f"Arabidopsis thaliana: token-wise HyenaDNA cross-entropy\n{region_label}",
        fontsize=17,
        fontweight="bold",
        y=0.975,
    )
    fig.subplots_adjust(left=0.025, right=0.91, top=0.80, bottom=0.18, wspace=0.06)
    colourbar_ax = fig.add_axes((0.925, 0.24, 0.014, 0.57))
    colourbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=colourbar_ax)
    colourbar.set_label("Cross-entropy (nats; lower is more predictable)", fontsize=10)
    fig.legend(
        handles=[Patch(facecolor="#a8a8a8", label="No preceding context")],
        loc="lower center",
        bbox_to_anchor=(0.47, 0.083),
        frameon=False,
    )
    fig.text(0.47, 0.045, FIGURE_EXPLANATION, ha="center", va="center", fontsize=9, wrap=True)

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    complete_sequence = read_fasta(args.fasta)
    sequence, region_label = select_evaluation_region(
        complete_sequence,
        args.evaluation_region,
        args.context_size,
        args.validation_fraction,
    )
    if args.max_tokens is not None:
        if args.max_tokens < 2:
            raise ValueError("--max-tokens must be at least 2.")
        sequence = sequence[: args.max_tokens]
        region_label += f" (first {len(sequence)} selected bases)"
    device = choose_device(args.device)
    print(f"Loaded {len(sequence):,} nucleotides from {args.fasta}")
    print(f"Using device: {device}")
    print(f"Loading {args.model} ...")

    tokenizer = AutoTokenizer.from_pretrained(args.model, trust_remote_code=True)
    causal_model = AutoModelForCausalLM.from_pretrained(
        args.model, trust_remote_code=True, return_dict=True
    ).to(device).eval()
    hidden_size = int(causal_model.config.d_model)
    projection = load_projection(args.projection_head, hidden_size, device)
    input_ids, nucleotide_ids = encode_sequence(tokenizer, sequence)

    print("Computing token-wise cross-entropies ...")
    base_scores, tuned_scores = compute_cross_entropies(
        sequence,
        input_ids,
        nucleotide_ids,
        causal_model,
        projection,
        device,
        args.context_size,
        args.stride,
    )
    base_mean, base_ppl = summary(base_scores)
    tuned_mean, tuned_ppl = summary(tuned_scores)
    print(f"Base:       mean CE={base_mean:.4f} nats, perplexity={base_ppl:.4f}")
    print(f"Fine-tuned: mean CE={tuned_mean:.4f} nats, perplexity={tuned_ppl:.4f}")

    make_figure(
        sequence,
        base_scores,
        tuned_scores,
        args.output,
        args.line_width,
        args.dpi,
        args.colour_percentile,
        region_label,
    )
    print(f"Saved figure to {args.output}")


if __name__ == "__main__":
    main()
