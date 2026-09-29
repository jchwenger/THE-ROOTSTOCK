#!/usr/bin/env bash
set -euo pipefail

# Train the same A/C/G/T projection head across the official HyenaDNA model sizes.
# Models run sequentially so only one frozen backbone occupies memory at a time.
models=(
  "LongSafari/hyenadna-tiny-1k-seqlen-hf"
  "LongSafari/hyenadna-tiny-1k-seqlen-d256-hf"
  "LongSafari/hyenadna-tiny-16k-seqlen-d128-hf"
  "LongSafari/hyenadna-small-32k-seqlen-hf"
  "LongSafari/hyenadna-medium-160k-seqlen-hf"
  "LongSafari/hyenadna-medium-450k-seqlen-hf"
  "LongSafari/hyenadna-large-1m-seqlen-hf"
)

mkdir -p heads runs

for model in "${models[@]}"; do
  model_id="${model##*/}"
  model_data="${model_id#hyenadna-}"
  model_data="${model_data/-seqlen-/-}"
  model_data="${model_data%-hf}"
  log_path="runs/train_hyenadna_${model_data}.log"
  metrics_path="runs/train_hyenadna_${model_data}.metrics.json"

  echo "Training ${model}"
  echo "Logging to ${log_path}"
  uv run --frozen python train_projection_head.py \
    --model-name "${model}" \
    --log-file "${log_path}" \
    --metrics-file "${metrics_path}" \
    --fasta-files 'data/circadian/*.fasta' \
    --fasta-test-files 'data/circadian-held-out/*.fasta' \
    --fasta-test-files 'data/non-circadian/*.fasta' \
    --stride 16 \
    --device auto
done
