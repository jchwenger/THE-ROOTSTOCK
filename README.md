# THE ROOTSTOCK

**Yvonne Wang**

> *A bioart installation in which the body of a living plant writes its own poem.*

*Arabidopsis thaliana* — the model organism of plant molecular biology — speaks through its genome. A piezoelectric sensor reads vibration from the plant's stem. That signal passes through a genomic language model (HyenaDNA), translating the CCA1 circadian clock gene into a stream of English words that accumulate, drift, and dissolve on screen. Touch makes the poem urgent. Stillness makes it slow. The plant is the author.

---

## Quick Start

> **Before starting:** upload `rootstock_sensor.ino` to the Arduino via Arduino IDE.
> The Arduino must be connected and Arduino IDE's Serial Monitor must be **closed**.

Open **three terminal tabs** and run each command in order:

**Terminal 1 — Frontend server**
```bash
cd /Users/mac/rootstock && python -m http.server 8080
```

**Terminal 2 — Main backend**
```bash
cd /Users/mac/rootstock && python rootstock.py
```

**Terminal 3 — (Optional) Sensor debug monitor**
```bash
cd /Users/mac/rootstock && python tools/debug_vibration.py
```

Then open your browser:

| Page | URL |
|---|---|
| Visualization | `http://localhost:8080` |
| WebSocket monitor (terminal) | `ws://localhost:8765` |
| OSC output for Max/MSP | `localhost:9000` |

Touch or breathe on the plant — words appear when vibration presence rises above 0.02. Stronger contact → faster lines, stress-response vocabulary.

### Running without the Arduino (DUMMY mode)

To rehearse or develop the piece without the sensor, set `DUMMY = True` at the top of `rootstock.py`, or pass the flag:

```bash
python rootstock.py --dummy
```

No serial port is opened. Instead, a simulator emits raw values shaped like the Arduino's output (25 Hz, 0–1023), and they go through the same noise floor, normalization and smoothing as real sensor data. It moves on its own between the three vibration levels, staying in each one for a random time:

| State | Simulated raw | Resulting presence | Effect |
|---|---|---|---|
| `still` | ~300 | ≈ 0.00 | below the 0.02 gate: silence |
| `touch` | ~430, flickering | ≈ 0.2 – 0.6 | circadian and stress words mixed, ~10 s cycles |
| `intense` | ~760, wide spread | ≈ 0.65 – 1.0 | mostly stress-response words, fast short lines |

To force a state from the backend terminal, type a key and then Enter: `0` still · `1` touch · `2` intense · `a` back to auto. Enter on its own still pauses or resumes. The states are defined in `DUMMY_STATES` in `rootstock.py`.

---

## Visual Logic

```
╔══════════════════════════════════════════════════════════════════════╗
║  PLANT  →  SENSOR  →  SIGNAL  →  MODEL  →  LANGUAGE  →  SCREEN     ║
╚══════════════════════════════════════════════════════════════════════╝

  Arabidopsis thaliana (CCA1 gene)
        │
        │  piezoelectric vibration sensor
        ▼
  [Arduino Uno]
    baseline subtraction → ×100 amplify → Serial 9600 baud (25 Hz)
        │
        ▼
  [vibration_thread]  ── Python ──────────────────────────────────────
    raw ADC value
      − NOISE_FLOOR (320)          ← dead zone, rejects idle noise
      ÷ recent_max (adaptive)      ← normalization to [0, 1]
      → presence score             ← exponential smoothing (0.0 – 1.0)
        │
        │   presence gates vocabulary AND generation speed:
        │
        ├── 0.0 (still)  →  circadian dominates   (sleep, return, night…)   20 s/cycle
        ├── 0.5 (touch)  →  mixed, photosynthesis (light, leaf, absorb…)    10 s/cycle
        │                   as a constant bridge
        └── 1.0 (intense)→  stress_response       (threshold, rupture…)      1 s/cycle
        │
        ▼
  [HyenaDNA]  LongSafari/hyenadna-tiny-1k-seqlen-hf
    SEED: "ATGGATCTCGAGAAGAGAAGAGTTTCAGAG"  (CCA1 first 30 bp)
      → generates 30 new nucleotides per cycle
      → appended to growing sequence (last 512 bp as context)
        │
        ▼
  [sequence_to_words]
    split into codons (3-bp windows)
      → lookup in codon_word_mapping.json (one word per gene function)
      → choose a reading:  gene_fn drawn ∝ FUNCTION_WEIGHTS[gene_fn](presence)
      → presence-weighted gate:  random() < FUNCTION_WEIGHTS[gene_fn](presence)
      → words either surface or remain silent
        │
        ▼
  [WebSocket :8765]  ─────────────────────────────────────────────────
        │                              │
        ▼                              ▼
  [index.html]                   [OSC :9000]
  DNA helix visualization        Max/MSP or any OSC receiver
  plankton float animation       /rootstock/line, /rootstock/presence, …
  presence meter
  background video layer
```

### How presence chooses the words

`codon_word_mapping.json` gives each codon one word per gene function: all 64 codons have a circadian word and a stress-response word, and 41 also have a photosynthesis word. For example, `ATG` can be read as a circadian word or as `depth` (stress response).

Every time a codon is translated, `sequence_to_words` does two things:

1. **Chooses a reading.** It picks one of the codon's gene functions at random, weighted by `FUNCTION_WEIGHTS` at the current presence. Circadian weighs more when the plant is still, stress response more when it is touched, and photosynthesis stays constant.
2. **Gates it.** The chosen word then surfaces only with that same probability; otherwise the codon stays silent.

So the same DNA reads differently depending on the plant's state. With the current mapping and weights, the words that surface split like this:

| Presence | Circadian | Photosynthesis | Stress response | Codons that surface |
|---|---|---|---|---|
| 0.0 (still) | 86% | 14% | 0% | 87% |
| 0.45 (touch) | 53% | 21% | 26% | 56% |
| 1.0 (intense) | 4% | 15% | 81% | 76% |

(Below presence 0.02 no generation happens at all; the "still" row describes the vocabulary just above that gate.) The dominant function of each line is sent as `function` in the `line` WebSocket message, and each word's own function is sent in its `codon` message.

---

## File Structure

```
rootstock/
├── rootstock.py              # Main backend: sensor → model → WebSocket
├── rootstock_sensor.ino      # Arduino sketch: ADC read, amplify, Serial output
├── index.html                # Browser visualization: helix, plankton, presence meter
├── tools/
│   ├── build_mapping.py      # Offline: fetch NCBI genes → build root-level mapping
│   └── debug_vibration.py    # Diagnostic: live sensor terminal dashboard
├── codon_word_mapping.json   # Pre-built codon → word table (commit this, don't re-run)
└── bg.mp4                    # Background video (add your own, not tracked by git)
```

---

## Requirements

```bash
pip install -r requirements.txt
```

| Dependency | Purpose |
|---|---|
| `pyserial` | Arduino serial communication |
| `websockets` | Real-time browser updates |
| `python-osc` | OSC output for Max/MSP |
| `transformers` | HyenaDNA model loading |
| `torch` | Neural network inference |
| `biopython` | NCBI gene fetch (`tools/build_mapping.py` only) |
| `sentence-transformers` | Semantic word embeddings (`tools/build_mapping.py` only) |
| `wordfreq` | Word frequency ranking (`tools/build_mapping.py` only) |

---

## WebSocket Interface

The backend broadcasts JSON messages on `ws://localhost:8765`. You can monitor or drive the installation from any terminal without opening a browser.

### Connect with wscat

```bash
# Install wscat (requires Node.js)
npm install -g wscat

# Connect — messages stream as the plant generates poetry
wscat -c ws://localhost:8765
```

### Connect with websocat

```bash
# Install via Homebrew (no Node.js required)
brew install websocat

# Connect
websocat ws://localhost:8765
```

### Connect with Python

```python
import asyncio, websockets, json

async def monitor():
    async with websockets.connect("ws://localhost:8765") as ws:
        async for raw in ws:
            msg = json.loads(raw)
            print(msg)

asyncio.run(monitor())
```

### Message types

| `type` | Fields | Description |
|---|---|---|
| `line` | `text`, `codons`, `function` | A new poetry line. `function` is one of `circadian`, `photosynthesis`, `stress_response`. |
| `codon` | `codon`, `word`, `function` | Single codon translated to a word. |
| `presence` | `value` | Vibration presence score (0.0 – 1.0), sent at 5 Hz. |
| `dna` | `sequence` | Current DNA sequence after each HyenaDNA extension. |
| `cycle` | `count` | Generation cycle number. |
| `status` | `message` | Backend status events (startup, sensor connect/disconnect). |

### Example session

```
$ wscat -c ws://localhost:8765
Connected (press CTRL+C to quit)
< {"type":"status","message":"vibration sensor connected"}
< {"type":"presence","value":0.031}
< {"type":"presence","value":0.187}
< {"type":"codon","codon":"ATG","word":"root","function":"circadian"}
< {"type":"codon","codon":"GAT","word":"return","function":"circadian"}
< {"type":"line","text":"root return","codons":[...],"function":"circadian"}
< {"type":"dna","sequence":"ATGGATCTC..."}
< {"type":"presence","value":0.412}
```

### Pretty-print with jq

```bash
# Show only poetry lines with their gene function
wscat -c ws://localhost:8765 | grep --line-buffered '"type":"line"' | jq '.text, .function'

# Watch presence score only
wscat -c ws://localhost:8765 | grep --line-buffered '"type":"presence"' | jq '.value'
```

---

## Troubleshooting

### `✗ Arduino not found — retrying in 2 s...`
- Check USB cable is plugged in
- Close Arduino IDE's Serial Monitor (it holds the port exclusively)
- Run `python tools/debug_vibration.py` to see which ports are detected
- On macOS, check `ls /dev/cu.*` for available devices

### `✗ WebSocket failed (port 8765 in use?)`
```bash
kill $(lsof -ti :8765)
python rootstock.py
```

### `vibration 0.000` — sensor not responding
- Run `python tools/debug_vibration.py` (stop rootstock.py first) and watch the **Raw** column
- Idle raw value should be ~300 (100× amp of ~3 ADC units of sensor noise)
- If raw stays 0, the Arduino sketch may not be uploaded correctly
- If raw is very high (> 500) at rest, `NOISE_FLOOR` in rootstock.py needs to be raised

### Words generate without touching the plant
- `NOISE_FLOOR` is too low — environmental vibration leaks through
- Increase `NOISE_FLOOR` in rootstock.py (try 340–380)
- Check that Arduino IDE is closed and no other process reads the port

### Words never appear even when touching
- `NOISE_FLOOR` is too high — signal is being fully filtered
- Lower `NOISE_FLOOR` (try 280–300)
- Or lower the gate threshold: `if presence < 0.02` → `if presence < 0.01`

### HyenaDNA loads slowly / crashes on M1
- Normal: first load downloads ~400 MB of model weights
- Subsequent runs use local cache (`~/.cache/huggingface/`)
- If it crashes with memory error, close other GPU-heavy applications

### `codon_word_mapping.json` not found
```bash
# Regenerate the mapping table (requires biopython, sentence-transformers)
# Edit tools/build_mapping.py first: set Entrez.email = "your@email.com"
python tools/build_mapping.py
```

---

## References

- **HyenaDNA**: Nguyen et al., *"HyenaDNA: Long-Range Genomic Sequence Modeling at Single Nucleotide Resolution,"* NeurIPS 2023. [arxiv.org/abs/2306.15794](https://arxiv.org/abs/2306.15794)
- **Sentence-BERT**: Reimers & Gurevych, *"Sentence-BERT,"* EMNLP 2019. [arxiv.org/abs/1908.10084](https://arxiv.org/abs/1908.10084)
- **Gene sources**: NCBI Nucleotide — CCA1 `NM_001035612`, Photosynthesis `AY091856`, Stress response `NM_124370`
- **wordfreq**: Robyn Speer et al. [github.com/rspeer/wordfreq](https://github.com/rspeer/wordfreq)

---

*THE ROOTSTOCK · Arabidopsis thaliana · CCA1 · © Yvonne Wang*
