---
title: "Mochi SFT on the published protein-sequence probes"
source_repository: "marin-community/marin"
task_issue: "https://github.com/marin-community/marin-logbooks/issues/73"
related_issues: ["https://github.com/marin-community/marin/issues/11723"]
research_branch: "codex/73-mochi-bioseq"
marin_base_revision: "4aab2cff9e0a7d105928ba0bdc575ee5c43c6610"
archive_template_revision: "fc1cb30da0967f6ab78ffd47444ed9144f50cba9"
status: active
---

# Mochi SFT on the published protein-sequence probes

## Question

Where does Mochi SFT step 400407 fall relative to the language models in
[Tu and Chen's study](https://xinmingtu.cn/blog/2026/claude-opus-5-5-protein-sequences/)?
The primary outcome is next-amino-acid accuracy on real proteins minus accuracy
on their shuffled controls, averaged across 167 protein families. Bootstrap
confidence intervals resample families, not individual prompts.

## Baseline and plan

Use [LLM-bioseq-probe at 1c0eb30](https://github.com/XinmingTu/LLM-bioseq-probe/tree/1c0eb30ec4aacc672386ec5e8a0d7e9f17e8e649)
for the exact released inputs, prompt constructors, shuffles, parsers and
published comparator outputs. Run instructed chat and untemplated FASTA
continuation on the main set. Supplement with the smaller-family set,
August-2026 proteins, pilot set, 158 published recall cases, DNA windows and
contact maps. The August set is not established to be after Mochi's cutoff.

The existing shared server reports `hero-sft-step400407`, max length 262144.
Its reported checkpoint is the export ending in
`hero-sft-fullab-step393k-v2/step-400407/2af3f828aa-cpu-a1`.
The endpoint is backed by an existing four-task serving job; this experiment
allocates no additional GPUs. Server identity is checked before and after
each suite, but this remains a shared serve rather than an immutable deployment.

Prompts use `marin-community/marin-tokenizer-tools` at
`3f6879acdaf3dae7874b2d48c7d11b9d8c6faf09`, rendered locally with
`enable_thinking=False`. Every request uses temperature 0,
`stop_token_ids=[128009,128001]`, and `skip_special_tokens=false`.
Token limits follow the authors' scripts: 4000 for chat protein/DNA,
24 for raw continuation, 8000 for recall, 32000 for contact maps.
The only text normalization maps Marin's thinking delimiters to the tags the
upstream parser recognizes; original response text is retained alongside it.

Transport failures get at most three attempts, each recorded. Non-answers are
never retried. Preserve all outputs, and distinguish errors and exclusions from
incorrect predictions. Report the main score under the published exclusion rule
and with non-answers counted as incorrect.

## Reproduction

Clone the upstream repository at the revision above. Download the pinned
tokenizer with `hf download`. Supply an endpoint URL ending in `/v1` in a local
credential file outside the repository. Set `PYTHONPATH` to the upstream `src`
directory, then run `run.py` with these arguments:

```bash
uv run --no-project --with google-genai --with anthropic --with openai \
  --with python-dotenv --with pandas --with numpy --with httpx \
  --with transformers --with biopython --with gemmi --with scipy --with matplotlib \
  python logbooks/73-mochi-bioseq/run.py \
  --upstream /path/to/LLM-bioseq-probe --endpoint-file /private/endpoint.txt \
  --tokenizer /path/to/pinned-tokenizer --output /path/to/results \
  --suite confirm --workers 8 --run
```

Omit `--run` for input validation without inference. Other suite names are
`raw`, `small`, `d_set`, `pilot`, `recall`, `dna`, and `contact`.
Keep each output directory for one source revision and request configuration.
The runner resumes existing question keys without overwriting records.

## Entries

### 2026-11-23 — Setup

The live model identity matches prior Mochi SFT evaluations in issue 11723.
Published comparator scores will be distinguished from newly measured Mochi
scores. No model result is available yet. The next step is a small prompt and
parser smoke check before the complete suite.
