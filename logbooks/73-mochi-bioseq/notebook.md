---
title: "Mochi SFT on the published protein-sequence probes"
source_repository: "marin-community/marin"
task_issue: "https://github.com/marin-community/marin-logbooks/issues/73"
related_issues: ["https://github.com/marin-community/marin/issues/11723"]
research_branch: "codex/73-mochi-bioseq"
marin_base_revision: "4aab2cff9e0a7d105928ba0bdc575ee5c43c6610"
archive_template_revision: "fc1cb30da0967f6ab78ffd47444ed9144f50cba9"
status: complete
---

# Mochi SFT on the published protein-sequence probes

Mochi SFT step 400407 shows a small protein-continuation signal. On the main
668-protein benchmark, the real-minus-shuffled next-residue accuracy gap is
**+1.31 percentage points in chat** (95% family-bootstrap CI −0.15 to +2.78)
and **+3.28 points in raw FASTA continuation** (+1.89 to +4.65). This places
Mochi near the study's lower-scoring general language models, well below the
published Opus 5.5, Fable 5.2 and strongest protein-language-model results.
Neither Mochi format produces an exact 10-residue continuation on the main set.
A separate, selected recall probe finds four cases with at least 50 exact
residues (longest 266), compared with 46 for Opus 5.5; three Mochi cases share
the same or overlapping protein tail.

| Model and format | Gap, percentage points | 95% CI | Source |
| --- | ---: | --- | --- |
| Claude Fable 5.2, chat | +45.09 | +42.51 to +47.63 | Published |
| ProGen3-3B, prefix | +32.54 | +29.31 to +35.69 | Published |
| ProtGPT3-10B, prefix | +27.66 | +24.19 to +31.14 | Published |
| Claude Opus 5.5, chat | +27.41 | +24.79 to +30.06 | Published |
| ESM C 600M, prefix | +12.96 | +10.72 to +15.24 | Published |
| Mochi SFT, raw FASTA | +3.28 | +1.89 to +4.65 | This run |
| GPT-6 Astra, chat | +2.85 | +1.45 to +4.24 | Published |
| Mochi SFT, chat | +1.31 | −0.15 to +2.78 | This run |
| Gemini 4 Carbon, chat | +0.92 | −0.37 to +2.25 | Published |

The metric averages each family's real-versus-shuffled accuracy difference;
it is not absolute sequence accuracy. Other models were not rerun. Their
scores come from the pinned upstream result tables. Prompting formats differ,
as labeled, and overlapping intervals do not establish a ranking.

All nine suites are complete: 29,002 scored question cases. DNA produces no
exact 20-base continuations; contact pair lists perform near their random
baseline, and the grid format yields no usable long-range contact predictions.

![Main benchmark comparison; orange points are newly measured Mochi results.](results/comparison.png)

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

Pfam groups related protein domains by sequence similarity. A domain annotation
marks the segment of a full protein assigned to that family; the shown prefix
can include sequence before that segment.

The main set has four proteins from each of 167 Pfam families and five fixed
cuts per protein, inside its annotated domain. Each cut is queried on both the
real sequence and a composition-preserving shuffle of the full protein:
668 × 5 × 2 = 6680 requests per format. Chat asks for the next ten amino acids;
the headline accuracy scores only the first predicted residue. Exact-10 counts
require all ten residues to match and the target ten-residue sequence to be
absent from the prefix.

The primary rule excludes non-answers independently in each condition. A
non-answer is a refusal or a response from which the upstream parser extracts
no continuation. Compute the first-residue accuracy within each family and
condition, then average the real-minus-shuffled differences equally over
families. All 167 families remain in both main estimates. The sensitivity
analyses either count non-answers as incorrect or retain only cuts with valid
answers in both conditions.

Raw parsing takes up to ten residues from the beginning of the continuation,
stopping at a header, markup, a blank line or prose, following upstream
`parse_raw`. A nonempty answer shorter than ten residues still contributes to
first-residue accuracy but cannot count as an exact-10 continuation.

Mochi is the 535B-parameter model in Marin's
[post-training program](https://github.com/marin-community/marin/issues/11535);
SFT means supervised fine-tuning. This evaluation uses its step-400407 SFT
checkpoint, also evaluated in
[the prior weakness analysis](https://github.com/marin-community/marin/issues/11723).
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

Score recomputation from the public response bundle requires no Marin access.
Fresh inference requires access to Marin's serving infrastructure and the
step-400407 checkpoint; this record does not publish the model weights.

Clone the upstream repository at the revision above. Download the pinned
tokenizer with `hf download`. Supply an endpoint URL ending in `/v1` in a local
credential file outside the repository. Set `PYTHONPATH` to the upstream `src`
directory, then run `run.py` with these arguments:

```bash
uv run --no-project --with google-genai --with anthropic --with openai \
  --with python-dotenv --with pandas --with numpy --with httpx \
  --with transformers --with jinja2 --with biopython --with gemmi --with scipy --with matplotlib \
  python logbooks/73-mochi-bioseq/run.py \
  --upstream /path/to/LLM-bioseq-probe --endpoint-file /private/endpoint.txt \
  --tokenizer /path/to/pinned-tokenizer --output /path/to/results \
  --suite confirm --workers 32 --run
```

Omit `--run` for input validation without inference. Other suite names are
`raw`, `small`, `d_set`, `pilot`, `recall`, `dna`, `contact`, and `contact_list`.
Keep each suite's output file and manifest together. Resuming requires the same
runner and request settings; a directory containing several suites retains each
suite's own provenance. The runner resumes existing question keys without
overwriting records. Raw continuation used eight workers; other suites used 32.

## Entries

### 2026-11-23 — Setup

The live model identity matches prior Mochi SFT evaluations in issue 11723.
Published comparator scores will be distinguished from newly measured Mochi
scores. No model result is available yet. The next step is a small prompt and
parser smoke check before the complete suite.

Publication note: the current account has read access but no push access to the
shared logbooks repository. Preserve the research branch in its GitHub fork and
link it from the coordinating issue.

### 2026-11-23 — Tokenization correction before scoring

The initial public proxy intermittently returned HTTP 503; a standard authenticated
Iris tunnel reaches the same serving endpoint. All 20 smoke requests then completed,
six reaching the 4000-token limit with repeated amino-acid text.

A token-count audit found that the rendered prompt's 179 tokens became 180 at
the completion API: the tokenizer prepended a second BOS token. Stop the
incomplete `confirm` and `small` runs from runner revision `43f074efe4`, retain
their outputs as invalidated setup evidence, and exclude them from all scores.
Set `add_special_tokens=false` for rendered chat prompts and verify the API
prompt count before restarting into a new output directory. Raw FASTA
completion uses the API's single automatic BOS and is unaffected; its original
run can be retained.

The analysis reproduces the published Opus 5.5 real accuracy, shuffled accuracy,
and family-averaged gap exactly (34.8607367%, 7.4617920%, +27.4062402 points).
Bootstrap bounds vary slightly with random-number ordering; the estimator and
10,000-family-resample design match the source. Permutation p-values use the
finite-simulation correction `(hits + 1) / (10000 + 1)`.

### 2026-11-23 — Protein suites complete

The corrected chat prompts match the server's token IDs exactly on three
audited examples (179, 85 and 83 tokens). All five protein suites completed
without an exhausted transport retry. Chat, small-family, newer-protein and
pilot requests used 32 concurrent workers; raw continuation used eight.

Pooled main chat accuracy is 10.06% on real prefixes and 8.78% on shuffled
prefixes among parsed answers. Their difference differs slightly from the
headline gap because the latter weights families equally. There are 185
excluded non-answers: 10 real and 175
shuffled. Counting them as incorrect increases the family-averaged gap to
+1.71 points (95% CI +0.30 to +3.14); retaining only cuts with both answers
valid gives +1.31 (−0.17 to +2.83). The main chat run has 286 token-limited
responses out of 6680. Raw continuation has 42 non-answers; counting them as
incorrect changes its gap to +3.26 (+1.89 to +4.61). All raw responses reach
the prescribed 24-token cap, which is expected for untemplated continuation.

Neither format passes the authors' H1 criterion: one-sided family sign-flip
`p < 0.001` and the 95% interval's lower bound above 10 percentage points.
Neither passes H2 (positive Spearman correlation between the gap and log10
Pfam family size on the main 167 families, one-sided `p < 0.01`) or H3 (at least 10 exact
continuations, over 100 times the independent-position expectation, with none
on shuffled controls). These are the published decision rules applied to Mochi,
not newly preregistered hypotheses. For H3, estimate each of the ten positions'
accuracy within a family, multiply those ten probabilities, multiply by the
number of valid real-prefix answers in that family, then sum over families.
This is the expected exact-10 count if positions were independent.
Family size is the upstream table's
`family_proteins` count of members in the biological family, not the four
proteins sampled into the benchmark.

The smaller-family set has 90 families containing 20–300 known proteins each.
The August set contains proteins first public in RefSeq on August 2–27, 2026,
whose closest pre-August relatives are 34–49% identical according to the
authors' selection. In the pilot, “well-studied” means a named family with
at least 100 PDB structures; “obscure” means a domain-of-unknown-function
family with no PDB structure and at most three Europe PMC hits for its Pfam
accession. These are the authors' proxies for literature attention.

Each row below uses the main benchmark's equal family weighting, independent
non-answer exclusions and family bootstrap. Family counts are the independent
units resampled for the confidence intervals.

| Supplementary chat set | Families | Questions | Gap, percentage points | 95% CI |
| --- | ---: | ---: | ---: | --- |
| Smaller families | 90 | 3600 | +0.57 | −1.46 to +2.61 |
| August-2026 proteins | 31 | 560 | −7.22 | −13.57 to −1.12 |
| Main proteins restricted to the same 31 families | 31 | 1240 | +2.02 | −1.07 to +5.08 |
| Well-studied pilot families | 19 | 1140 | +8.59 | +4.06 to +13.19 |
| Obscure pilot families | 18 | 1080 | −1.59 | −4.38 to +1.08 |

Across the main and smaller-family sets (257 families), the correlation of
the accuracy gap with log family size is Spearman ρ = 0.050, two-sided
`p = 0.421`. The pilot's stronger performance on well-studied families is
descriptive: these are different families, so this comparison does not isolate
the causal effect of literature exposure. The newer-protein set is small
(56 proteins, 31 families), and Mochi's training cutoff has not been verified.
These probes do not establish which sequences were in its training data.

The analysis reproduces the upstream Opus 5.5 gap, all 12 exact continuations,
and the independent-position expectation of 0.5706195603 exactly. The raw
runner is preserved at `43f074efe4`; the corrected chat runner is at
`26493a821c`. Contact-pair-list support is added at `80ebf04db2` without changing
protein or DNA requests. Per-suite manifests identify the executed runner by
SHA-256, tokenizer revision, input hashes and package versions.


### 2026-11-23 — DNA and long recall complete

The DNA probe uses 4500 windows: 500 from each of nine species/region cells,
with a real and a dinucleotide-preserving shuffled version of each window.
The prompt shows 200 bases and asks for the next 20. As upstream, score the
first base and estimate the paired real-minus-shuffled gap only where both
answers parse; bootstrap windows within each cell. There are 9000 responses,
including 190 that reached the 4000-token cap.

| DNA cell | Paired windows | Gap, percentage points | 95% CI |
| --- | ---: | ---: | --- |
| Bacteria/archaea, coding | 452 | −0.22 | −5.97 to +5.53 |
| Human, coding | 457 | +6.78 | +1.09 to +12.47 |
| Human, intergenic | 460 | +1.52 | −4.57 to +7.61 |
| Human, intron | 459 | −3.05 | −8.71 to +2.61 |
| Human, promoter | 452 | −3.10 | −8.63 to +2.22 |
| Mouse, coding | 457 | +4.16 | −1.75 to +10.07 |
| Mouse, intergenic | 459 | −1.09 | −6.75 to +4.79 |
| Mouse, intron | 456 | 0.00 | −5.92 to +5.92 |
| Mouse, promoter | 462 | −0.43 | −6.28 to +5.19 |

None passes the published D1 rule (one-sided sign-flip `p < 0.001` and lower
95% confidence bound above five percentage points), and there are no exact
20-base continuations in either condition, so none passes D2 (at least three
exact real-prefix continuations and none on shuffled controls). Human coding
DNA has a positive interval, but these nine intervals are not adjusted for
multiple comparisons. For comparison, the published Opus 5.5 human-coding gap
is +30.4 points (+24.4 to +36.4), with 25 exact 20-base continuations; its
mouse-coding gap is +28.66 (+23.05 to +34.27), with five exact continuations.

Long recall uses the same 158 pilot cut points reported for the comparator
models. The authors selected these because at least one tested model had
previously continued ten residues correctly; selection is fixed independently
of Mochi. Ask for 200 residues or the remainder if the protein ends, with an
8000-token limit. Count consecutive correct residues immediately after the
cut, following the original parser's removal of any echoed prefix suffix.
The score can exceed 200 if the model writes more than requested.

| Model | At least 20 residues | At least 50 | At least 100 | At least 200 | Longest |
| --- | ---: | ---: | ---: | ---: | ---: |
| Mochi SFT, this run | 4 | 4 | 1 | 1 | 266 |
| Claude Opus 5.5, published | 71 | 46 | 27 | 15 | 314 |
| Gemini 4 Carbon, published | 74 | 42 | 23 | 16 | 201 |
| GPT-5.5, published | 0 | 0 | 0 | 0 | 2 |

Mochi's four successes are not four independent proteins. Three cover the
same or overlapping tail in two PF00019 proteins: `K7C4V8` at cuts 330 and
336, and `A0A1S3ERL3` at cut 330. They reproduce 68, 62 and 68 residues,
respectively, through the end of the protein. The two 68-residue tails are
identical. The fourth is `A0A8C2LMQ0`, a Chinese hamster PF00339 protein,
at cut 101: 266 consecutive residues match before the first mismatch.
The response first echoes 41 residues from the shown prefix; those are removed
and are not counted in the 266. All four exact continuations occur literally
in the saved response and are absent from the shown prefix.

This is evidence of occasional long exact recall, consistent with sequence
memorization, despite weak average performance on the main benchmark. It does
not identify which training source supplied the sequence. The selected recall
set is enriched for cases other models could recall, so 4/158 is not an
estimate of the fraction of arbitrary proteins Mochi has memorized. Of 158
responses, 154 parse and 93 reach the 8000-token limit.

The contact-list client timeout was increased from 600 to 1800 seconds to
accommodate its 32000-token allowance. The 23 completed replies from runner
`80ebf04db2` were preserved and the remaining requests used `1b65240b2d`;
request payloads were unchanged. Separating the final grid run from the
sequential launcher also ended its recall child during terminal cleanup.
Recall resumed with the identical runner and settings, retaining all 77
completed replies. The archive records both resumptions and the retained
responses; canceled in-flight calls are not included among completed cases.


### 2026-11-23 — Contact pair lists complete

The structural probe asks for contacts from the sequence alone, using 52
proteins: six famous structures and 46 recent structures. A true contact has
C-beta atoms (C-alpha for glycine) less than 8 Å apart. The headline score
counts predicted contacts between residues separated by at least 24 positions
in the sequence. Pair-list parsing and scoring reuse the original functions:
ignore pairs outside the sequence or separated by fewer than six positions,
deduplicate pairs, and score the upper triangle. A parsed list is considered
valid even if some lines are rejected. The alternative grid format requires
exactly an L-by-L binary matrix for a protein of length L, with no padding or
cropping.

All 52 Mochi lists parse. Their long-range predictions are close to chance:

| Group | Correct / predicted contacts | Precision | Expected correct, random | Expected correct, separation-matched |
| --- | ---: | ---: | ---: | ---: |
| Six famous proteins | 93 / 1502 | 6.19% | 88.30 | 88.64 |
| 46 recent structures | 3075 / 114734 | 2.68% | 3133.29 | 3131.21 |
| Two author-labeled “clean” cases, a subset of recent structures | 115 / 6342 | 1.81% | 115.41 | 113.97 |

The random expectation preserves each protein's number of predicted contacts
and uses its true long-range contact density. The separation-matched expectation
also preserves each predicted pair's sequence separation. Mochi labels 114734
of 216244 possible long-range pairs as contacts on recent structures (53.1%),
compared with only 964 calls in the published Opus 5.5 lists. Opus gets 119 of
those correct (12.34%; random expectation 29.44), and Fable 5.2 gets 319/1775
correct (17.97%; random expectation 54.22). These methods choose different
numbers of pairs; precision alone is not a comparison at equal recall or an
equal prediction budget. Mochi's near-baseline performance gives no evidence
of useful structural prediction in this configuration.

Nine replies reach the 32000-token cap. Fourteen requests receive a gateway
HTTP 504 and succeed on the recorded second attempt; none exhausts retries.
The grid results are pending. The authors' novelty categories, including
“clean,” are retained for comparison but do not establish novelty relative to
Mochi's unverified training cutoff. Here “clean” means no earlier PDB relative
at an E-value below 0.001, with sequence first public on or after July 1, 2026
in the authors' UniProtKB, UniParc and organism-record checks.


### 2026-11-23 — Grid results and final audit

All 52 grid requests completed. Only four responses satisfy the strict matrix
format: one of six famous proteins and three of 46 recent structures. All four
contain zero long-range contact predictions, so their long-range precision is
undefined, not zero. The other 48 are excluded as malformed matrices. Forty
responses reached the 32000-token cap. This format does not yield a usable
long-range structural prediction from Mochi under the published protocol.

The completed run contains 29,002 question cases across nine suites, with no
exhausted transport retries. This count excludes setup, invalidated double-BOS
runs, canceled in-flight calls and transport retries. There are 30 recorded
retries among scored cases: 16 raw-continuation retries and 14 contact-list
HTTP 504 retries. The full audit checks unique question coverage, targets,
metadata, regenerated prompt hashes and parser outputs. Model identity
projections agree at every retained start and end snapshot. This checks the
server's reported identity, not an immutable deployment guarantee.

The main conclusion is limited to this checkpoint with thinking disabled and
the stated token budgets. Mochi shows a modest raw-continuation advantage over
shuffled sequences and occasional long exact recall, but it is far below the
published leading Claude and protein-language-model scores. The DNA and
contact probes do not add evidence of broad biological sequence or structural
prediction ability. Comparator models were not rerun; supplementary intervals
are not corrected for multiple comparisons, and Mochi's cutoff is unknown.


The compressed-response reproduction produces identical summaries for all nine
suites and identical score tables (ignoring trailing whitespace). Original
responses, contact truth, public identity projections, per-suite manifests and
setup/resumption evidence are preserved in the
[response bundle](https://github.com/timodonnell/marin-logbooks/releases/download/mochi-bioseq-73-2026-11-23/mochi-bioseq-step400407-2026-11-23.tar.gz).

Bundle SHA-256: `26c9a98a568fc34b6b1612f6a4d7743a4a9894c4d7ab47b2cb619c2aefa0446c`.
The bundle contains an internal `SHA256SUMS` file. Extract it and pass its
`results` directory to `analyze.py` with the same dependency command as above:

```bash
python logbooks/73-mochi-bioseq/analyze.py \
  --upstream /path/to/LLM-bioseq-probe \
  --output /path/to/mochi-bioseq-step400407/results
```

This reads compressed responses directly and requires no model endpoint.
The analysis source is preserved at
[77abdc0f92](https://github.com/timodonnell/marin-logbooks/tree/77abdc0f92cdbebcdbc40b9d39ba49198cb823a7/logbooks/73-mochi-bioseq).
