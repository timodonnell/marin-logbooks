# Copyright The Marin Authors
# SPDX-License-Identifier: Apache-2.0

"""Run the pinned bioseq probes against an existing Mochi completion endpoint."""

import argparse
import csv
import gzip
import hashlib
import importlib.metadata
import json
import logging
import os
import runpy
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import httpx
from bioseq_probe import dna, probe
from transformers import AutoTokenizer

UPSTREAM_REVISION = "1c0eb30ec4aacc672386ec5e8a0d7e9f17e8e649"
MODEL = "hero-sft-step400407"
TOKENIZER_REVISION = "3f6879acdaf3dae7874b2d48c7d11b9d8c6faf09"
PROTEIN_SETS = {
    "confirm": "confirm_sequences.tsv",
    "raw": "confirm_sequences.tsv",
    "small": "small_family_sequences.tsv",
    "d_set": "d_set_sequences.tsv",
    "pilot": "sequences_pilot.tsv",
}
CONTACT_SUITES = ("contact", "contact_list")
LOG = logging.getLogger(__name__)


@dataclass(frozen=True)
class Question:
    key: str
    prefix: str
    prompt: str
    target: str
    metadata: dict
    max_tokens: int


def table(path: Path) -> list[dict]:
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def questions(upstream: Path, suite: str, output: Path) -> list[Question]:
    result = []
    if suite in PROTEIN_SETS:
        for row in table(upstream / "data" / PROTEIN_SETS[suite]):
            for cut in map(int, row["cuts"].split(",")):
                for condition, sequence in (
                    ("real", row["sequence"]),
                    ("shuffled", probe.shuffled(row["sequence"], row["accession"])),
                ):
                    prefix = sequence[:cut]
                    text = probe.prompt(prefix)
                    if suite == "raw":
                        text = ">sequence\n" + probe.fasta(prefix)
                        if len(prefix) % 60 == 0:
                            text += "\n"
                    result.append(
                        Question(
                            f"{row['accession']}:{cut}:{condition}",
                            prefix,
                            text,
                            sequence[cut : cut + probe.K],
                            {"accession": row["accession"], "pfam": row["pfam"], "cut": cut, "condition": condition},
                            24 if suite == "raw" else 4000,
                        )
                    )
        return result
    if suite == "dna":
        for row in table(upstream / "data/dna_confirm_windows.tsv"):
            for condition, sequence in (
                ("real", row["window"]),
                ("shuffled", dna.dinucleotide_shuffle(row["window"], row["locus"])),
            ):
                prefix = sequence[: dna.PREFIX]
                result.append(
                    Question(
                        f"{row['locus']}:{condition}",
                        prefix,
                        dna.prompt(prefix),
                        sequence[dna.PREFIX : dna.PREFIX + dna.K],
                        {
                            "locus": row["locus"],
                            "species": row["species"],
                            "region": row["region"],
                            "condition": condition,
                        },
                        4000,
                    )
                )
        return result
    if suite == "recall":
        sequences = {r["accession"]: r["sequence"] for r in table(upstream / "data/sequences_pilot.tsv")}
        # The published common set fixes selection independently of Mochi's answers.
        with gzip.open(upstream / "results/outputs/recall_length/claude-opus-5-5.jsonl.gz", "rt") as stream:
            cases = [json.loads(line)["id"] for line in stream]
        assert len(cases) == len(set(cases)) == 158
        for key in cases:
            accession, cut_text = key.split(":")
            cut = int(cut_text)
            prefix = sequences[accession][:cut]
            text = probe.prompt(prefix).replace(
                "exactly the next 10 amino acids", "the next 200 amino acids (or fewer if the protein ends)"
            )
            result.append(
                Question(
                    key, prefix, text, sequences[accession][cut:], {"id": key, "accession": accession, "cut": cut}, 8000
                )
            )
        return result
    assert suite in CONTACT_SUITES
    contact_prompt = "list" if suite == "contact_list" else "grid"
    # The upstream loader verifies each sequence against the structure used for truth.
    contact = runpy.run_path(str(upstream / "scripts/probe/contact_map.py"))
    truth_path = output / "contact_truth.json"
    if truth_path.exists():
        truth = json.loads(truth_path.read_text())
    else:
        original = Path.cwd()
        try:
            os.chdir(upstream)
            truth = contact["load_truth"]()
        finally:
            os.chdir(original)
        truth_path.write_text(json.dumps(truth))
    for name, row in truth.items():
        sequence = row["sequence"]
        result.append(
            Question(
                name,
                sequence,
                contact["PROMPTS"][contact_prompt](sequence),
                "",
                {"protein": name, "pdb": row["pdb"], "length": len(sequence), "prompt": contact_prompt},
                32000,
            )
        )
    return result


def answer(client: httpx.Client, base_url: str, tokenizer, suite: str, q: Question) -> dict:
    if suite == "raw":
        rendered = q.prompt
    else:
        rendered = tokenizer.apply_chat_template(
            [{"role": "user", "content": q.prompt}],
            tokenize=False,
            add_generation_prompt=True,
            enable_thinking=False,
        )
    payload = {
        "model": MODEL,
        "prompt": rendered,
        "temperature": 0.0,
        "max_tokens": q.max_tokens,
        "stop_token_ids": [128009, 128001],
        "skip_special_tokens": False,
        "add_special_tokens": suite == "raw",
    }
    rec = {
        "model": "mochi-sft",
        "suite": suite,
        "key": q.key,
        **q.metadata,
        "target": q.target,
        "prompt_sha256": hashlib.sha256(rendered.encode()).hexdigest(),
    }
    started = time.monotonic()
    # Retry only transport errors, never a model non-answer. Preserve all attempts.
    attempts = []
    for attempt in range(3):
        try:
            response = client.post(base_url + "/completions", json=payload)
        except httpx.TransportError as exc:
            attempts.append({"error_type": type(exc).__name__})
        else:
            attempts.append({"status": response.status_code})
            if response.status_code == 200:
                data = response.json()
                choice = data["choices"][0]
                raw_text = choice["text"]
                text = raw_text.replace("<|start_think|>", "<think>").replace("<|end_think|>", "</think>")
                rec.update(raw_text=raw_text, text=text, stop=choice.get("finish_reason"), usage=data.get("usage", {}))
                break
            if response.status_code not in (429, 500, 502, 503, 504):
                raise RuntimeError(f"Non-retryable completion status {response.status_code}")
        if attempt < 2:
            time.sleep(2**attempt)
    else:
        rec["error"] = "transport retries exhausted"
    rec.update(attempts=attempts, elapsed=time.monotonic() - started)
    if "error" in rec:
        return rec
    text, stopped = rec["text"], probe.cut_off(rec)
    if suite in CONTACT_SUITES:
        return rec
    if suite == "raw":
        pred = probe.parse_raw(text)
    elif suite == "dna":
        pred = dna.parse(text, q.prefix, cut=stopped)
    else:
        pred = probe.parse(text, q.prefix, prefer_len=20 if suite == "recall" else 1, cut=stopped)
    rec["pred"] = pred
    if suite not in ("dna", "recall"):
        rec.update(probe.score(pred, q.target))
    if suite == "recall":
        run = 0
        while run < min(len(pred), len(q.target)) and pred[run] == q.target[run]:
            run += 1
        rec.update(exact_run=run, written=len(pred), truth_available=len(q.target))
    return rec


def identity(models: dict) -> list[tuple]:
    return [(m["id"], m["root"], m["max_model_len"]) for m in models["data"]]


def models_snapshot(client: httpx.Client, base_url: str) -> dict:
    # A capability URL is a credential; never include HTTP exceptions containing it.
    response = client.get(base_url + "/models")
    if response.status_code != 200:
        raise RuntimeError(f"Model discovery returned HTTP {response.status_code}")
    return response.json()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--upstream", type=Path, required=True)
    ap.add_argument("--endpoint-file", type=Path, required=True)
    ap.add_argument("--tokenizer", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--suite", choices=[*PROTEIN_SETS, "dna", "recall", *CONTACT_SUITES], required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int)
    ap.add_argument("--run", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    logging.getLogger("httpx").setLevel(logging.WARNING)
    args.upstream = args.upstream.resolve()
    revision = subprocess.check_output(["git", "-C", str(args.upstream), "rev-parse", "HEAD"], text=True).strip()
    assert revision == UPSTREAM_REVISION, revision
    args.output.mkdir(parents=True, exist_ok=True)
    manifest = {
        "upstream_revision": revision,
        "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "tokenizer_revision": TOKENIZER_REVISION,
        "temperature": 0.0,
        "enable_thinking": False,
        "add_special_tokens": args.suite == "raw",
        "workers": args.workers,
        "inputs_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((args.upstream / "data").glob("*.tsv"))
        },
        "tokenizer_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(args.tokenizer.glob("*")) if p.is_file()
        },
        "packages": {name: importlib.metadata.version(name) for name in ("httpx", "transformers", "numpy", "pandas")},
    }
    manifest_path = args.output / f"{args.suite}-manifest.json"
    if manifest_path.exists():
        assert json.loads(manifest_path.read_text()) == manifest, "resume configuration changed"
    else:
        manifest_path.write_text(json.dumps(manifest, indent=2))
    items = questions(args.upstream, args.suite, args.output)[: args.limit]
    assert len(items) == len({q.key for q in items})
    LOG.info("%s: %d questions; run=%s", args.suite, len(items), args.run)
    if not args.run:
        return
    base_url = args.endpoint_file.read_text().strip().rstrip("/")
    tokenizer = AutoTokenizer.from_pretrained(args.tokenizer)
    out = args.output / f"{args.suite}.jsonl"
    previous = [json.loads(line) for line in out.read_text().splitlines()] if out.exists() else []
    assert len(previous) == len({r["key"] for r in previous}), "duplicate records"
    done = {r["key"] for r in previous}
    todo = [q for q in items if q.key not in done]
    request_timeout = 1800 if args.suite in CONTACT_SUITES else 600
    with httpx.Client(timeout=request_timeout) as client:
        before = models_snapshot(client, base_url)
        assert before["data"][0]["id"] == MODEL
        (args.output / f"{args.suite}-models-start.json").write_text(json.dumps(before, indent=2))
        with out.open("a") as stream, ThreadPoolExecutor(args.workers) as pool:
            futures = [pool.submit(answer, client, base_url, tokenizer, args.suite, q) for q in todo]
            errors = 0
            for i, future in enumerate(as_completed(futures), 1):
                record = future.result()
                stream.write(json.dumps(record) + "\n")
                stream.flush()
                errors += "error" in record
                if i % 100 == 0 or i == len(todo):
                    LOG.info("%s: %d/%d new replies, %d transport failures", args.suite, i, len(todo), errors)
        after = models_snapshot(client, base_url)
        (args.output / f"{args.suite}-models-end.json").write_text(json.dumps(after, indent=2))
        assert identity(before) == identity(after), "server identity changed"


if __name__ == "__main__":
    main()
