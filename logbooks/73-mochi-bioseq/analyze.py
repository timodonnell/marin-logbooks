# Copyright The Marin Authors
# SPDX-License-Identifier: Apache-2.0

"""Audit saved Mochi replies and compare with the pinned published results."""

import argparse
import ast
import gzip
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from bioseq_probe import dna, probe
from run import PROTEIN_SETS, questions
from scipy.stats import hypergeom, spearmanr


def records(path: Path) -> list[dict]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as stream:
        return [json.loads(line) for line in stream]


def interval(values: np.ndarray) -> dict:
    rng = np.random.default_rng(2026)
    boot = values[rng.integers(0, len(values), (10_000, len(values)))].mean(1)
    low, high = np.percentile(boot, [2.5, 97.5])
    flips = (rng.choice([-1, 1], (10_000, len(values))) * values).mean(1)
    return {
        "gap": float(values.mean() * 100),
        "low": float(low * 100),
        "high": float(high * 100),
        "p_one_sided": float((1 + (flips >= values.mean()).sum()) / 10_001),
    }


def protein_summary(frame: pd.DataFrame) -> dict:
    valid = frame[frame.valid]
    fam = valid.groupby(["pfam", "condition"]).hit.mean().unstack().dropna()
    all_fam = frame.groupby(["pfam", "condition"]).hit.mean().unstack().dropna()
    pair_keys = ["accession", "cut"]
    paired = frame.groupby(pair_keys).valid.transform("all")
    paired_fam = frame[paired].groupby(["pfam", "condition"]).hit.mean().unstack().dropna()
    return {
        "questions": len(frame),
        "valid": len(valid),
        "families": len(fam),
        "real": float(valid[valid.condition == "real"].hit.mean() * 100),
        "shuffled": float(valid[valid.condition == "shuffled"].hit.mean() * 100),
        **interval((fam.real - fam.shuffled).to_numpy()),
        "nonanswers_wrong": interval((all_fam.real - all_fam.shuffled).to_numpy()),
        "both_answers_valid": interval((paired_fam.real - paired_fam.shuffled).to_numpy()),
        "exact10_real": int(valid[(valid.condition == "real")].all10.sum()),
        "exact10_shuffled": int(valid[(valid.condition == "shuffled")].all10.sum()),
        "positive_families": int((fam.real > fam.shuffled).sum()),
        "nonanswers_real": int((~frame[frame.condition == "real"].valid).sum()),
        "nonanswers_shuffled": int((~frame[frame.condition == "shuffled"].valid).sum()),
    }


def audit_suite(upstream: Path, output: Path, suite: str) -> tuple[list[dict], list]:
    saved = records(output / f"{suite}.jsonl")
    expected = questions(upstream, suite, output)
    by_key = {q.key: q for q in expected}
    assert len(saved) == len(expected) == len({d["key"] for d in saved}), suite
    assert set(by_key) == {d["key"] for d in saved}, suite
    assert not any("error" in d for d in saved), f"transport errors in {suite}"
    for d in saved:
        q = by_key[d["key"]]
        assert d["target"] == q.target
        if suite == "contact":
            continue
        if suite == "raw":
            pred = probe.parse_raw(d["text"])
        elif suite == "dna":
            pred = dna.parse(d["text"], q.prefix, cut=probe.cut_off(d))
        else:
            pred = probe.parse(d["text"], q.prefix, prefer_len=20 if suite == "recall" else 1, cut=probe.cut_off(d))
        assert pred == d["pred"], (suite, d["key"])
        d.update(
            valid=bool(pred) and d.get("stop") != "refusal",
            hit=bool(pred) and d.get("stop") != "refusal" and pred[0] == q.target[0],
            all10=pred[:10] == q.target[:10] and q.target[:10] not in q.prefix,
            all20=pred[:20] == q.target[:20] and q.target[:20] not in q.prefix,
        )
    return saved, expected


def contacts(upstream: Path, output: Path, saved: list[dict]) -> dict:
    # Reuse pure upstream functions without running its plotting script at import time.
    source = ast.parse((upstream / "scripts/analysis/contact_map.py").read_text())
    definitions = ast.Module(
        body=[
            node for node in source.body if isinstance(node, ast.FunctionDef) and node.name in ("grid", "band", "score")
        ],
        type_ignores=[],
    )
    namespace = {"np": np, "re": re, "hypergeom": hypergeom}
    exec(compile(definitions, "upstream_contact_functions", "exec"), namespace)
    truth = json.loads((output / "contact_truth.json").read_text())
    rows = []
    for d in saved:
        t = truth[d["protein"]]
        true = np.array([[int(c) for c in line] for line in t["map"]])
        pred, note = namespace["grid"](d["text"], len(true))
        row = {
            "protein": d["protein"],
            "pdb": d["pdb"],
            "group": "famous" if t["set"] == "famous" else "recent",
            "valid": pred is not None,
            "note": note,
        }
        if pred is not None:
            row.update({k: v for k, v in namespace["score"](pred, true, 24, 10**6).items() if k != "sep_probs"})
        rows.append(row)
    frame = pd.DataFrame(rows)
    frame.to_csv(output / "contact_scores.tsv", sep="\t", index=False)
    result = {}
    for name, group in frame.groupby("group"):
        valid = group[group.valid]
        n_pred = int(valid.n_pred.sum()) if "n_pred" in valid else 0
        correct = int(valid.n_correct.sum()) if "n_correct" in valid else 0
        result[name] = {
            "asked": len(group),
            "valid_grids": len(valid),
            "predicted_long_range": n_pred,
            "correct_long_range": correct,
            "precision": correct / n_pred if n_pred else None,
        }
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--upstream", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--suites", nargs="+", default=[*PROTEIN_SETS, "dna", "recall", "contact"])
    args = ap.parse_args()
    results, frames = {}, {}
    for suite in args.suites:
        saved, _expected = audit_suite(args.upstream, args.output, suite)
        frame = pd.DataFrame(saved)
        if suite in PROTEIN_SETS:
            results[suite] = protein_summary(frame)
            frames[suite] = frame
        elif suite == "recall":
            results[suite] = {
                "questions": len(frame),
                "valid": int(frame.valid.sum()),
                "longest": int(frame.exact_run.max()),
                **{f"at_least_{n}": int((frame.exact_run >= n).sum()) for n in (20, 50, 100, 200)},
            }
        elif suite == "dna":
            results[suite] = {}
            for (species, region), group in frame.groupby(["species", "region"]):
                valid = group[group.valid]
                pairs = valid.pivot(index="locus", columns="condition", values="hit").dropna().astype(float)
                results[suite][f"{species} {region}"] = {
                    "questions": len(group),
                    "valid": len(valid),
                    "paired_windows": len(pairs),
                    "real": float(valid[valid.condition == "real"].hit.mean() * 100),
                    "shuffled": float(valid[valid.condition == "shuffled"].hit.mean() * 100),
                    **interval((pairs.real - pairs.shuffled).to_numpy()),
                    "exact20_real": int(valid[valid.condition == "real"].all20.sum()),
                    "exact20_shuffled": int(valid[valid.condition == "shuffled"].all20.sum()),
                }
        else:
            results[suite] = contacts(args.upstream, args.output, saved)
        results[suite]["transport_attempts"] = sum(len(d["attempts"]) for d in saved)
        results[suite]["truncated_replies"] = int(frame.stop.isin(probe.CUT_OFF).sum())
        results[suite]["thinking_replies"] = sum("<|start_think|>" in d.get("raw_text", "") for d in saved)
        results[suite]["completion_tokens"] = sum(d.get("usage", {}).get("completion_tokens", 0) for d in saved)
    if "confirm" in frames and "small" in frames:
        valid = pd.concat([frames["confirm"], frames["small"]]).query("valid")
        fam = valid.groupby(["pfam", "condition"]).hit.mean().unstack().dropna()
        seq = pd.concat([pd.read_csv(args.upstream / "data" / PROTEIN_SETS[s], sep="\t") for s in ("confirm", "small")])
        size = seq.groupby("pfam").family_proteins.first().reindex(fam.index)
        rho, pvalue = spearmanr(fam.real - fam.shuffled, np.log10(size))
        results["family_size"] = {"families": len(fam), "rho": float(rho), "p_two_sided": float(pvalue)}
        fam.assign(size=size, gap=100 * (fam.real - fam.shuffled)).to_csv(args.output / "family_gaps.tsv", sep="\t")
    if "confirm" in frames and "d_set" in frames:
        families = set(frames["d_set"].pfam)
        results["main_same_d_families"] = protein_summary(frames["confirm"][frames["confirm"].pfam.isin(families)])
    if "pilot" in frames:
        seq = pd.read_csv(args.upstream / "data/sequences_pilot.tsv", sep="\t").set_index("accession")
        frame = frames["pilot"].copy()
        frame["literature"] = frame.accession.map(seq.lit_cell)
        results["pilot_by_literature"] = {name: protein_summary(group) for name, group in frame.groupby("literature")}
    (args.output / "summary.json").write_text(json.dumps(results, indent=2, allow_nan=False))
    if "confirm" in results:
        published = pd.read_csv(args.upstream / "results/confirm_tests.tsv", sep="\t")
        plm = pd.read_csv(args.upstream / "results/plm_summary.tsv", sep="\t")
        comparison = [
            {
                "name": row.model,
                "mode": "chat (published)",
                "gap": row.D_points,
                "low": row.D_ci_low,
                "high": row.D_ci_high,
            }
            for row in published.itertuples()
            if row.model != "claude-fable-5"
        ]
        comparison.extend(
            {"name": row["name"], "mode": "protein LM (published)", "gap": row["d"], "low": row["lo"], "high": row["hi"]}
            for row in plm.to_dict("records")
            if row["mode"] == "prefix"
        )
        for suite in ("confirm", "raw"):
            if suite in results:
                row = results[suite]
                comparison.append(
                    {
                        "name": "Mochi SFT " + ("chat" if suite == "confirm" else "raw"),
                        "mode": "measured",
                        **{k: row[k] for k in ("gap", "low", "high")},
                    }
                )
        comp = pd.DataFrame(comparison).sort_values("gap")
        comp.to_csv(args.output / "comparison.tsv", sep="\t", index=False)
        fig, ax = plt.subplots(figsize=(9, max(5, len(comp) * 0.27)))
        for y, row in enumerate(comp.itertuples()):
            measured = row.mode == "measured"
            ax.errorbar(
                row.gap,
                y,
                xerr=[[row.gap - row.low], [row.high - row.gap]],
                fmt="o",
                color="#b45309" if measured else "#64748b",
                capsize=2,
                markersize=6 if measured else 4,
            )
        ax.set_yticks(range(len(comp)), comp.name)
        ax.axvline(0, color="#cbd5e1", linewidth=1)
        ax.set_xlabel("Real - shuffled next-residue accuracy (percentage points; 95% family-bootstrap CI)")
        ax.set_title("Mochi SFT on the 668-protein benchmark")
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(args.output / "comparison.png", dpi=180)
        fig.savefig(args.output / "comparison.svg")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
