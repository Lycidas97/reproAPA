"""Split the existing L3 recurrence permutations into proximal and distal counts."""
import ast
import copy
import hashlib
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TABLES = ROOT / "notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables"
SIM_SOURCE = ROOT / "notes/2026-10-02-null-model-controls/recurrence.py"
COMMON_SOURCE = ROOT / "notes/2026-10-02-null-model-controls/common.py"
SEED_SOURCE = ROOT / "analysis_code/pas_rerun/downstream_runtime.py"
COLS = ["eligible_units", "series_pattern", "can_negative", "can_positive",
        "mean_GEX_bin", "APA_coverage_bin", "baseline_bin", "PAS_number_bin"]
DIRECTIONS = ["proximal", "distal"]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def numerical_function(path, name):
    tree = ast.parse(Path(path).read_text())
    functions = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name]
    assert len(functions) == 1, (path, name)
    namespace = dict(np=np, pd=pd, hashlib=hashlib)
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(path), "exec"), namespace)
    return namespace[name]


groups = numerical_function(COMMON_SOURCE, "groups")
bh = numerical_function(COMMON_SOURCE, "bh")
seed = numerical_function(SEED_SOURCE, "seed")


def directional_simulator():
    """Change only the recorded statistic; preserve the reviewed permutation loop."""
    tree = ast.parse(SIM_SOURCE.read_text())
    function = copy.deepcopy(next(n for n in tree.body
                                 if isinstance(n, ast.FunctionDef) and n.name == "simulate"))
    expected = "result.append((np.maximum(plus, minus) >= 2).sum(1))"
    candidates = [n for n in ast.walk(function)
                  if isinstance(n, ast.Expr) and ast.unparse(n) == expected]
    assert len(candidates) == 1, "Reviewed simulator output expression changed"
    candidates[0].value = ast.parse(
        "result.append(np.column_stack(((minus >= 2).sum(1), (plus >= 2).sum(1))))"
    ).body[0].value
    namespace = dict(np=np)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[function], type_ignores=[])),
                 str(SIM_SOURCE), "exec"), namespace)
    return namespace["simulate"]


simulate = directional_simulator()


def observed_members(frame):
    labels = frame.pivot(index="gene_id", columns="unit", values="event")
    assert labels.shape[1] <= 3, "Both directions may recur with four or more studies"
    proximal = labels.eq(-1).sum(axis=1).ge(2)
    distal = labels.eq(1).sum(axis=1).ge(2)
    assert not (proximal & distal).any()
    members = labels.index[proximal | distal]
    symbols = frame.drop_duplicates("gene_id").set_index("gene_id").gene_symbol
    return pd.DataFrame(dict(
        gene_id=members,
        gene_symbol=symbols.reindex(members).to_numpy(),
        direction=np.where(proximal.reindex(members), "proximal", "distal"),
        proximal_support=labels.eq(-1).sum(axis=1).reindex(members).to_numpy(),
        distal_support=labels.eq(1).sum(axis=1).reindex(members).to_numpy(),
        eligible_studies=labels.notna().sum(axis=1).reindex(members).to_numpy(),
    ))


def prepare_units(frame):
    genes = pd.Index(sorted(frame.gene_id.unique()))
    units, margin_rows = [], []
    probability = {d: [] for d in DIRECTIONS}
    for unit, part in frame.groupby("unit", sort=True):
        part = part.reset_index(drop=True)
        gi = genes.get_indexer(part.gene_id)
        cells = groups(part, COLS)
        units.append((part, gi, cells))
        probabilities = np.zeros((len(genes), 2))
        for j, ii in enumerate(cells):
            events = part.event.iloc[ii].to_numpy()
            counts = [(events == sign).sum() for sign in [-1, 0, 1]]
            probabilities[gi[ii], 0] = counts[0] / len(ii)
            probabilities[gi[ii], 1] = counts[2] / len(ii)
            margin_rows.append(dict(unit=unit, stratum=j, genes=len(ii),
                                    proximal=counts[0], non_event=counts[1], distal=counts[2]))
        for i, direction in enumerate(DIRECTIONS):
            probability[direction].append(probabilities[:, i])
    expected = {}
    for direction in DIRECTIONS:
        # Exact Poisson-binomial tail per gene; unavailable studies have probability zero.
        distribution = np.zeros((len(genes), len(units) + 1))
        distribution[:, 0] = 1
        for j, p in enumerate(probability[direction]):
            previous = distribution.copy()
            distribution[:, 0] = previous[:, 0] * (1 - p)
            for k in range(1, j + 2):
                distribution[:, k] = previous[:, k] * (1 - p) + previous[:, k - 1] * p
        expected[direction] = distribution[:, 2:].sum()
    return genes, units, expected, pd.DataFrame(margin_rows)


def save(frame, name):
    compression = {"method": "gzip", "mtime": 0} if name.endswith(".gz") else None
    frame.to_csv(HERE / name, sep="\t", index=False, na_rep="NA",
                 float_format="%.17g", compression=compression)


def main():
    inputs = [HERE / "protocol.json", TABLES / "recurrence_controls/recurrence_frozen_covariates.tsv.gz",
              TABLES / "recurrence_controls/recurrence_ladder.tsv",
              TABLES / "recurrence_controls/recurrence_draws.npz",
              TABLES / "recurrence/recurrence_members.tsv.gz", SIM_SOURCE, COMMON_SOURCE, SEED_SOURCE]
    hashes = {str(p.relative_to(ROOT)): digest(p) for p in inputs}
    audit_path = ROOT / "notes/2026-10-03-polyapipe-integrated-reference/statistics_results/v2_primary/recurrence_controls/audit.json"
    audit = json.loads(audit_path.read_text())
    assert digest(SIM_SOURCE) in audit["inputs"].values()
    assert digest(COMMON_SOURCE) in audit["inputs"].values()
    archive_path = ROOT / "table/Supplementary-Source-Data-20261009.zip"
    with zipfile.ZipFile(archive_path) as archive:
        for path in inputs[2:3] + inputs[4:5]:
            archived = archive.read("project/" + path.relative_to(ROOT).as_posix())
            assert hashlib.sha256(archived).hexdigest() == digest(path)
    frame = pd.read_csv(inputs[1], sep="\t", float_precision="round_trip")
    assert not frame.duplicated(["species", "contrast", "threshold", "unit", "gene_id"]).any()
    assert frame.metric.eq("WPP").all() and frame.event.isin([-1, 0, 1]).all()
    assert set(frame.threshold) == {0.02, 0.05}
    ladder = pd.read_csv(inputs[2], sep="\t", float_precision="round_trip")
    reference = np.load(inputs[3])
    original_members = pd.read_csv(inputs[4], sep="\t", float_precision="round_trip")
    original_members = original_members[original_members.metric.eq("WPP") &
                                        original_members.series_matching &
                                        original_members.variant.eq("primary")]
    rows, compositions, rosters, margins, verifications, draws = [], [], [], [], [], {}
    for (species, contrast, threshold), part in frame.groupby(["species", "contrast", "threshold"], sort=True):
        context = dict(species=species, contrast=contrast, threshold=threshold)
        members = observed_members(part)
        original = original_members[original_members.species.eq(species) &
                                    original_members.contrast.eq(contrast) &
                                    original_members.threshold.eq(threshold)]
        observed_pairs = set(zip(members.gene_id, members.direction))
        original_pairs = set(zip(original.gene_id, original.direction.map(
            {"negative": "proximal", "positive": "distal"})))
        assert observed_pairs == original_pairs
        genes, units, analytic, cells = prepare_units(part)
        run_seed = seed("recurrence_control", species, contrast, threshold, "L3_full")
        null = simulate(units, len(genes), 10000, np.random.default_rng(run_seed))
        key = "|".join(map(str, (species, contrast, threshold, "L3_full")))
        assert np.array_equal(null.sum(axis=1), reference[key]), f"Original draw replay failed: {key}"
        old = ladder[ladder.species.eq(species) & ladder.contrast.eq(contrast) &
                     ladder.threshold.eq(threshold) & ladder.model.eq("L3_full")].iloc[0]
        assert len(members) == old.observed
        assert np.isclose(null.sum(axis=1).mean(), old.expected, atol=1e-12, rtol=0)
        counts = members.direction.value_counts().reindex(DIRECTIONS, fill_value=0)
        for i, direction in enumerate(DIRECTIONS):
            dnull, observed = null[:, i], int(counts[direction])
            mean, sd = float(dnull.mean()), float(dnull.std(ddof=1))
            assert abs(mean - analytic[direction]) <= 6 * sd / np.sqrt(len(dnull)) + 1e-10
            tail = int((dnull >= observed).sum())
            rows.append(context | dict(direction=direction, observed=observed,
                null_mean=mean, analytic_null_mean=analytic[direction],
                null_p025=float(np.quantile(dnull, .025)), null_p975=float(np.quantile(dnull, .975)),
                null_sd=sd, fold_enrichment=observed / mean, excess=observed - mean,
                excess_fraction_of_total=(observed - mean) / (len(members) - old.expected),
                p_enrichment=(1 + tail) / (len(dnull) + 1), null_tail_count=tail,
                permutations=len(dnull), seed=run_seed, measurable_genes=len(genes),
                contributing_studies=len(units)))
            draws[key + "|" + direction] = dnull
        assert (null.sum(axis=1) > 0).all()
        fraction = null[:, 1] / null.sum(axis=1)
        observed_fraction = float(counts["distal"] / len(members))
        p_lower = (1 + (fraction <= observed_fraction).sum()) / (len(fraction) + 1)
        p_upper = (1 + (fraction >= observed_fraction).sum()) / (len(fraction) + 1)
        compositions.append(context | dict(observed_distal_fraction=observed_fraction,
            null_mean_distal_fraction=float(fraction.mean()),
            null_distal_fraction_p025=float(np.quantile(fraction, .025)),
            null_distal_fraction_p975=float(np.quantile(fraction, .975)),
            distal_fraction_difference=observed_fraction - fraction.mean(),
            p_two_sided=min(1, 2 * min(p_lower, p_upper)), p_lower=p_lower, p_upper=p_upper))
        rosters.append(members.assign(**context))
        margins.append(cells.assign(**context))
        verifications.append(context | dict(observed_roster_exact_match=True,
            original_total_draws_exact_match=True, matched_draws=len(null),
            observed_total=len(members), replay_total_null_mean=float(null.sum(axis=1).mean()),
            reference_total_null_mean=old.expected, analytic_means_verified=True))
        print(key, "observed", counts.to_dict(), "null", null.mean(axis=0).tolist(), flush=True)
    result = pd.DataFrame(rows)
    result["q_threshold_8"] = result.groupby("threshold").p_enrichment.transform(bh)
    result["q_joint_16"] = bh(result.p_enrichment)
    composition = pd.DataFrame(compositions)
    composition["q_threshold_4"] = composition.groupby("threshold").p_two_sided.transform(bh)
    composition["q_joint_8"] = bh(composition.p_two_sided)
    save(result, "directional_recurrence.tsv")
    save(composition, "direction_composition.tsv")
    save(pd.concat(rosters, ignore_index=True), "observed_members.tsv.gz")
    save(pd.concat(margins, ignore_index=True), "study_stratum_event_margins.tsv.gz")
    save(pd.DataFrame(verifications), "verification.tsv")
    np.savez_compressed(HERE / "directional_null_draws.npz", **draws)
    assert all(digest(ROOT / p) == value for p, value in hashes.items()), "An input changed during analysis"
    outputs = [p for p in HERE.iterdir() if p.name in ["directional_recurrence.tsv", "direction_composition.tsv",
               "observed_members.tsv.gz", "study_stratum_event_margins.tsv.gz", "verification.tsv", "directional_null_draws.npz"]]
    manifest = dict(status="PASS", input_sha256=hashes, output_sha256={p.name: digest(p) for p in outputs},
        script_sha256=digest(Path(__file__)), original_audit_path=str(audit_path.relative_to(ROOT)),
        original_audit_sha256=digest(audit_path), archive_path=str(archive_path.relative_to(ROOT)),
        archive_sha256=digest(archive_path), python=sys.version, python_executable=sys.executable,
        numpy_version=np.__version__, pandas_version=pd.__version__, total_exactly_replayed_draws=80000,
        statistic_change="Only result.append in the reviewed simulate function; permutation loop and random seeds preserved.",
        units="Genes; strata are within study. Null intervals are permutation reference ranges, not confidence intervals.")
    (HERE / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
