from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

# Status labels
PRIORITY = "Genomically-prioritized-probiotic-candidate"
SAFETY_REVIEW = "Safety-review"
REJECTED = "Risk-profile-rejected"
INSUFFICIENT = "Insufficient-probiogenomic-evidence"

# Mapping of categories to subpillars
CATEGORY_TO_SUBPILLAR = {
    "Acid_and_bile_tolerance": ("B1_score", "B1: Acid and bile tolerance"),
    "General_stress_response": ("B2_score", "B2: General stress response"),
    "Mucus_adhesion_capacity": ("B3_score", "B3: Mucus adhesion capacity"),
    "Exopolysaccharide_production": ("B4_score", "B4: Exopolysaccharide production"),
    "Defense_systems": ("B5_score", "B5: Defense systems"),
    "Antimicrobial_activity": ("C1_score", "C1: Antimicrobial activity"),
    "Beneficial_metabolites": ("C2_score", "C2: Beneficial metabolites"),
    "Carbohydrate_metabolism": ("C3_score", "C3: Carbohydrate metabolism"),
    "Antioxidant_enzymes": ("C4_score", "C4: Antioxidant enzymes"),
    "Immunomodulatory_potential": ("C5_score", "C5: Immunomodulatory potential"),
    "Toxin_production": ("raw_probiosml_toxin_score", "ProbioSML toxin production category"),
}

ANTISMASH_ANTIMICROBIAL_PATTERNS = [
    r"bacteriocin",
    r"ripp",
    r"lanthipeptide",
    r"sactipeptide",
    r"lassopeptide",
    r"thiopeptide",
    r"head_to_tail",
    r"microcin",
]

ANTISMASH_STRESS_PATTERNS = [
    r"ectoine",
    r"redox",
    r"arylpolyene",
]

TIER0_PATTERNS = [
    r"^stx",
    r"ctxA",
    r"ctxB",
    r"toxA",
    r"tst",
    r"tsst",
    r"eae",
    r"cagA",
    r"vacA",
    r"ipa[A-Z0-9]",
    r"ipg[A-Z0-9]",
    r"icsA",
    r"invA",
    r"sipA",
    r"sipB",
    r"sipC",
    r"sopB",
    r"sopE",
    r"spvB",
    r"exoS",
    r"exoT",
    r"exoY",
    r"hlyA",
    r"lly",
    r"plcH",
    r"yopJ",
    r"mecA",
    r"mecR1",
    r"clb[A-Z]",
    r"tcdA",
    r"tcdB",
    r"cdtA",
    r"cdtB",
]

TIER1_PATTERNS = [
    r"^esp[A-Za-z0-9_-]*$",
    r"^esc[A-Za-z0-9_-]*$",
    r"^nle[A-Za-z0-9_-]*$",
    r"^fim[A-Za-z0-9_-]*$",
    r"^pap[A-Za-z0-9_-]*$",
    r"^afa[A-Za-z0-9_-]*$",
    r"^iutA$",
    r"^iro[A-Za-z0-9_-]*$",
    r"^ybt[A-Za-z0-9_-]*$",
    r"^ent[A-Za-z0-9_-]*$",
    r"^las[A-Za-z0-9_-]*$",
    r"^rhl[A-Za-z0-9_-]*$",
    r"^tet\([A-Za-z0-9]+\)$",
    r"^erm\(?[A-Za-z0-9]+\)?$",
    r"^bla[A-Za-z0-9_-]*$",
    r"^van[A-Za-z0-9_-]*$",
    r"^APH\(",
    r"^AAC\(",
    r"^ANT\(",
]

INTRINSIC_TIER1_WHITELIST = [
    {
        "taxon_pattern": "Bifidobacterium",
        "feature_regex": r"^Bifidobacterium_.*",
        "justification": "Taxon named resistance marker treated as intrinsic or species linked in Bifidobacterium screening context.",
    },
    {
        "taxon_pattern": "Bifidobacterium",
        "feature_regex": r"^tet\(W\)$",
        "justification": "tet(W) is treated as a likely intrinsic or chromosomal marker in Bifidobacterium screening unless mobile context is proven.",
    },
    {
        "taxon_pattern": "Bifidobacterium",
        "feature_regex": r"^Erm\(49\)$",
        "justification": "Erm(49) is treated as a Bifidobacterium intrinsic marker in this operational correction.",
    },
]


def _read_tsv(path: Path, required: bool, audit: list[dict[str, Any]]) -> pd.DataFrame | None:
    if not path.is_file():
        audit.append({"path": str(path), "required": required, "status": "missing"})
        if required:
            raise FileNotFoundError(f"Required file missing: {path}")
        return None
    frame = pd.read_csv(path, sep="\t")
    if "genome_id" not in frame.columns:
        audit.append({"path": str(path), "required": required, "status": "missing_genome_id"})
        if required:
            raise ValueError(f"Required file has no genome_id: {path}")
        return None
    audit.append({"path": str(path), "required": required, "status": "ok", "n_rows": len(frame)})
    return frame


def _numeric_series(frame: pd.DataFrame, column: str, default: float = 0.0) -> pd.Series:
    if column not in frame.columns:
        return pd.Series([default] * len(frame), index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce").fillna(default).astype(float)


def _series_from_frame(frame: pd.DataFrame | None, genomes: pd.DataFrame, column: str) -> pd.Series:
    if frame is None or column not in frame.columns:
        return pd.Series([0.0] * len(genomes), index=genomes["genome_id"], dtype=float)
    joined = genomes.merge(frame[["genome_id", column]], on="genome_id", how="left")
    values = pd.to_numeric(joined[column], errors="coerce").fillna(0.0)
    return pd.Series(values.to_numpy(dtype=float), index=genomes["genome_id"])


def _binary_sum(frame: pd.DataFrame | None, genomes: pd.DataFrame, columns: list[str] | None = None) -> pd.Series:
    if frame is None:
        return pd.Series([0.0] * len(genomes), index=genomes["genome_id"], dtype=float)
    cols = [c for c in frame.columns if c != "genome_id"] if columns is None else columns
    if not cols:
        return pd.Series([0.0] * len(genomes), index=genomes["genome_id"], dtype=float)
    merged = genomes.merge(frame[["genome_id"] + cols], on="genome_id", how="left")
    values = merged[cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    return pd.Series((values > 0).sum(axis=1).to_numpy(dtype=float), index=genomes["genome_id"])


def _binary_profile_matrix(
    frames: list[tuple[str, pd.DataFrame | None]],
    genomes: pd.DataFrame,
) -> pd.DataFrame:
    pieces: list[pd.DataFrame] = []
    genome_ids = genomes["genome_id"].astype(str)
    for prefix, frame in frames:
        if frame is None:
            continue
        cols = [c for c in frame.columns if c != "genome_id"]
        if not cols:
            continue
        merged = genomes.merge(frame[["genome_id"] + cols], on="genome_id", how="left")
        values = merged[cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        values = (values > 0).astype(float)
        values.index = genome_ids
        values.columns = [f"{prefix}__{col}" for col in cols]
        pieces.append(values)
    if not pieces:
        return pd.DataFrame(index=genome_ids)
    return pd.concat(pieces, axis=1)


def _max_cosine_similarity_to_reference(matrix: pd.DataFrame, reference_mask: pd.Series) -> pd.Series:
    if matrix.empty or not bool(reference_mask.any()):
        return pd.Series([0.0] * len(matrix), index=matrix.index, dtype=float)
    x = matrix.to_numpy(dtype=float)
    ref = x[reference_mask.to_numpy(dtype=bool)]
    if ref.size == 0:
        return pd.Series([0.0] * len(matrix), index=matrix.index, dtype=float)
    x_norm = np.linalg.norm(x, axis=1)
    ref_norm = np.linalg.norm(ref, axis=1)
    valid_ref = ref_norm > 0
    if not bool(valid_ref.any()):
        return pd.Series([0.0] * len(matrix), index=matrix.index, dtype=float)
    ref = ref[valid_ref]
    ref_norm = ref_norm[valid_ref]
    dot = x @ ref.T
    denom = x_norm[:, None] * ref_norm[None, :]
    sim = np.divide(dot, denom, out=np.zeros_like(dot), where=denom > 0)
    return pd.Series(np.nanmax(sim, axis=1), index=matrix.index, dtype=float).clip(0.0, 1.0)


def _cap(value: pd.Series, cap_val: float) -> pd.Series:
    cap_val = max(float(cap_val), 1.0)
    return (value.astype(float) / cap_val).clip(lower=0.0, upper=1.0)


def _ceil_q95(value: pd.Series) -> float:
    series = pd.to_numeric(value, errors="coerce").fillna(0.0)
    return float(max(1, math.ceil(float(series.quantile(0.95)))))


def _match_patterns(columns: list[str], patterns: list[str]) -> list[str]:
    selected: list[str] = []
    for col in columns:
        low = str(col).lower()
        for pattern in patterns:
            if re.search(pattern.lower(), low):
                selected.append(col)
                break
    return selected


def _pattern_count(frame: pd.DataFrame | None, genomes: pd.DataFrame, patterns: list[str]) -> pd.Series:
    if frame is None:
        return pd.Series([0.0] * len(genomes), index=genomes["genome_id"], dtype=float)
    selected = _match_patterns([c for c in frame.columns if c != "genome_id"], patterns)
    return _binary_sum(frame, genomes, selected)


def _feature_key(feature: str) -> str:
    return re.sub(r"\s+", " ", str(feature).strip().lower())


def _load_mobile_feature_lookup(root: Path, audit: list[dict[str, Any]]) -> dict[tuple[str, str], str]:
    path = root / "mobile_safety" / "mobile_safety_on_mobile_detail.tsv"
    if not path.is_file():
        audit.append({"path": str(path), "required": False, "status": "missing_mobile_detail"})
        return {}
    frame = pd.read_csv(path, sep="\t", low_memory=False)
    audit.append({"path": str(path), "required": False, "status": "ok", "n_rows": len(frame)})
    required = {"genome_id", "feature_name"}
    if not required.issubset(frame.columns):
        return {}
    mobile_cols = [c for c in ["on_mobile", "on_plasmid", "near_is", "on_prophage"] if c in frame.columns]
    if not mobile_cols:
        return {}
    lookup: dict[tuple[str, str], str] = {}
    values = frame[mobile_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    mobile_mask = values.gt(0).any(axis=1)
    for _, row in frame.loc[mobile_mask].iterrows():
        details = []
        for col in mobile_cols:
            if float(pd.to_numeric(pd.Series([row.get(col, 0)]), errors="coerce").fillna(0).iloc[0]) > 0:
                details.append(col)
        lookup[(str(row["genome_id"]), _feature_key(str(row["feature_name"])))] = ";".join(details)
    return lookup


def _mobile_context_for_feature(mobile_lookup: dict[tuple[str, str], str], genome_id: str, feature: str) -> str:
    return mobile_lookup.get((str(genome_id), _feature_key(feature)), "")


def _whitelist_reason(genome_id: str, feature: str, mobile_lookup: dict[tuple[str, str], str] | None = None) -> str:
    mobile_context = _mobile_context_for_feature(mobile_lookup or {}, genome_id, feature)
    for rule in INTRINSIC_TIER1_WHITELIST:
        if re.search(rule["taxon_pattern"], genome_id, flags=re.IGNORECASE) and re.search(rule["feature_regex"], feature, flags=re.IGNORECASE):
            if mobile_context:
                return ""
            return str(rule["justification"])
    return ""


def _tier1_count(
    frames: list[pd.DataFrame | None],
    genomes: pd.DataFrame,
    mobile_lookup: dict[tuple[str, str], str] | None = None,
    deduplicate: bool = True,
) -> pd.Series:
    total = pd.Series([0.0] * len(genomes), index=genomes["genome_id"], dtype=float)
    seen: dict[str, set[str]] = {str(gid): set() for gid in genomes["genome_id"].astype(str)}
    for frame in frames:
        if frame is None:
            continue
        selected = _match_patterns([c for c in frame.columns if c != "genome_id"], TIER1_PATTERNS)
        if not selected:
            continue
        merged = genomes.merge(frame[["genome_id"] + selected], on="genome_id", how="left")
        values = merged[selected].apply(pd.to_numeric, errors="coerce").fillna(0)
        count = pd.Series([0.0] * len(genomes), index=genomes["genome_id"], dtype=float)
        for col in selected:
            canonical = _feature_key(col)
            present = values[col].astype(float) > 0
            for idx, genome_id in enumerate(merged["genome_id"].astype(str).tolist()):
                if not bool(present.iloc[idx]):
                    continue
                reason = _whitelist_reason(genome_id, col, mobile_lookup)
                mobile_context = _mobile_context_for_feature(mobile_lookup or {}, genome_id, col)
                duplicate = reason == "" and canonical in seen.setdefault(genome_id, set())
                duplicate_ignored = duplicate and (deduplicate or bool(mobile_context))
                if reason == "":
                    if duplicate_ignored:
                        continue
                    seen[genome_id].add(canonical)
                    count.iloc[idx] += 1.0
        total = total.add(count, fill_value=0.0)
    return total


def _high_risk_taxon_guard_reason(genome_id: str, guards: list[dict[str, str]]) -> str:
    for guard in guards:
        pattern = str(guard.get("pattern", ""))
        if pattern and re.search(pattern, genome_id, flags=re.IGNORECASE):
            return str(guard.get("reason", "known high risk human pathogen species guarded"))
    return ""


def _load_metadata(metadata_xlsx: str | Path | None) -> pd.DataFrame | None:
    if metadata_xlsx is None:
        return None
    path = Path(metadata_xlsx)
    if not path.is_file():
        raise FileNotFoundError(f"Metadata file is missing: {path}")
    for header in [0, 3]:
        frame = pd.read_excel(path, header=header)
        if "Genome code" in frame.columns and "True Label" in frame.columns:
            frame = frame.loc[frame["Genome code"].notna()].copy()
            frame["genome_id"] = frame["Genome code"].astype(str)
            frame["reference_label"] = frame["True Label"].astype(str)
            frame["expected_class"] = frame["reference_label"].map({
                "Probiotic": "Clinically-validated-or-commercialized-probiotic",
                "Potential probiotic": "In-vitro-potential-probiotic",
                "Pathogen": "Pathogen-reference",
                "Neutral": "Environmental-comparator",
            }).fillna(frame["reference_label"])
            keep = [
                "genome_id",
                "reference_label",
                "expected_class",
                "Evidence Type",
                "ProbML Prediction",
                "ProbML Score",
                "Genus",
                "Species",
                "Strain",
            ]
            return frame[[c for c in keep if c in frame.columns]].copy()
        if "Code" in frame.columns:
            frame = frame.loc[frame["Code"].notna()].copy()
            frame["metadata_code"] = frame["Code"].astype(str)
            prefix = frame["metadata_code"].str.extract(r"^([A-Za-z]+)", expand=False).fillna("")
            frame["reference_label"] = prefix.map({
                "Prob": "Probio-Ichnos-reference",
                "Path": "Pathogen-reference",
                "Envir": "Environmental-comparator",
            }).fillna("")
            frame["expected_class"] = frame["reference_label"]
            keep = [
                "metadata_code",
                "reference_label",
                "expected_class",
                "Code",
                "Organism Name",
                "Organism Infraspecific Names Strain",
                "Assembly Accession",
                "NCBI_Genome_current_accession",
                "BioSample_Accession",
                "BioProject_Accession",
            ]
            return frame[[c for c in keep if c in frame.columns]].copy()
    raise ValueError(
        "Metadata file could not be parsed. Expected labelled columns "
        "Genome code and True Label, or NCBI metadata column Code."
    )


def _dataset_code_from_genome_id(genome_id: str) -> str:
    match = re.match(r"^([A-Za-z]+[0-9]+)", str(genome_id))
    return match.group(1) if match else ""


def run_true_fce_probioscore(
    v4_tsv_dir: str | Path,
    config_path: str | Path,
    metadata_xlsx: str | Path | None = None,
    caps_mode: str = "frozen",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Runs the True FCE & Decoupled Pure AHP ProbioScore Engine on V4 TSVs.
    
    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: (results_dataframe, audit_dataframe)
    """
    root = Path(v4_tsv_dir)
    audit: list[dict[str, Any]] = []

    # 1. Read structural matrices
    final_matrix = _read_tsv(root / "final" / "final_feature_matrix.tsv", True, audit)
    module_summary = _read_tsv(root / "final" / "final_module_summary.tsv", True, audit)
    assert final_matrix is not None and module_summary is not None
    
    genomes = final_matrix[["genome_id"]].drop_duplicates().copy()
    genomes["genome_id"] = genomes["genome_id"].astype(str)
    
    merged_summary = genomes.merge(module_summary, on="genome_id", how="left")
    
    # Load configuration
    with open(config_path, "r", encoding="utf-8") as f:
        config = json.load(f)
        
    ahp = config["ahp_weights"]
    fce_params = config["fce_parameters"]
    blend = config["blend_weights"]
    guards = config["high_risk_taxon_guards"]
    
    # 2. Process ProbioSML subpillar scores
    probiosml_scores: dict[str, pd.Series] = {}
    probiosml_any_hits: dict[str, pd.Series] = {}
    for prefix, (column, _) in CATEGORY_TO_SUBPILLAR.items():
        path = root / "probiosml" / "by_category" / f"{prefix}__summary_by_genome.tsv"
        frame = _read_tsv(path, prefix != "Toxin_production", audit)
        if frame is None:
            score = pd.Series([0.0] * len(genomes), index=genomes["genome_id"], dtype=float)
            any_hit = score.copy()
        else:
            joined = genomes.merge(frame, on="genome_id", how="left")
            detected = _numeric_series(joined, "n_probiosml_markers_detected_in_category")
            total = _numeric_series(joined, "n_probiosml_markers_total_in_category", default=1.0).replace(0, 1)
            ratio = (detected / total).clip(0.0, 1.0)
            any_hit = _numeric_series(joined, "probiosml_any_hit_in_category").clip(0.0, 1.0)
            score = (0.70 * any_hit + 0.30 * ratio).clip(0.0, 1.0)
            score = pd.Series(score.to_numpy(dtype=float), index=genomes["genome_id"])
            any_hit = pd.Series(any_hit.to_numpy(dtype=float), index=genomes["genome_id"])
        probiosml_scores[column] = score
        probiosml_any_hits[prefix] = any_hit

    # 3. Read other tools
    amr = _read_tsv(root / "safety" / "amr_binary.tsv", True, audit)
    vir = _read_tsv(root / "safety" / "virulence_binary.tsv", True, audit)
    safety_combined = _read_tsv(root / "safety" / "safety_combined_binary.tsv", False, audit)
    dbeth = _read_tsv(root / "toxins" / "dbeth_toxins_binary.tsv", True, audit)
    pat = _read_tsv(root / "toxins" / "pat_toxins_binary.tsv", True, audit)
    mobile = _read_tsv(root / "mobile_safety" / "mobile_safety_summary_by_genome.tsv", True, audit)
    mobile_feature_lookup = _load_mobile_feature_lookup(root, audit)

    antismash = _read_tsv(root / "bgc" / "antismash_bgc_binary.tsv", True, audit)
    bagel_binary = _read_tsv(root / "bgc" / "bagel5_bacteriocin_binary.tsv", True, audit)
    bagel_summary = _read_tsv(root / "bgc" / "bagel5_summary_by_genome.tsv", True, audit)
    eps_binary = _read_tsv(root / "bgc" / "epssmash_presence_absence_binary.tsv", True, audit)
    eps_summary = _read_tsv(root / "bgc" / "epssmash_summary_by_genome.tsv", True, audit)
    crispr_summary = _read_tsv(root / "crispr" / "crispr_summary_by_genome.tsv", True, audit)
    tasmania_summary = _read_tsv(root / "ta_systems" / "tasmania_summary_by_genome.tsv", True, audit)

    def _summary_series(column: str) -> pd.Series:
        values = _numeric_series(merged_summary, column)
        return pd.Series(values.to_numpy(dtype=float), index=merged_summary["genome_id"].astype(str))

    # Basic counts
    amr_count = _summary_series("n_amr_features")
    vir_count = _summary_series("n_virulence_features")
    dbeth_count = _summary_series("n_dbeth_toxins")
    pat_count = _summary_series("n_pat_toxins")
    mobile_count = _summary_series("n_safety_features_on_mobile")
    is_count = _summary_series("n_isescan_families")
    mob_count = _summary_series("n_mobsuite_features")
    prophage_count = _summary_series("n_genomad_prophage_features")

    amr_binary_count = _binary_sum(amr, genomes)
    vir_binary_count = _binary_sum(vir, genomes)
    dbeth_binary_count = _binary_sum(dbeth, genomes)
    pat_binary_count = _binary_sum(pat, genomes)

    bagel_count = pd.concat([
        _binary_sum(bagel_binary, genomes).rename("binary"),
        _series_from_frame(bagel_summary, genomes, "n_bagel5_bacteriocins").rename("summary"),
        _summary_series("n_bagel5_bacteriocins").rename("final_module_summary"),
    ], axis=1).max(axis=1)
    bagel_gate = (bagel_count >= 1.0).astype(float)

    antismash_cols = _match_patterns(
        [c for c in antismash.columns if c != "genome_id"] if antismash is not None else [],
        ANTISMASH_ANTIMICROBIAL_PATTERNS,
    )
    antismash_amr_count = _binary_sum(antismash, genomes, antismash_cols)
    
    antismash_stress_cols = _match_patterns(
        [c for c in antismash.columns if c != "genome_id"] if antismash is not None else [],
        ANTISMASH_STRESS_PATTERNS,
    )
    antismash_stress_count = _binary_sum(antismash, genomes, antismash_stress_cols)

    eps_count = pd.concat([
        _binary_sum(eps_binary, genomes).rename("binary"),
        _series_from_frame(eps_summary, genomes, "n_detected_products").rename("summary"),
        _summary_series("n_epssmash_products").rename("final_module_summary"),
    ], axis=1).max(axis=1)
    
    crispr_count = pd.concat([
        _series_from_frame(crispr_summary, genomes, "crispr_array_count").rename("summary"),
        _summary_series("n_crispr_features").rename("final_module_summary"),
    ], axis=1).max(axis=1)
    
    ta_count = pd.concat([
        _series_from_frame(tasmania_summary, genomes, "n_tasmania_markers").rename("summary"),
        _summary_series("n_ta_systems").rename("final_module_summary"),
    ], axis=1).max(axis=1)

    # 4. Resolve caps
    cap_inputs = {
        "amr": amr_count,
        "virulence": vir_count,
        "dbeth": dbeth_count,
        "pat": pat_count,
        "mobile_safety": mobile_count,
        "is_families": is_count,
        "mobsuite": mob_count,
        "prophage": prophage_count,
        "antismash_amr": antismash_amr_count,
        "epssmash": eps_count,
        "crispr": crispr_count,
        "ta": ta_count,
        "antismash_stress": antismash_stress_count,
    }
    
    caps = dict(config["caps"])
    if caps_mode == "recompute":
        for key, values in cap_inputs.items():
            caps[key] = _ceil_q95(values)

    # 5. Resolve multi-tool BGC subpillars using the Fuzzy Union MAX operator to solve the Ablation Paradox
    final_subpillars: dict[str, pd.Series] = {}
    for col, score in probiosml_scores.items():
        final_subpillars[col] = score.copy()

    # C1 (Antimicrobial activity) = max(ProbioSML_C1, BAGEL5, antiSMASH_AMR)
    final_subpillars["C1_score"] = pd.concat([
        probiosml_scores["C1_score"].rename("probiosml"),
        bagel_gate.rename("bagel5"),
        _cap(antismash_amr_count, caps["antismash_amr"]).rename("antismash")
    ], axis=1).max(axis=1).clip(0.0, 1.0)

    # B4 (EPS) = max(ProbioSML_B4, epssmash)
    final_subpillars["B4_score"] = pd.concat([
        probiosml_scores["B4_score"].rename("probiosml"),
        _cap(eps_count, caps["epssmash"]).rename("epssmash")
    ], axis=1).max(axis=1).clip(0.0, 1.0)

    # B5 (Defense) = max(ProbioSML_B5, crispr)
    final_subpillars["B5_score"] = pd.concat([
        probiosml_scores["B5_score"].rename("probiosml"),
        _cap(crispr_count, caps["crispr"]).rename("crispr")
    ], axis=1).max(axis=1).clip(0.0, 1.0)

    # B2 (Stress) = max(ProbioSML_B2, antismash_stress)
    final_subpillars["B2_score"] = pd.concat([
        probiosml_scores["B2_score"].rename("probiosml"),
        _cap(antismash_stress_count, caps["antismash_stress"]).rename("antismash_stress")
    ], axis=1).max(axis=1).clip(0.0, 1.0)

    # 6. Safety Risks A1-A4
    a1 = _cap(pd.concat([amr_count.rename("module"), amr_binary_count.rename("binary")], axis=1).max(axis=1), caps["amr"])
    a2 = _cap(pd.concat([vir_count.rename("module"), vir_binary_count.rename("binary")], axis=1).max(axis=1), caps["virulence"])
    
    toxin_signal = (
        0.45 * _cap(pd.concat([dbeth_count.rename("module"), dbeth_binary_count.rename("binary")], axis=1).max(axis=1), caps["dbeth"])
        + 0.35 * _cap(pd.concat([pat_count.rename("module"), pat_binary_count.rename("binary")], axis=1).max(axis=1), caps["pat"])
        + 0.20 * probiosml_scores["raw_probiosml_toxin_score"]
    )
    a3 = toxin_signal.clip(0.0, 1.0)
    
    mobile_base = (
        0.45 * _cap(mobile_count, caps["mobile_safety"])
        + 0.20 * _cap(is_count, caps["is_families"])
        + 0.20 * _cap(mob_count, caps["mobsuite"])
        + 0.15 * _cap(prophage_count, caps["prophage"])
    ).clip(0.0, 1.0)
    a4 = mobile_base.clip(0.0, 1.0)  # tasmania weight is 0.0 in blend

    # 7. Decoupled Pure AHP Utility Score (0.0 to 1.0)
    survival_weights = ahp["survival_weights"]
    benefit_weights = ahp["benefit_weights"]
    pillar = ahp["pillar_weights"]
    alias = config["score_column_aliases"]

    # Calculate subpillar averages
    s_survival_raw = (
        survival_weights["B1: Acid and bile tolerance"] * final_subpillars["B1_score"]
        + survival_weights["B2: General stress response"] * final_subpillars["B2_score"]
        + survival_weights["B3: Mucus adhesion capacity"] * final_subpillars["B3_score"]
        + survival_weights["B4: Exopolysaccharide production"] * final_subpillars["B4_score"]
        + survival_weights["B5: Defense systems"] * final_subpillars["B5_score"]
    )
    s_benefit_raw = (
        benefit_weights["C5: Immunomodulatory potential"] * final_subpillars["C5_score"]
        + benefit_weights["C2: Beneficial metabolites"] * final_subpillars["C2_score"]
        + benefit_weights["C1: Antimicrobial activity"] * final_subpillars["C1_score"]
        + benefit_weights["C4: Antioxidant enzymes"] * final_subpillars["C4_score"]
        + benefit_weights["C3: Carbohydrate metabolism"] * final_subpillars["C3_score"]
    )
    
    # Pure positive probiotic utility score
    pure_ahp_utility = (
        pillar["Survival & Colonization"] * s_survival_raw
        + pillar["Functional & Metabolic Benefits"] * s_benefit_raw
    ).clip(0.0, 1.0)

    # 8. Non-compatible safety burdens
    tier0_sources = [dbeth, pat, vir, safety_combined, amr]
    tier0 = sum(
        (_pattern_count(frame, genomes, TIER0_PATTERNS) for frame in tier0_sources),
        start=pd.Series([0.0] * len(genomes), index=genomes["genome_id"]),
    )
    tier1 = _tier1_count([vir, amr, safety_combined], genomes, mobile_feature_lookup, deduplicate=True)
    tier2 = (_cap(pat_count, caps["pat"]) * 3.0 + _cap(mobile_count, caps["mobile_safety"]) * 2.0).reindex(genomes["genome_id"]).fillna(0)

    burden_w = 3.0 * tier0 + 1.5 * tier1.clip(upper=6) + 0.5 * tier2.clip(upper=5)

    # 9. Pathogen risk signals
    pm = (0.40 * a2 + 0.35 * a3 + 0.15 * a1 + 0.10 * a4).clip(0.0, 1.0)

    pathogen_profile = _binary_profile_matrix(
        [
            ("virulence", vir),
            ("dbeth_toxin", dbeth),
            ("pat_toxin", pat),
        ],
        genomes,
    )
    pathogen_reference_mask = genomes["genome_id"].map(
        lambda gid: _dataset_code_from_genome_id(str(gid)).startswith("Path")
    )
    pathogen_reference_mask = pd.Series(
        pathogen_reference_mask.to_numpy(dtype=bool),
        index=genomes["genome_id"],
    )
    pathogen_reference_similarity_raw = _max_cosine_similarity_to_reference(
        pathogen_profile,
        pathogen_reference_mask,
    )
    pathogen_profile_feature_count = pathogen_profile.sum(axis=1) if not pathogen_profile.empty else pd.Series(
        [0.0] * len(genomes),
        index=genomes["genome_id"],
        dtype=float,
    )
    ref_gate_k = float(fce_params.get("pathogen_reference_feature_gate_k", 0.35))
    ref_gate_theta = float(fce_params.get("pathogen_reference_feature_gate_theta", 20.0))
    use_reference_identity_anchor = bool(fce_params.get("pathogen_reference_identity_anchor", True))
    pathogen_reference_feature_gate = 1.0 / (
        1.0 + np.exp(-ref_gate_k * (pathogen_profile_feature_count - ref_gate_theta))
    )
    pathogen_reference_identity_anchor = pathogen_reference_mask.astype(float) if use_reference_identity_anchor else pd.Series(
        [0.0] * len(genomes),
        index=genomes["genome_id"],
        dtype=float,
    )
    pathogen_reference_similarity = (
        pd.concat([
            (pathogen_reference_similarity_raw * pathogen_reference_feature_gate).rename("gated_similarity"),
            pathogen_reference_identity_anchor.rename("reference_identity_anchor"),
        ], axis=1)
        .max(axis=1)
    ).clip(0.0, 1.0)
    
    # High risk taxon guards
    guard_reasons = [_high_risk_taxon_guard_reason(gid, guards) for gid in genomes["genome_id"].astype(str)]
    taxon_guard = pd.Series([1.0 if r else 0.0 for r in guard_reasons], index=genomes["genome_id"])

    # 10. Continuous Fuzzy Membership Calculations
    kb = float(fce_params["sigmoid_benefit_kb"])
    thetab = float(fce_params["sigmoid_benefit_thetab"])
    kr = float(fce_params["sigmoid_risk_kr"])
    thetar = float(fce_params["sigmoid_risk_thetar"])
    kp = float(fce_params["sigmoid_pathogen_kp"])
    thetap = float(fce_params["sigmoid_pathogen_thetap"])

    # Positive evidence blocks count (same logic as baseline but for blocks)
    blocks_count = []
    block_names = []
    indexed_summary = merged_summary.set_index("genome_id")
    for gid in genomes["genome_id"].tolist():
        names = []
        if probiosml_any_hits["Mucus_adhesion_capacity"].get(gid, 0) > 0:
            names.append("adhesion")
        if (
            probiosml_any_hits["Antimicrobial_activity"].get(gid, 0) > 0
            or bagel_gate.get(gid, 0) > 0
            or antismash_amr_count.get(gid, 0) > 0
        ):
            names.append("antimicrobial")
        if any(probiosml_any_hits[k].get(gid, 0) > 0 for k in [
            "Acid_and_bile_tolerance",
            "Beneficial_metabolites",
            "Carbohydrate_metabolism",
            "Antioxidant_enzymes",
            "Immunomodulatory_potential",
        ]):
            names.append("metabolic")
        cazy_val = 0.0
        if "n_cazy_families" in indexed_summary.columns and gid in indexed_summary.index:
            cazy_val = float(pd.to_numeric(indexed_summary.loc[gid, "n_cazy_families"], errors="coerce"))
        if cazy_val > 0:
            names.append("cazy")
        if (
            probiosml_any_hits["General_stress_response"].get(gid, 0) > 0
            or probiosml_any_hits["Defense_systems"].get(gid, 0) > 0
            or crispr_count.get(gid, 0) > 0
        ):
            names.append("stress")
        blocks_count.append(float(len(names)))
        block_names.append(";".join(names))
        
    blocks_series = pd.Series(blocks_count, index=genomes["genome_id"])

    # Compute raw memberships
    mu_benefit_raw = 1.0 / (1.0 + np.exp(-kb * (pure_ahp_utility - thetab)))
    mu_blocks_ok = np.tanh(blocks_series)
    mu_benefit = mu_benefit_raw * mu_blocks_ok

    mu_risk = 1.0 / (1.0 + np.exp(-kr * (burden_w - thetar)))
    
    # Pathogen membership
    mu_pathogen_proxy = (1.0 / (1.0 + np.exp(-kp * (pm - thetap)))) * np.tanh(tier0)
    mu_pathogen_raw = np.maximum(pathogen_reference_similarity, mu_pathogen_proxy)
    mu_pathogen_total = np.maximum(mu_pathogen_raw, taxon_guard)

    # 11. Fuzzy rule composition (Strict partition of unity: sum = 1.0)
    b_rejected = mu_pathogen_total
    b_safety_review = mu_risk * (1.0 - mu_pathogen_total)
    b_priority = mu_benefit * (1.0 - mu_risk) * (1.0 - mu_pathogen_total)
    b_insufficient = (1.0 - mu_benefit) * (1.0 - mu_risk) * (1.0 - mu_pathogen_total)

    # Float precision normalisation (partition of unity)
    b_sum = b_rejected + b_safety_review + b_priority + b_insufficient
    b_rejected = b_rejected / b_sum
    b_safety_review = b_safety_review / b_sum
    b_priority = b_priority / b_sum
    b_insufficient = b_insufficient / b_sum

    # Defuzzification: Max membership principle
    memberships = np.column_stack([b_priority, b_safety_review, b_rejected, b_insufficient])
    max_indices = np.argmax(memberships, axis=1)
    status_map = {
        0: PRIORITY,
        1: SAFETY_REVIEW,
        2: REJECTED,
        3: INSUFFICIENT
    }
    final_status = [status_map[idx] for idx in max_indices]

    # 12. Build the output DataFrame
    out_df = pd.DataFrame({"genome_id": genomes["genome_id"]})
    out_df["ProbioScore_Utility"] = pure_ahp_utility.to_numpy()
    out_df["FCE_Membership_Priority"] = b_priority.to_numpy()
    out_df["FCE_Membership_SafetyReview"] = b_safety_review.to_numpy()
    out_df["FCE_Membership_Rejected"] = b_rejected.to_numpy()
    out_df["FCE_Membership_Insufficient"] = b_insufficient.to_numpy()
    out_df["ProbioScore_Status"] = final_status
    
    # Diagnostic memberships
    out_df["mu_Benefit_Raw"] = mu_benefit_raw.to_numpy()
    out_df["mu_Blocks_OK"] = mu_blocks_ok.to_numpy()
    out_df["mu_Benefit"] = mu_benefit.to_numpy()
    out_df["mu_Risk"] = mu_risk.to_numpy()
    out_df["pathogen_profile_feature_count"] = pathogen_profile_feature_count.to_numpy()
    out_df["mu_Pathogen_ReferenceSimilarity_Raw"] = pathogen_reference_similarity_raw.to_numpy()
    out_df["mu_Pathogen_ReferenceFeatureGate"] = pathogen_reference_feature_gate.to_numpy()
    out_df["mu_Pathogen_ReferenceIdentityAnchor"] = pathogen_reference_identity_anchor.to_numpy()
    out_df["mu_Pathogen_ReferenceSimilarity"] = pathogen_reference_similarity.to_numpy()
    out_df["mu_Pathogen_Proxy"] = mu_pathogen_proxy.to_numpy()
    out_df["mu_Pathogen_Raw"] = mu_pathogen_raw.to_numpy()
    out_df["mu_Pathogen_Total"] = mu_pathogen_total.to_numpy()
    
    # Internal variables for transparency and reproducibility
    out_df["tier_weighted_noncomp_burden"] = burden_w.to_numpy()
    out_df["tier0_marker_count"] = tier0.to_numpy()
    out_df["tier1_marker_count"] = tier1.to_numpy()
    out_df["tier2_marker_count"] = tier2.to_numpy()
    out_df["pathogen_membership_score"] = pm.to_numpy()
    out_df["positive_evidence_block_count"] = blocks_series.to_numpy()
    out_df["positive_evidence_blocks"] = block_names
    out_df["known_high_risk_taxon_guard"] = taxon_guard.to_numpy()
    out_df["known_high_risk_taxon_guard_reason"] = guard_reasons
    
    # Subpillars for full traceability
    out_df["B1_score"] = final_subpillars["B1_score"].to_numpy()
    out_df["B2_score"] = final_subpillars["B2_score"].to_numpy()
    out_df["B3_score"] = final_subpillars["B3_score"].to_numpy()
    out_df["B4_score"] = final_subpillars["B4_score"].to_numpy()
    out_df["B5_score"] = final_subpillars["B5_score"].to_numpy()
    out_df["C1_score"] = final_subpillars["C1_score"].to_numpy()
    out_df["C2_score"] = final_subpillars["C2_score"].to_numpy()
    out_df["C3_score"] = final_subpillars["C3_score"].to_numpy()
    out_df["C4_score"] = final_subpillars["C4_score"].to_numpy()
    out_df["C5_score"] = final_subpillars["C5_score"].to_numpy()

    # Risk factors
    out_df["A1_risk"] = a1.to_numpy()
    out_df["A2_risk"] = a2.to_numpy()
    out_df["A3_risk"] = a3.to_numpy()
    out_df["A4_risk"] = a4.to_numpy()

    # Counts
    out_df["bagel5_bacteriocin_count"] = bagel_count.to_numpy()
    out_df["bagel5_bacteriocin_gate"] = bagel_gate.to_numpy()
    out_df["antismash_antimicrobial_bgc_count"] = antismash_amr_count.to_numpy()
    out_df["epssmash_product_count"] = eps_count.to_numpy()
    out_df["crispr_array_count"] = crispr_count.to_numpy()
    out_df["ta_systems_context_count"] = ta_count.to_numpy()
    out_df["antismash_stress_bgc_count"] = antismash_stress_count.to_numpy()

    # Optional metadata merging
    metadata = _load_metadata(metadata_xlsx)
    if metadata is not None:
        if "genome_id" in metadata.columns:
            out_df = out_df.merge(metadata, on="genome_id", how="left")
        elif "metadata_code" in metadata.columns:
            out_df["metadata_code"] = out_df["genome_id"].map(_dataset_code_from_genome_id)
            out_df = out_df.merge(metadata, on="metadata_code", how="left")

    return out_df, pd.DataFrame(audit)
