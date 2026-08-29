from __future__ import annotations
import pandas as pd

def collapse_fastq_pairs(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse paired-ended FASTQ rows into single rows anchored on R1.
    If R2 exists but is archived, drop R2 and mark single_end="true".

    Adds columns:
      fastq_1, fastq_2, url_r2, md5sum_r2, file_size_r2, file_accession_r2, single_end
    """
    if "file_format" not in df.columns:
        return df.copy()

    fq_mask = df["file_format"].astype(str).str.lower().eq("fastq")
    other = df[~fq_mask].copy()
    fq = df[fq_mask].copy()
    if fq.empty:
        return df.copy()

    def _path_from_row(r):
        lp = str(r.get("local_path","") or "").strip()
        return lp if lp else str(r.get("url",""))

    def _is_paired(run_type): return "paired" in str(run_type).lower()
    def _pe(x):
        s = str(x).strip()
        return "1" if s == "1" else ("2" if s == "2" else "")
    by_acc = {acc: row for acc, row in fq.set_index("file_accession").iterrows()}

    fq["_paired"] = fq["run_type"].apply(_is_paired)
    fq["_pe"] = fq["paired_end"].apply(_pe)
    fq["_path"] = fq.apply(_path_from_row, axis=1)
    fq["_mate"] = fq.get("paired_accession","")

    anchors = pd.concat([fq[fq["_paired"] & fq["_pe"].eq("1")], fq[~fq["_paired"]]], ignore_index=True)

    out = []
    for _, r1 in anchors.iterrows():
        rec = r1.to_dict()
        rec["fastq_1"] = r1["_path"]
        r2 = None
        if r1["_paired"]:
            macc = r1["_mate"]
            if macc and macc in by_acc:
                r2 = by_acc[macc]
            else:
                cand = fq[(fq["experiment_accession"] == r1["experiment_accession"]) & (fq["_paired"]) & fq["_pe"].eq("2")]
                if "biological_replicates" in fq.columns:
                    cand = cand[cand["biological_replicates"] == r1.get("biological_replicates","")]
                if "technical_replicates" in fq.columns:
                    cand = cand[cand["technical_replicates"] == r1.get("technical_replicates","")]
                if not cand.empty:
                    backptr = cand.get("paired_accession","") == r1["file_accession"]
                    r2 = cand[backptr].iloc[0] if backptr.any() else cand.iloc[0]

        if r2 is not None and str(r2.get("file_status", r2.get("status",""))).lower() != "archived":
            rec["fastq_2"] = _path_from_row(r2)
            rec["url_r2"] = r2.get("url","")
            rec["md5sum_r2"] = r2.get("md5sum","")
            rec["file_size_r2"] = r2.get("file_size","")
            rec["file_accession_r2"] = r2.get("file_accession","") or getattr(r2, "name", "")
            rec["single_end"] = "false"
        else:
            rec["fastq_2"] = ""
            rec["url_r2"] = ""
            rec["md5sum_r2"] = ""
            rec["file_size_r2"] = ""
            rec["file_accession_r2"] = ""
            rec["single_end"] = "true"

        out.append(rec)

    collapsed = pd.DataFrame(out)
    for c in ["_paired","_pe","_path","_mate"]:
        if c in collapsed.columns:
            collapsed.drop(columns=c, inplace=True, errors="ignore")

    combined = pd.concat([other, collapsed], ignore_index=True, sort=False)
    sort_cols = [c for c in ["is_control","experiment_accession","biological_replicates","technical_replicates"] if c in combined.columns]
    if sort_cols:
        combined = combined.sort_values(sort_cols).reset_index(drop=True)
    return combined


# ---------------------------------------------------------------------------
# Control-presence filter
# ---------------------------------------------------------------------------

_VALID_CONTROL_PRESENCE_MODES = frozenset({"any", "present", "none"})


def filter_control_presence(df: pd.DataFrame, mode: str) -> pd.DataFrame:
    """Retain cases (and their matched control rows) based on whether ENCODE
    metadata records a control relationship for the experiment.

    This reflects only the experiment-level relationship reported by ENCODE
    (``matched_control_experiments``). It does not depend on which file
    types, assemblies, or statuses were requested, so the result is stable
    regardless of ``--file-type`` or other file filters applied earlier in
    the pipeline.

    Parameters
    ----------
    df:
        DataFrame produced by :func:`experiments_to_df` (or its wrappers),
        containing at least the columns ``is_control``,
        ``experiment_accession``, and ``matched_control_experiments``.
    mode:
        One of:

        - ``any``     – keep all cases; preserves current behaviour (no-op).
        - ``present`` – retain only cases whose ENCODE metadata contains at
          least one linked control experiment (``matched_control_experiments``
          is non-empty).
        - ``none``    – retain only cases with no linked control.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame that always includes the matched control rows for
        every retained case.

    Raises
    ------
    ValueError
        If ``mode`` is not one of the accepted values.
    """
    if mode not in _VALID_CONTROL_PRESENCE_MODES:
        raise ValueError(
            f"Invalid control-presence mode {mode!r}. "
            f"Choose from: {sorted(_VALID_CONTROL_PRESENCE_MODES)}."
        )

    if mode == "any":
        return df.copy()

    if df.empty:
        return df.copy()

    cases = df[~df["is_control"]].copy()
    controls = df[df["is_control"]].copy()

    def _nonempty(series_value) -> bool:
        """Return True when a string cell is non-empty."""
        return bool(str(series_value).strip())

    if mode == "present":
        keep_mask = cases["matched_control_experiments"].apply(_nonempty)

    elif mode == "none":
        keep_mask = ~cases["matched_control_experiments"].apply(_nonempty)

    retained_cases = cases[keep_mask]

    # Preserve control rows whose experiment is referenced by a retained case.
    retained_ctrl_exps: set[str] = set()
    for val in retained_cases["matched_control_experiments"]:
        for acc in str(val).split(","):
            acc = acc.strip()
            if acc:
                retained_ctrl_exps.add(acc)

    retained_controls = controls[
        controls["experiment_accession"].isin(retained_ctrl_exps)
    ]

    combined = pd.concat([retained_cases, retained_controls], ignore_index=True, sort=False)
    sort_cols = [
        c for c in ["is_control", "experiment_accession", "file_accession"]
        if c in combined.columns
    ]
    if sort_cols:
        combined = combined.sort_values(sort_cols).reset_index(drop=True)
    return combined
