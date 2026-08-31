"""Summary generation for ENCODEfetch manifest DataFrames.

Provides three outputs from the same underlying statistics:
  - A compact Rich table printed to the terminal.
  - A self-contained HTML dashboard saved alongside the manifest.
  - A machine-readable JSON file for automated workflows.
"""

from __future__ import annotations

import base64
import html as html_mod
import importlib.resources
import json
from datetime import datetime
from pathlib import Path
from statistics import median
from typing import Dict, List, Optional, Tuple

import pandas as pd
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

_LOGO_DATA_URI: Optional[str] = None


def _logo_data_uri() -> str:
    """Base64 data URI for the ENCODEfetch logo, embedded so the HTML report
    stays a single portable file.

    Reads ``encodefetch/assets/logo.png`` as a package resource, so it works
    for a real pip install too, not just a source checkout. In this repo
    that file is a symlink to ``docs/img/logo.png`` (single source of
    truth); ``setuptools`` resolves it to real bytes when building the
    wheel, so installed users get an actual file, not a dangling link.
    """
    global _LOGO_DATA_URI
    if _LOGO_DATA_URI is None:
        try:
            data = (
                importlib.resources.files("encodefetch")
                .joinpath("assets", "logo.png")
                .read_bytes()
            )
            _LOGO_DATA_URI = "data:image/png;base64," + base64.b64encode(data).decode(
                "ascii"
            )
        except (FileNotFoundError, ModuleNotFoundError, OSError):
            _LOGO_DATA_URI = ""
    return _LOGO_DATA_URI


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _safe_col(df: pd.DataFrame, col: str) -> pd.Series:
    """Return the column if it exists, otherwise an empty string series."""
    if col in df.columns:
        return df[col].astype(str).replace({"nan": "", "None": "", "<NA>": ""})
    return pd.Series([""] * len(df), index=df.index)


def _combine_cols(
    df: pd.DataFrame, col_a: str, col_b: str, sep: str = " "
) -> pd.Series:
    """Join two columns (e.g. ``donor_age`` + ``donor_age_units`` -> "120 day"),
    leaving the row empty when the first column is empty rather than
    producing a dangling unit like " day"."""
    a = _safe_col(df, col_a)
    b = _safe_col(df, col_b)
    combined = (a + sep + b).where(a != "", "")
    return combined.str.strip()


def _format_bytes(n: float) -> str:
    """Human-readable file size."""
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} PB"


def _value_counts(series: pd.Series, top_n: int = 10) -> List[Tuple[str, int]]:
    """Return top-N (value, count) pairs for non-empty values."""
    filtered = series[series.astype(str).str.strip() != ""]
    counts = filtered.value_counts().head(top_n)
    return [(str(k), int(v)) for k, v in counts.items()]


def _unique_file_count(df: pd.DataFrame) -> int:
    """Count unique file accessions including R2 mates."""
    accessions = set()
    if "file_accession" in df.columns:
        accessions.update(df["file_accession"].dropna().unique())
    if "file_accession_r2" in df.columns:
        r2 = df["file_accession_r2"].dropna()
        accessions.update(r2[r2.astype(str).str.strip() != ""].unique())
    accessions.discard("")
    return len(accessions)


def _download_size(df: pd.DataFrame) -> float:
    """Estimated download size from unique files, including R2."""
    total = pd.to_numeric(_safe_col(df, "file_size"), errors="coerce").dropna().sum()
    if "file_size_r2" in df.columns:
        total += pd.to_numeric(df["file_size_r2"], errors="coerce").dropna().sum()
    return float(total)


def _replication_stats(subset: pd.DataFrame) -> Dict:
    """Compute per-experiment replication statistics for a case or control subset."""
    if subset.empty or "experiment_accession" not in subset.columns:
        return {
            "one_bio_rep": 0,
            "two_plus_bio_rep": 0,
            "missing_rep_metadata": 0,
            "replication_type_counts": [],
            "unique_bio_rep_groups": 0,
            "unique_tech_rep_groups": 0,
            "files_per_tech_rep": {"min": 0, "max": 0, "median": 0},
        }

    exps = subset.drop_duplicates(subset="experiment_accession")

    # Bio replicate count per experiment
    bio_counts = pd.to_numeric(_safe_col(exps, "bio_replicate_count"), errors="coerce")
    one_bio = int((bio_counts == 1).sum())
    two_plus = int((bio_counts >= 2).sum())
    missing = int(bio_counts.isna().sum())

    # Replication type
    rep_type_counts = _value_counts(_safe_col(exps, "replication_type"))

    # Unique replicate groups across all files
    bio_reps = set()
    tech_reps = set()
    for _, row in subset.iterrows():
        exp = str(row.get("experiment_accession", ""))
        br = str(row.get("biological_replicates", "")).strip()
        tr = str(row.get("technical_replicates", "")).strip()
        if br:
            for b in br.split(","):
                b = b.strip()
                if b:
                    bio_reps.add((exp, b))
        if tr:
            for t in tr.split(","):
                t = t.strip()
                if t:
                    tech_reps.add((exp, t))

    # Files per technical replicate
    files_per_tr = {}
    for _, row in subset.iterrows():
        exp = str(row.get("experiment_accession", ""))
        tr = str(row.get("technical_replicates", "")).strip()
        if tr:
            for t in tr.split(","):
                t = t.strip()
                if t:
                    key = (exp, t)
                    files_per_tr[key] = files_per_tr.get(key, 0) + 1

    counts_list = list(files_per_tr.values()) if files_per_tr else [0]

    return {
        "one_bio_rep": one_bio,
        "two_plus_bio_rep": two_plus,
        "missing_rep_metadata": missing,
        "replication_type_counts": rep_type_counts,
        "unique_bio_rep_groups": len(bio_reps),
        "unique_tech_rep_groups": len(tech_reps),
        "files_per_tech_rep": {
            "min": min(counts_list),
            "max": max(counts_list),
            "median": round(median(counts_list), 1),
        },
    }


# ---------------------------------------------------------------------------
# Statistics extraction — single source of truth
# ---------------------------------------------------------------------------


def compute_summary(df: pd.DataFrame) -> Dict:
    """Compute all summary statistics from a manifest DataFrame.

    Returns a dict consumed by ``print_terminal_summary``,
    ``write_html_summary``, and ``write_json_summary``.
    """
    cases = df[~df["is_control"]] if "is_control" in df.columns else df
    controls = df[df["is_control"]] if "is_control" in df.columns else pd.DataFrame()

    # Unique experiment counts
    n_case_exps = (
        cases["experiment_accession"].nunique()
        if "experiment_accession" in cases.columns
        else 0
    )
    n_control_exps = (
        controls["experiment_accession"].nunique()
        if not controls.empty and "experiment_accession" in controls.columns
        else 0
    )

    # Accession lists (used to name experiments in the terminal summary when few enough)
    case_accessions = (
        sorted(cases["experiment_accession"].dropna().unique().tolist())
        if "experiment_accession" in cases.columns
        else []
    )
    control_accessions = (
        sorted(controls["experiment_accession"].dropna().unique().tolist())
        if not controls.empty and "experiment_accession" in controls.columns
        else []
    )

    # Unique file counts (deduplicated by accession)
    n_unique_case_files = _unique_file_count(cases)
    n_unique_control_files = _unique_file_count(controls)

    # Manifest row counts (may differ from unique file counts due to shared controls)
    n_total_rows = len(df)
    n_case_rows = len(cases)
    n_control_rows = len(controls)

    # Download size
    total_download_size = _download_size(df)
    case_download_size = _download_size(cases)
    control_download_size = _download_size(controls) if not controls.empty else 0.0

    # Control availability
    cases_with_ctrl = 0
    cases_without_ctrl = 0
    cases_unusable_ctrl = 0
    # Cases whose ENCODE metadata records a control link at all (matches the
    # --control-presence=present semantics: metadata link only, independent
    # of whether a control file row survived the current filters).
    n_cases_with_control_link = 0
    if "matched_control_experiments" in cases.columns and not cases.empty:
        ctrl_exps_in_df = (
            set(controls["experiment_accession"].unique())
            if not controls.empty
            else set()
        )
        per_exp = cases.drop_duplicates(subset="experiment_accession")
        for _, row in per_exp.iterrows():
            linked = str(row.get("matched_control_experiments", "")).strip()
            if not linked:
                cases_without_ctrl += 1
            else:
                n_cases_with_control_link += 1
                linked_accs = {a.strip() for a in linked.split(",") if a.strip()}
                if linked_accs & ctrl_exps_in_df:
                    cases_with_ctrl += 1
                else:
                    cases_unusable_ctrl += 1

    # Breakdowns. Fields that describe the *sample/experiment* (constant for
    # every file in that experiment - organism, biosample, lab, etc.) are
    # counted per unique CASE experiment, not per file row and not
    # including controls: a 20-file experiment shouldn't outweigh a 2-file
    # one just because it has more replicates, and a control's own organism/
    # biosample/lab shouldn't pad "your" experiment counts (the control is
    # already surfaced separately, in the Assay panel and the Replication
    # section). Fields that are genuinely properties of individual *files*
    # (format, output type, assembly, status - an experiment can
    # legitimately mix fastq+bam+bed, or released+archived files) stay
    # counted per row across case and control alike.
    exps_cases = (
        cases.drop_duplicates(subset="experiment_accession")
        if "experiment_accession" in cases.columns
        else cases
    )

    assay_counts = _value_counts(_safe_col(df, "assay_title"))
    file_format_counts = _value_counts(_safe_col(df, "file_format"))
    organism_counts = _value_counts(_safe_col(exps_cases, "organism"))
    biosample_counts = _value_counts(_safe_col(exps_cases, "biosample_term_name"))
    target_counts = _value_counts(_safe_col(exps_cases, "target_label"))
    lab_counts = _value_counts(_safe_col(exps_cases, "lab"))
    award_counts = _value_counts(_safe_col(exps_cases, "award"))
    assembly_counts = _value_counts(_safe_col(df, "assembly"))
    output_type_counts = _value_counts(_safe_col(df, "output_type"))
    classification_counts = _value_counts(_safe_col(exps_cases, "classification"))
    platform_counts = _value_counts(_safe_col(exps_cases, "platform"))
    status_counts = _value_counts(_safe_col(df, "file_status"))
    sex_counts = _value_counts(_safe_col(exps_cases, "donor_sex"))
    run_type_counts = _value_counts(_safe_col(df, "run_type"))
    donor_age_counts = _value_counts(
        _combine_cols(exps_cases, "donor_age", "donor_age_units")
    )
    donor_life_stage_counts = _value_counts(_safe_col(exps_cases, "donor_life_stage"))
    donor_ethnicity_counts = _value_counts(_safe_col(exps_cases, "donor_ethnicity"))
    perturbed_counts = _value_counts(_safe_col(exps_cases, "perturbed"))

    # Control type (e.g. "input library") only exists on control rows - it
    # has no case-side counterpart, so it's computed from controls alone
    # rather than forced into the case-only pattern above.
    exps_controls = (
        controls.drop_duplicates(subset="experiment_accession")
        if not controls.empty and "experiment_accession" in controls.columns
        else controls
    )
    control_type_counts = _value_counts(_safe_col(exps_controls, "control_type"))

    # Date range
    dates = _safe_col(df, "date_released")
    dates_valid = dates[dates.astype(str).str.strip() != ""]
    date_min = str(dates_valid.min()) if not dates_valid.empty else ""
    date_max = str(dates_valid.max()) if not dates_valid.empty else ""

    # Replication
    case_replication = _replication_stats(cases)
    control_replication = _replication_stats(controls)

    return {
        "n_case_exps": n_case_exps,
        "n_control_exps": n_control_exps,
        "case_accessions": case_accessions,
        "control_accessions": control_accessions,
        "n_unique_case_files": n_unique_case_files,
        "n_unique_control_files": n_unique_control_files,
        "n_total_rows": n_total_rows,
        "n_case_rows": n_case_rows,
        "n_control_rows": n_control_rows,
        "total_download_size": total_download_size,
        "case_download_size": case_download_size,
        "control_download_size": control_download_size,
        "cases_with_ctrl": cases_with_ctrl,
        "cases_without_ctrl": cases_without_ctrl,
        "cases_unusable_ctrl": cases_unusable_ctrl,
        "n_cases_with_control_link": n_cases_with_control_link,
        "assay_counts": assay_counts,
        "file_format_counts": file_format_counts,
        "organism_counts": organism_counts,
        "biosample_counts": biosample_counts,
        "target_counts": target_counts,
        "lab_counts": lab_counts,
        "award_counts": award_counts,
        "assembly_counts": assembly_counts,
        "output_type_counts": output_type_counts,
        "classification_counts": classification_counts,
        "platform_counts": platform_counts,
        "status_counts": status_counts,
        "sex_counts": sex_counts,
        "run_type_counts": run_type_counts,
        "donor_age_counts": donor_age_counts,
        "donor_life_stage_counts": donor_life_stage_counts,
        "donor_ethnicity_counts": donor_ethnicity_counts,
        "perturbed_counts": perturbed_counts,
        "control_type_counts": control_type_counts,
        "date_min": date_min,
        "date_max": date_max,
        "case_replication": case_replication,
        "control_replication": control_replication,
    }


# ---------------------------------------------------------------------------
# Terminal summary (Rich)
# ---------------------------------------------------------------------------


def print_terminal_summary(
    df: pd.DataFrame,
    console: Optional[Console] = None,
    query_params: Optional[Dict[str, str]] = None,
) -> None:
    """Print a compact Rich summary to the terminal.

    Parameters
    ----------
    df:
        The manifest DataFrame.
    console:
        Rich console to print to (defaults to a new one).
    query_params:
        Optional ordered mapping of the CLI filters that produced ``df``
        (e.g. ``{"Accessions": "...", "File type": "fastq"}``), rendered at
        the top of the panel so the summary is traceable back to the
        command that generated it.
    """
    if console is None:
        console = Console()

    stats = compute_summary(df)

    # Overview panel
    overview = Table.grid(padding=(0, 2))
    overview.add_column(style="bold cyan", justify="right")
    overview.add_column()

    if query_params:
        for label, value in query_params.items():
            overview.add_row(label, value)
        overview.add_row("", "")

    overview.add_row("Case experiments", str(stats["n_case_exps"]))
    overview.add_row("Control experiments", str(stats["n_control_exps"]))
    if stats["n_case_exps"] > 0:
        overview.add_row(
            "Cases with control",
            f"{stats['n_cases_with_control_link']} / {stats['n_case_exps']}",
        )
    overview.add_row("Unique case files", str(stats["n_unique_case_files"]))
    overview.add_row("Unique control files", str(stats["n_unique_control_files"]))
    overview.add_row("Manifest rows", str(stats["n_total_rows"]))
    if stats["total_download_size"] > 0:
        overview.add_row(
            "Estimated download",
            f"{_format_bytes(stats['total_download_size'])} "
            f"(Case: {_format_bytes(stats['case_download_size'])}, "
            f"Control: {_format_bytes(stats['control_download_size'])})",
        )
    if stats["date_min"]:
        overview.add_row(
            "Release date range", f"{stats['date_min']} → {stats['date_max']}"
        )

    console.print(
        Panel(overview, title="[bold]ENCODEfetch Summary[/bold]", border_style="cyan")
    )


# ---------------------------------------------------------------------------
# JSON summary
# ---------------------------------------------------------------------------


def write_json_summary(
    df: pd.DataFrame, out_path: Path, query_params: Optional[Dict[str, str]] = None
) -> Path:
    """Write a machine-readable JSON summary alongside the manifest.

    Parameters
    ----------
    df:
        The manifest DataFrame.
    out_path:
        Destination file path (e.g. ``outdir / "summary.json"``).
    query_params:
        Optional mapping of the CLI filters that produced ``df``, recorded
        under the ``"query"`` key so the summary is traceable back to the
        command that generated it.

    Returns
    -------
    Path
        The path that was written.
    """
    stats = compute_summary(df)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Convert tuples to dicts for JSON
    serializable = {}
    for k, v in stats.items():
        if isinstance(v, list) and v and isinstance(v[0], tuple):
            serializable[k] = {label: count for label, count in v}
        else:
            serializable[k] = v

    if query_params:
        serializable["query"] = dict(query_params)
    serializable["generated_at"] = datetime.now().isoformat()
    out_path.write_text(json.dumps(serializable, indent=2), encoding="utf-8")
    return out_path


# ---------------------------------------------------------------------------
# HTML summary
# ---------------------------------------------------------------------------

_HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>ENCODEfetch &middot; Manifest Summary</title>
<style>
:root {{
  --paper:#fff;
  --ink:#172033;
  --muted:#667085;
  --line:#d9dee7;
  --line-strong:#aeb7c5;
  --wash:#f5f7f9;
  --navy:#0c2e62;
  --blue:#1679d3;
  --blue-tint:#e5f0fb;
  --teal:#0f8b8d;
  --teal-tint:#e3f5f4;
}}
*{{box-sizing:border-box}}
html{{color-scheme:light}}
body{{
  margin:0;
  background:var(--paper);
  color:var(--ink);
  font-family:Inter,ui-sans-serif,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
  font-size:16px;
  line-height:1.5;
  -webkit-font-smoothing:antialiased;
}}
main{{width:min(100% - 40px,1080px);margin:0 auto;padding:32px 0 48px}}
.mono,.number,.count{{font-family:"SFMono-Regular",Consolas,"Liberation Mono",monospace;font-variant-numeric:tabular-nums}}
.brand-bar{{background:#fff;border-bottom:1px solid var(--line)}}
.brand-bar-inner{{width:min(100% - 40px,1080px);height:56px;margin:0 auto;display:flex;align-items:center}}
.brand{{color:var(--navy);font-size:19px;font-weight:700;letter-spacing:-.02em}}
.brand-logo{{display:block;height:28px;width:auto}}
.hero{{background:linear-gradient(105deg,#092653 0%,#0d3977 57%,#087ec7 100%);color:#fff}}
.hero-inner{{width:min(100% - 40px,1080px);max-width:800px;margin:0 auto;padding:40px 0 36px}}
.hero h1{{margin:0;color:#fff;font-size:clamp(30px,4.6vw,44px);font-weight:700;letter-spacing:-.03em;line-height:1.05}}
.hero .generated{{margin-top:14px;color:#cfe0f7;font-size:13px}}
.hero .generated .mono{{margin-left:5px;color:#fff}}
.query{{padding:17px 19px;border:1px solid var(--line);border-radius:10px;background:#fff}}
.query h2{{margin:0 0 12px;color:var(--muted);font-size:13px;font-weight:700;letter-spacing:.08em;text-transform:uppercase}}
.eyebrow{{display:flex;align-items:center;gap:9px;margin:0 0 16px;color:var(--navy);font-size:13px;font-weight:700;letter-spacing:.08em;text-transform:uppercase}}
.eyebrow::before{{content:"";width:4px;height:16px;border-radius:2px;background:linear-gradient(180deg,var(--navy),var(--blue));flex-shrink:0}}
.group-tag{{display:inline-flex;align-items:center;margin-bottom:12px;padding:4px 12px;border-radius:9999px;font-size:11px;font-weight:700;letter-spacing:.06em;text-transform:uppercase}}
.group-tag.shared{{background:var(--wash);color:var(--muted)}}
.group-tag.exp{{background:var(--blue-tint);color:var(--blue)}}
.group-tag.file{{background:var(--teal-tint);color:var(--teal)}}
table.kv{{width:100%;border-collapse:collapse}}
table.kv tr+tr{{border-top:1px solid var(--line)}}
table.kv th,table.kv td{{padding:9px 0;text-align:left;font-size:15px;font-weight:400}}
table.kv th{{width:34%;color:var(--muted)}}
table.kv td{{font-family:inherit}}
.metrics{{display:grid;grid-template-columns:repeat(5,1fr);margin-top:18px;border:1px solid var(--line);border-radius:12px;overflow:hidden}}
.metric{{padding:22px}}
.metric:not(:first-child){{border-left:1px solid var(--line)}}
.metric .number{{color:var(--navy);font-size:36px;font-weight:600;line-height:1;letter-spacing:-.04em}}
.metric .label{{margin-top:9px;color:var(--muted);font-size:14px}}
.facts{{width:100%;margin-top:28px;border-collapse:collapse}}
.facts tr+tr{{border-top:1px solid var(--line)}}
.facts th,.facts td{{padding:11px 0;text-align:left;font-size:15px;font-weight:400}}
.facts th{{width:28%;color:var(--muted)}}
.facts td{{font-family:inherit}}
section{{margin-top:44px}}
.section-note{{max-width:680px;margin:-8px 0 18px;color:var(--muted);font-size:15px}}
.replication{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));border:1px solid var(--line);border-radius:12px;overflow:hidden}}
.replication article+article{{border-left:1px solid var(--line)}}
.replication article:only-child{{grid-column:1 / -1}}
.replication h3{{display:flex;justify-content:space-between;align-items:center;margin:0;padding:14px 16px;background:var(--wash);border-bottom:1px solid var(--line);font-size:15px}}
.status{{color:var(--muted);font-size:12px;font-weight:700;letter-spacing:.05em;text-transform:uppercase}}
.replication table{{width:100%;border-collapse:collapse}}
.replication tr+tr{{border-top:1px solid var(--line)}}
.replication th,.replication td{{padding:10px 16px;font-size:14px;font-weight:400}}
.replication th{{color:var(--muted);text-align:left}}
.replication td{{text-align:right}}
.metadata{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));border-top:1px solid var(--line-strong);margin-bottom:24px}}
table.kv + .metadata{{margin-top:16px}}
.metadata:last-child{{margin-bottom:0}}
.metadata article{{padding:17px 0 19px;border-bottom:1px solid var(--line)}}
.metadata article:nth-child(odd){{padding-right:26px}}
.metadata article:nth-child(even){{padding-left:26px;border-left:1px solid var(--line)}}
.metadata article:only-child{{grid-column:1 / -1;padding-right:0}}
.metadata h3{{margin:0 0 12px;color:var(--muted);font-size:12px;font-weight:700;letter-spacing:.07em;text-transform:uppercase}}
.bar{{margin-top:10px}}
.bar:first-of-type{{margin-top:0}}
.bar-label{{display:flex;justify-content:space-between;gap:10px;margin-bottom:6px;font-size:14px}}
.bar-label span:first-child{{overflow:hidden;text-overflow:ellipsis;white-space:nowrap}}
.count{{color:var(--muted)}}
.track{{height:3px;background:#e8ecf1}}
.fill{{display:block;height:100%;background:var(--blue)}}
footer{{display:flex;justify-content:space-between;gap:16px;flex-wrap:wrap;margin-top:48px;padding-top:16px;border-top:1px solid var(--line-strong);color:var(--muted);font-size:13px}}
footer a{{color:var(--navy);text-decoration:none}}
footer a:hover{{text-decoration:underline}}
@media(max-width:700px){{
  main{{width:min(100% - 28px,1080px)}}
  .brand-bar-inner{{width:min(100% - 28px,1080px)}}
  .hero-inner{{width:min(100% - 28px,1080px);padding:28px 0 24px}}
  .metrics{{grid-template-columns:1fr}}
  .metric:not(:first-child){{border-left:0;border-top:1px solid var(--line)}}
  .replication,.metadata{{grid-template-columns:1fr}}
  .replication article+article{{border-left:0;border-top:1px solid var(--line)}}
  .metadata article:nth-child(odd),.metadata article:nth-child(even){{padding:16px 0 18px;border-left:0}}
  .facts th,table.kv th{{width:42%}}
}}
@media print{{
  main{{width:100%;padding:0}}
  .hero{{background:var(--navy)}}
  section{{break-inside:avoid}}
  a{{color:inherit!important}}
}}
</style>
</head>
<body>
<div class="brand-bar">
  <div class="brand-bar-inner">
    {logo_html}
  </div>
</div>
<div class="hero">
  <div class="hero-inner">
    <h1>Manifest Summary</h1>
    <div class="generated">Generated<span class="mono">{timestamp}</span></div>
  </div>
</div>
<main>

  <aside class="query" aria-labelledby="query-title">
    <h2 id="query-title">Selection criteria</h2>
    <table class="kv">
      <tbody>
        {query_items}
      </tbody>
    </table>
  </aside>

  <section class="metrics" aria-label="Report totals">
    {metrics}
  </section>

  <table class="facts" aria-label="Report details">
    <tbody>
      {facts_rows}
    </tbody>
  </table>

  {replication_section}

  {breakdown_section}

  <footer>
    <span>ENCODEfetch manifest report</span>
    <a href="https://github.com/khan-lab/ENCODEfetch">github.com/khan-lab/ENCODEfetch</a>
  </footer>
</main>
</body>
</html>
"""

_PRETTY_UPPER_KEYS = {"File type"}
_PRETTY_CAP_KEYS = {"Status", "Control presence", "Assembly", "Perturbed"}


def _pretty_query_value(key: str, value: str) -> str:
    """Title-case CLI filter values that were passed through in raw lowercase."""
    if key in _PRETTY_UPPER_KEYS:
        if value == "any":
            return "Any"
        return ", ".join(v.strip().upper() for v in value.split(","))
    if key in _PRETTY_CAP_KEYS:
        return value.capitalize()
    return value


def _kv_row_html(label: str, value_html: str, mono_value: bool = False) -> str:
    """One aligned key/value row for a ``table.kv``/``table.facts``: the key
    column is a fixed width so every value in the table starts at the same
    x position, however long the keys are."""
    cls = ' class="mono"' if mono_value else ""
    return f"<tr><th>{html_mod.escape(label)}</th><td{cls}>{value_html}</td></tr>"


def _query_items_html(query_params: Dict[str, str]) -> str:
    rows = [
        _kv_row_html(
            k,
            html_mod.escape(_pretty_query_value(k, str(v))),
            mono_value=(k == "Accessions"),
        )
        for k, v in query_params.items()
    ]
    return "\n        ".join(rows)


def _metric_html(value, label: str) -> str:
    return f'<div class="metric"><div class="number">{value}</div><div class="label">{html_mod.escape(label)}</div></div>'


def _pretty_date(raw: str) -> Tuple[str, str]:
    """Return (iso, display) for a date string, falling back to the raw value."""
    for fmt in ("%Y-%m-%d",):
        try:
            dt = datetime.strptime(raw, fmt)
            return raw, dt.strftime("%-d %b %Y")
        except ValueError:
            continue
    return raw, raw


def _bar_html(label: str, count: int, max_count: int) -> str:
    pct = (count / max_count * 100) if max_count else 0
    return (
        f'<div class="bar"><div class="bar-label">'
        f'<span>{html_mod.escape(label or "(empty)")}</span><span class="count">{count}</span></div>'
        f'<div class="track"><span class="fill" style="width:{pct:.1f}%"></span></div></div>'
    )


def _metadata_article_html(title: str, items: List[Tuple[str, int]]) -> str:
    max_count = max((c for _, c in items), default=0)
    bars = "".join(_bar_html(val, count, max_count) for val, count in items)
    return f"<article><h3>{html_mod.escape(title)}</h3>{bars}</article>"


def _replication_status(rep: Dict, is_control: bool = False) -> Tuple[int, int, str]:
    """Return (two_plus, total, status_text) for a replication stat dict.

    ``is_control`` softens the wording for the "none replicated" case: ENCODE
    only requires 2+ biological replicates for the case (target) experiment,
    since that's what needs reproducibility QC. A control/input experiment
    is routinely a single, deeply-sequenced sample - often reused across
    several case experiments - so calling it "Unreplicated" would read as a
    defect it isn't; "Single replicate" states the same fact neutrally.
    """
    total = rep["one_bio_rep"] + rep["two_plus_bio_rep"] + rep["missing_rep_metadata"]
    two_plus = rep["two_plus_bio_rep"]
    if total == 0:
        status = ""
    elif two_plus == total:
        status = "Replicated"
    elif two_plus > 0:
        status = "Mixed"
    else:
        status = "Single replicate" if is_control else "Unreplicated"
    return two_plus, total, status


def _replication_article_html(title: str, rep: Dict, is_control: bool = False) -> str:
    two_plus, total, status = _replication_status(rep, is_control=is_control)
    if total == 0:
        return ""
    rows = [
        ("2+ biological replicates", f"{two_plus} / {total}"),
        ("Biological replicate groups", rep["unique_bio_rep_groups"]),
        ("Technical replicate groups", rep["unique_tech_rep_groups"]),
    ]
    rows_html = "".join(
        f'<tr><th>{html_mod.escape(label)}</th><td class="mono">{v}</td></tr>'
        for label, v in rows
    )
    return (
        f'<article><h3>{html_mod.escape(title)}<span class="status">{status}</span></h3>'
        f"<table><tbody>{rows_html}</tbody></table></article>"
    )


def write_html_summary(
    df: pd.DataFrame, out_path: Path, query_params: Optional[Dict[str, str]] = None
) -> Path:
    """Write a self-contained HTML report summarising the manifest.

    Single-page, static, dependency-free (no external JS or fonts —
    bioinformaticians frequently open these reports on clusters or offline
    machines, so every visualization is plain CSS/HTML with system fonts).
    A plain, print-friendly report layout: a "Selection criteria" list, a
    totals row, a facts table, then Replication and Metadata index sections
    as tables/bar lists - not cards.

    Parameters
    ----------
    df:
        The manifest DataFrame.
    out_path:
        Destination file path (e.g. ``outdir / "summary.html"``).
    query_params:
        Optional mapping of the CLI filters that produced ``df``, rendered
        under "Selection criteria" near the top of the page. Metadata index
        categories that would just re-state one of these filters (e.g. an
        Assay panel when ``--assay-title`` was explicitly given) are
        skipped, since selection criteria already shows it; categories that
        are genuine properties of the returned metadata (not something you
        searched for) are always kept, even when every row happens to share
        one value.

    Returns
    -------
    Path
        The path that was written.
    """
    stats = compute_summary(df)
    query_params = query_params or {}

    n_case = stats["n_case_exps"]
    n_control = stats["n_control_exps"]

    query_items = _query_items_html(query_params)

    # Metrics: experiment counts and file counts, each split case/control, so
    # "how many experiments" and "how many files" are both stated plainly and
    # neither is mistaken for the other (several files from one experiment,
    # e.g. replicate FASTQs sharing a biosample, isn't several experiments).
    metrics = "\n    ".join(
        [
            _metric_html(
                n_case, "Case experiments" if n_case != 1 else "Case experiment"
            ),
            _metric_html(
                n_control,
                "Control experiments" if n_control != 1 else "Control experiment",
            ),
            _metric_html(stats["n_unique_case_files"], "Case files"),
            _metric_html(stats["n_unique_control_files"], "Control files"),
            _metric_html(stats["n_total_rows"], "Manifest rows"),
        ]
    )

    # Facts table: control coverage sentence, download size, release span.
    facts_rows = []
    n_with_ctrl = stats["n_cases_with_control_link"]
    if n_case > 0:
        verb = "has" if n_with_ctrl == 1 else "have"
        facts_rows.append(
            _kv_row_html(
                "Control coverage",
                f"{n_with_ctrl} of {n_case} case experiments {verb} a linked control",
            )
        )
    if stats["total_download_size"] > 0:
        facts_rows.append(
            _kv_row_html(
                "Estimated download",
                f'<span class="mono">{_format_bytes(stats["total_download_size"])}</span> total '
                f'({_format_bytes(stats["case_download_size"])} case, '
                f'{_format_bytes(stats["control_download_size"])} control)',
            )
        )
    if stats["date_min"]:
        min_iso, min_disp = _pretty_date(stats["date_min"])
        max_iso, max_disp = _pretty_date(stats["date_max"])
        facts_rows.append(
            _kv_row_html(
                "Release span",
                f'<time datetime="{min_iso}">{min_disp}</time> to <time datetime="{max_iso}">{max_disp}</time>',
            )
        )

    # Replication section. ENCODE only expects 2+ biological replicates for
    # the case (target) experiment, where reproducibility QC matters - not
    # for the control, which is routinely a single deeply-sequenced sample
    # (often reused across several case experiments). So the note below only
    # holds the case side to that threshold, and calls out the control's
    # different standard rather than lumping both into one ratio.
    case_two_plus, case_total, _ = _replication_status(stats["case_replication"])
    ctrl_two_plus, ctrl_total, _ = _replication_status(
        stats["control_replication"], is_control=True
    )
    rep_articles = [
        html
        for html in (
            _replication_article_html("Case experiment", stats["case_replication"]),
            _replication_article_html(
                "Control experiment", stats["control_replication"], is_control=True
            ),
        )
        if html
    ]
    if case_total:
        rep_note = f"{case_two_plus} of {case_total} case experiment(s) meet the two-biological-replicate threshold."
    else:
        rep_note = f"{ctrl_two_plus} of {ctrl_total} control experiment(s) have 2+ biological replicates."
    replication_section = (
        f'<section aria-labelledby="replication-title">'
        f'<h2 class="eyebrow" id="replication-title">Replication</h2>'
        f'<p class="section-note">{rep_note}</p>'
        f'<div class="replication">{"".join(rep_articles)}</div>'
        f"</section>"
        if rep_articles
        else ""
    )

    # Metadata breakdown. Skip any category that just echoes a filter already
    # shown in the selection criteria above - it's a searched-for value, not
    # a discovered property of the metadata (e.g. --status always yields one
    # File Status value; that's not worth a bar).
    file_type_filter = query_params.get("File type", "any")
    # Assay is handled separately from the other breakdowns: case and control
    # experiments are always tagged with categorically different assay
    # titles (e.g. "TF ChIP-seq" vs "Control ChIP-seq"), and counting by
    # manifest row conflates "16 rows of TF ChIP-seq" with "16 TF ChIP-seq
    # experiments" - it's really 1 case experiment and 1 control experiment.
    # So this counts unique experiments, not file rows, and labels each
    # value by which side it's on.
    cases_df = df[~df["is_control"]] if "is_control" in df.columns else df
    controls_df = df[df["is_control"]] if "is_control" in df.columns else pd.DataFrame()
    case_assay_counts = (
        _value_counts(
            _safe_col(
                cases_df.drop_duplicates(subset="experiment_accession"), "assay_title"
            )
        )
        if not cases_df.empty
        else []
    )
    control_assay_counts = (
        _value_counts(
            _safe_col(
                controls_df.drop_duplicates(subset="experiment_accession"),
                "assay_title",
            )
        )
        if not controls_df.empty
        else []
    )
    assay_items = []
    if "Assay title" not in query_params:
        assay_items += [(f"Case: {val}", count) for val, count in case_assay_counts]
    assay_items += [(f"Control: {val}", count) for val, count in control_assay_counts]

    # "exp" fields are constant per experiment (counted per unique
    # experiment); "file" fields genuinely vary within one experiment
    # (counted per manifest row - e.g. one experiment can mix fastq+bam, or
    # single- and paired-end runs across replicates). Labeling the two
    # groups separately means a reader never has to wonder why a "file"
    # category's numbers don't sum to the experiment count - they were
    # never meant to.
    breakdown_configs = [
        ("Organism", stats["organism_counts"], "Organism" not in query_params, "exp"),
        (
            "Biosample",
            stats["biosample_counts"],
            "Biosample" not in query_params,
            "exp",
        ),
        ("Target", stats["target_counts"], "Target label" not in query_params, "exp"),
        ("Lab", stats["lab_counts"], True, "exp"),
        ("Award", stats["award_counts"], True, "exp"),
        ("Platform", stats["platform_counts"], True, "exp"),
        ("Donor Sex", stats["sex_counts"], True, "exp"),
        ("Donor Age", stats["donor_age_counts"], True, "exp"),
        ("Donor Life Stage", stats["donor_life_stage_counts"], True, "exp"),
        ("Donor Ethnicity", stats["donor_ethnicity_counts"], True, "exp"),
        (
            "Perturbed",
            stats["perturbed_counts"],
            "Perturbed" not in query_params,
            "exp",
        ),
        ("Classification", stats["classification_counts"], True, "exp"),
        ("Control Type", stats["control_type_counts"], True, "exp"),
        ("File Format", stats["file_format_counts"], file_type_filter == "any", "file"),
        ("Run Type", stats["run_type_counts"], True, "file"),
        ("Assembly", stats["assembly_counts"], "Assembly" not in query_params, "file"),
        ("Output Type", stats["output_type_counts"], True, "file"),
        ("File Status", stats["status_counts"], "Status" not in query_params, "file"),
    ]
    # A field with only one distinct value isn't a "breakdown" - it's just a
    # fact about the sample (e.g. donor sex, organism, biosample can't vary
    # within a single experiment, and often don't even across a handful of
    # related ones). A bar chart that's always 100% one bar tells you
    # nothing a plain fact doesn't already say, so those go in one compact
    # table instead of a whole panel; only fields that actually vary keep
    # the bar-list treatment. Assay always keeps the bar-list treatment
    # (even with 1 item per side) since the case/control labeling itself is
    # the point.
    kept = [
        (title, items, basis)
        for title, items, show, basis in breakdown_configs
        if items and show
    ]
    varying = [(title, items, basis) for title, items, basis in kept if len(items) > 1]
    # Constants split by basis too: a single-valued experiment field (e.g.
    # every experiment here happens to be female) is a "shared across
    # experiments" fact, not "shared across files" - conflating them in one
    # list would repeat the exact mixing problem the exp/file split above
    # already fixes for the varying fields.
    constant_exp = [
        (title, items[0][0])
        for title, items, basis in kept
        if len(items) == 1 and basis == "exp"
    ]
    constant_file = [
        (title, items[0][0])
        for title, items, basis in kept
        if len(items) == 1 and basis == "file"
    ]

    exp_articles = [_metadata_article_html(t, i) for t, i, b in varying if b == "exp"]
    if assay_items:
        exp_articles.insert(0, _metadata_article_html("Assay", assay_items))
    file_articles = [_metadata_article_html(t, i) for t, i, b in varying if b == "file"]

    def _group_html(
        label: str,
        kind: str,
        constant_items: List[Tuple[str, str]],
        articles: List[str],
    ) -> str:
        if not constant_items and not articles:
            return ""
        parts = [f'<h3 class="group-tag {kind}">{html_mod.escape(label)}</h3>']
        if constant_items:
            shared_rows = "\n        ".join(
                _kv_row_html(title, html_mod.escape(value))
                for title, value in constant_items
            )
            parts.append(f'<table class="kv">{shared_rows}</table>')
        if articles:
            parts.append(f'<div class="metadata">{"".join(articles)}</div>')
        return "".join(parts)

    exp_group_html = _group_html("Per experiment", "exp", constant_exp, exp_articles)
    file_group_html = _group_html("Per file", "file", constant_file, file_articles)
    breakdown_section = (
        f'<section aria-labelledby="metadata-title">'
        f'<h2 class="eyebrow" id="metadata-title">Metadata index</h2>'
        f"{exp_group_html}"
        f"{file_group_html}"
        f"</section>"
        if (exp_group_html or file_group_html)
        else ""
    )

    logo_uri = _logo_data_uri()
    logo_html = (
        f'<img class="brand-logo" src="{logo_uri}" alt="ENCODEfetch">'
        if logo_uri
        else '<div class="brand">ENCODEfetch</div>'
    )

    html_content = _HTML_TEMPLATE.format(
        logo_html=logo_html,
        timestamp=datetime.now().strftime("%Y-%m-%d %H:%M"),
        query_items=query_items,
        metrics=metrics,
        facts_rows="\n      ".join(facts_rows),
        replication_section=replication_section,
        breakdown_section=breakdown_section,
    )

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html_content, encoding="utf-8")
    return out_path
