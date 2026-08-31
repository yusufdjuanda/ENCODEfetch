"""Tests for encodefetch.summary.

Covers the report as it actually works today (nextflow/MultiQC-style plain
report: no Chart.js, no cards - tables and bar-lists), including the fixes
made during the --summary redesign:

  - compute_summary: experiment/file counts, download sizes, control
    availability, replication stats (case AND control)
  - compute_summary: "per experiment" breakdown fields (Organism, Biosample,
    Lab, Award, Platform, Classification, Donor Sex, Target) are counted per
    unique CASE experiment, not per file row and not including controls
  - compute_summary: "per file" breakdown fields (File Format, Run Type,
    Output Type, Assembly, File Status) stay counted per manifest row
  - print_terminal_summary: Rich terminal rendering, with and without
    query_params
  - write_html_summary: query_params plumbing (Selection criteria table,
    value pretty-printing, skip-if-redundant-with-query logic), the Assay
    panel's per-experiment Case:/Control: split, the Per experiment / Per
    file grouping with shared-constant-facts attached to the right group,
    replication status wording (Replicated/Mixed/Single replicate), the
    control-coverage/estimated-download/release-span facts, HTML escaping,
    and the embedded logo
  - write_json_summary: valid JSON output, query recorded under "query"
  - Edge cases: empty DataFrames, no control experiment at all
"""

import json

import pandas as pd

from encodefetch.summary import (
    _format_bytes,
    _pretty_date,
    _pretty_query_value,
    _replication_status,
    compute_summary,
    print_terminal_summary,
    write_html_summary,
    write_json_summary,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _rich_manifest_df():
    """Realistic manifest DataFrame with cases, a control, paired FASTQs,
    and replication data.

    Deliberately shaped so experiment-level and file-level counting give
    DIFFERENT answers, so tests can catch a regression back to file-row
    counting for experiment-level fields:
      - ENCSRCASE1 contributes 2 file rows (paired FASTQ mates on separate
        rows) but is 1 experiment - Run Type ("paired-ended") should read 2
        (file-level), Lab ("ENCODE Lab A") should read 1 (experiment-level).
      - ENCSRCASE1 and ENCSRCASE2 are both "TF ChIP-seq" - the Assay panel
        should read 2 experiments, not 3 file rows (CASE1 has 2 rows).
    """
    return pd.DataFrame(
        [
            {
                "experiment_accession": "ENCSRCASE1",
                "is_control": False,
                "matched_control_experiments": "ENCSRCTRL1",
                "controlled_by_files": "ENCFFCTRL1",
                "file_accession": "ENCFF001R1",
                "file_accession_r2": "ENCFF001R2",
                "file_format": "fastq",
                "file_size": 1_000_000_000,
                "file_size_r2": 1_000_000_000,
                "file_status": "released",
                "assay_title": "TF ChIP-seq",
                "target_label": "CTCF",
                "organism": "Homo sapiens",
                "biosample_term_name": "K562",
                "lab": "ENCODE Lab A",
                "award": "ENCODE4",
                "assembly": "",
                "output_type": "reads",
                "run_type": "paired-ended",
                "classification": "cell line",
                "platform": "Illumina NovaSeq 6000",
                "donor_sex": "female",
                "donor_age": "30",
                "donor_age_units": "year",
                "donor_life_stage": "adult",
                "donor_ethnicity": "",
                "perturbed": False,
                "date_released": "2023-01-15",
                "bio_replicate_count": 2,
                "biological_replicates": "1",
                "technical_replicates": "1_1",
                "replication_type": "isogenic",
            },
            {
                "experiment_accession": "ENCSRCASE1",
                "is_control": False,
                "matched_control_experiments": "ENCSRCTRL1",
                "controlled_by_files": "ENCFFCTRL1",
                "file_accession": "ENCFF002R1",
                "file_accession_r2": "ENCFF002R2",
                "file_format": "fastq",
                "file_size": 1_200_000_000,
                "file_size_r2": 1_200_000_000,
                "file_status": "released",
                "assay_title": "TF ChIP-seq",
                "target_label": "CTCF",
                "organism": "Homo sapiens",
                "biosample_term_name": "K562",
                "lab": "ENCODE Lab A",
                "award": "ENCODE4",
                "assembly": "",
                "output_type": "reads",
                "run_type": "paired-ended",
                "classification": "cell line",
                "platform": "Illumina NovaSeq 6000",
                "donor_sex": "female",
                "donor_age": "30",
                "donor_age_units": "year",
                "donor_life_stage": "adult",
                "donor_ethnicity": "",
                "perturbed": False,
                "date_released": "2023-01-15",
                "bio_replicate_count": 2,
                "biological_replicates": "2",
                "technical_replicates": "2_1",
                "replication_type": "isogenic",
            },
            {
                "experiment_accession": "ENCSRCASE2",
                "is_control": False,
                "matched_control_experiments": "ENCSRCTRL_MISSING",
                "controlled_by_files": "",
                "file_accession": "ENCFF003",
                "file_format": "fastq",
                "file_size": 500_000_000,
                "file_status": "released",
                "assay_title": "TF ChIP-seq",
                "target_label": "MYC",
                "organism": "Homo sapiens",
                "biosample_term_name": "K562",
                "lab": "ENCODE Lab B",
                "award": "ENCODE4",
                "assembly": "",
                "output_type": "reads",
                "run_type": "single-ended",
                "classification": "cell line",
                "platform": "Illumina NovaSeq 6000",
                "donor_sex": "male",
                "donor_age": "30",
                "donor_age_units": "year",
                "donor_life_stage": "adult",
                "donor_ethnicity": "",
                "perturbed": True,
                "date_released": "2024-06-01",
                "bio_replicate_count": 1,
                "biological_replicates": "1",
                "technical_replicates": "1_1",
                "replication_type": "unreplicated",
            },
            {
                "experiment_accession": "ENCSRCASE3_NOCONTROL",
                "is_control": False,
                "matched_control_experiments": "",
                "controlled_by_files": "",
                "file_accession": "ENCFF004",
                "file_format": "fastq",
                "file_size": 300_000_000,
                "file_status": "released",
                "assay_title": "total RNA-seq",
                "target_label": "",
                "organism": "Homo sapiens",
                "biosample_term_name": "K562",
                "lab": "ENCODE Lab B",
                "award": "ENCODE4",
                "assembly": "",
                "output_type": "reads",
                "run_type": "single-ended",
                "classification": "cell line",
                "platform": "Illumina NovaSeq 6000",
                "donor_sex": "male",
                "donor_age": "30",
                "donor_age_units": "year",
                "donor_life_stage": "adult",
                "donor_ethnicity": "",
                "perturbed": False,
                "date_released": "2024-06-01",
                "bio_replicate_count": "",
                "biological_replicates": "",
                "technical_replicates": "",
                "replication_type": "",
            },
            {
                "experiment_accession": "ENCSRCTRL1",
                "is_control": True,
                "matched_control_experiments": "",
                "controlled_by_files": "",
                "file_accession": "ENCFFCTRL1",
                "file_format": "fastq",
                "file_size": 600_000_000,
                "file_status": "released",
                "assay_title": "Control ChIP-seq",
                "target_label": "",
                "organism": "Homo sapiens",
                "biosample_term_name": "K562",
                "lab": "ENCODE Lab A",
                "award": "ENCODE4",
                "assembly": "",
                "output_type": "reads",
                "run_type": "single-ended",
                "classification": "cell line",
                "platform": "Illumina NovaSeq 6000",
                "donor_sex": "female",
                "donor_age": "30",
                "donor_age_units": "year",
                "donor_life_stage": "adult",
                "donor_ethnicity": "",
                "perturbed": False,
                "control_type": "input library",
                "date_released": "2023-01-10",
                "bio_replicate_count": 1,
                "biological_replicates": "1",
                "technical_replicates": "1_1",
                "replication_type": "unreplicated",
            },
        ]
    )


def _single_experiment_no_control_df():
    """One case experiment, no control at all - exercises the ":only-child"
    layout path (a lone replication article / metadata group)."""
    return pd.DataFrame(
        [
            {
                "experiment_accession": "ENCSRSOLO",
                "is_control": False,
                "matched_control_experiments": "",
                "file_accession": "ENCFF900",
                "file_format": "fastq",
                "file_size": 100_000_000,
                "file_status": "released",
                "assay_title": "total RNA-seq",
                "target_label": "",
                "organism": "Homo sapiens",
                "biosample_term_name": "HepG2",
                "lab": "ENCODE Lab A",
                "award": "ENCODE4",
                "assembly": "",
                "output_type": "reads",
                "run_type": "single-ended",
                "classification": "cell line",
                "platform": "Illumina NovaSeq 6000",
                "donor_sex": "male",
                "date_released": "2024-03-01",
                "bio_replicate_count": 1,
                "biological_replicates": "1",
                "technical_replicates": "1_1",
                "replication_type": "unreplicated",
            }
        ]
    )


def _empty_df():
    return pd.DataFrame(
        columns=[
            "experiment_accession",
            "is_control",
            "matched_control_experiments",
            "file_accession",
            "file_format",
            "file_size",
        ]
    )


# ---------------------------------------------------------------------------
# compute_summary: counts, sizes, control availability
# ---------------------------------------------------------------------------


def test_compute_summary_experiment_and_accession_counts():
    stats = compute_summary(_rich_manifest_df())
    assert stats["n_case_exps"] == 3
    assert stats["n_control_exps"] == 1
    assert stats["case_accessions"] == [
        "ENCSRCASE1",
        "ENCSRCASE2",
        "ENCSRCASE3_NOCONTROL",
    ]
    assert stats["control_accessions"] == ["ENCSRCTRL1"]


def test_compute_summary_unique_file_counts():
    stats = compute_summary(_rich_manifest_df())
    # Cases: ENCFF001R1, ENCFF001R2, ENCFF002R1, ENCFF002R2, ENCFF003, ENCFF004 = 6
    assert stats["n_unique_case_files"] == 6
    assert stats["n_unique_control_files"] == 1
    assert stats["n_total_rows"] == 5


def test_compute_summary_shared_control_deduplication():
    """A duplicate control file row must not inflate unique control counts."""
    df = _rich_manifest_df()
    dup_ctrl = df.iloc[4:5].copy()
    combined = pd.concat([df, dup_ctrl], ignore_index=True)
    stats = compute_summary(combined)
    assert stats["n_control_exps"] == 1
    assert stats["n_unique_control_files"] == 1
    assert stats["n_total_rows"] == 6


def test_compute_summary_download_sizes():
    stats = compute_summary(_rich_manifest_df())
    assert stats["case_download_size"] == 5_200_000_000
    assert stats["control_download_size"] == 600_000_000
    assert stats["total_download_size"] == 5_800_000_000


def test_compute_summary_control_availability():
    stats = compute_summary(_rich_manifest_df())
    # ENCSRCASE1 -> ENCSRCTRL1 is present in the df -> "usable"
    assert stats["cases_with_ctrl"] == 1
    # ENCSRCASE2 -> ENCSRCTRL_MISSING is linked but not present in the df
    assert stats["cases_unusable_ctrl"] == 1
    # ENCSRCASE3_NOCONTROL has no link at all
    assert stats["cases_without_ctrl"] == 1
    # n_cases_with_control_link is present-mode: metadata link only, so it
    # counts CASE1 (linked+present) AND CASE2 (linked, even though the
    # control file row is missing) - both "have a linked control" per
    # --control-presence=present semantics.
    assert stats["n_cases_with_control_link"] == 2


def test_compute_summary_empty_df():
    stats = compute_summary(_empty_df())
    assert stats["n_total_rows"] == 0
    assert stats["n_case_exps"] == 0
    assert stats["n_unique_case_files"] == 0
    assert stats["total_download_size"] == 0.0
    assert stats["case_accessions"] == []


# ---------------------------------------------------------------------------
# compute_summary: replication stats, case AND control
# ---------------------------------------------------------------------------


def test_compute_summary_case_replication_stats():
    stats = compute_summary(_rich_manifest_df())
    case_rep = stats["case_replication"]
    assert case_rep["one_bio_rep"] == 1  # ENCSRCASE2
    assert case_rep["two_plus_bio_rep"] == 1  # ENCSRCASE1
    assert case_rep["missing_rep_metadata"] == 1  # ENCSRCASE3_NOCONTROL
    rep_types = dict(case_rep["replication_type_counts"])
    assert rep_types.get("isogenic") == 1
    assert rep_types.get("unreplicated") == 1
    assert case_rep["unique_bio_rep_groups"] == 3
    assert case_rep["unique_tech_rep_groups"] == 3
    assert case_rep["files_per_tech_rep"] == {"min": 1, "max": 1, "median": 1}


def test_compute_summary_control_replication_stats():
    stats = compute_summary(_rich_manifest_df())
    ctrl_rep = stats["control_replication"]
    assert ctrl_rep["one_bio_rep"] == 1
    assert ctrl_rep["two_plus_bio_rep"] == 0
    assert ctrl_rep["missing_rep_metadata"] == 0


def test_replication_status_control_reads_single_replicate_not_unreplicated():
    """ENCODE doesn't require controls to be replicated - the status word
    for a control with 0/1 two-plus-rep experiments must not read as a
    deficiency the way it would for a case."""
    stats = compute_summary(_rich_manifest_df())
    two_plus, total, status = _replication_status(
        stats["control_replication"], is_control=True
    )
    assert (two_plus, total) == (0, 1)
    assert status == "Single replicate"


def test_replication_status_case_reads_unreplicated_when_none_replicated():
    rep = {"one_bio_rep": 2, "two_plus_bio_rep": 0, "missing_rep_metadata": 0}
    _, _, status = _replication_status(rep, is_control=False)
    assert status == "Unreplicated"


def test_replication_status_reads_replicated_when_all_meet_threshold():
    rep = {"one_bio_rep": 0, "two_plus_bio_rep": 3, "missing_rep_metadata": 0}
    _, _, status = _replication_status(rep)
    assert status == "Replicated"


def test_replication_status_reads_mixed():
    rep = {"one_bio_rep": 1, "two_plus_bio_rep": 1, "missing_rep_metadata": 0}
    _, _, status = _replication_status(rep)
    assert status == "Mixed"


def test_replication_status_empty_group_has_no_status():
    rep = {"one_bio_rep": 0, "two_plus_bio_rep": 0, "missing_rep_metadata": 0}
    two_plus, total, status = _replication_status(rep)
    assert (two_plus, total, status) == (0, 0, "")


# ---------------------------------------------------------------------------
# compute_summary: per-experiment vs. per-file breakdown basis
# ---------------------------------------------------------------------------


def test_breakdown_experiment_level_fields_are_case_only_and_deduplicated():
    """Organism/Biosample/Lab/Award/Platform/Classification/Donor Sex/Target
    must count unique CASE experiments - not file rows, and not including
    the control (whose own organism/lab/etc. would otherwise pad the count,
    e.g. 3 total instead of the true count of experiments you're studying).
    """
    stats = compute_summary(_rich_manifest_df())
    # Lab varies across the 3 case experiments: Lab A (CASE1) vs Lab B
    # (CASE2, CASE3) - 2 unique experiments each, NOT 4 file rows for Lab A.
    lab = dict(stats["lab_counts"])
    assert lab == {"ENCODE Lab A": 1, "ENCODE Lab B": 2}
    # The control is Lab A too, but must not appear in this case-only count
    # (it would make Lab A read 2 instead of 1 if controls were included).
    # Donor Sex likewise: female (CASE1) vs male (CASE2, CASE3).
    sex = dict(stats["sex_counts"])
    assert sex == {"male": 2, "female": 1}
    # Target is case-only and deduplicated too.
    target = dict(stats["target_counts"])
    assert target == {"CTCF": 1, "MYC": 1}


def test_breakdown_file_level_fields_stay_per_row_including_control():
    """Run Type genuinely varies within one experiment (paired vs single
    end can differ file to file), so it must stay a raw row count across
    case AND control - not deduplicated by experiment. ENCSRCASE1 alone
    contributes 2 "paired-ended" rows (its two FASTQ mates on separate
    manifest rows) despite being 1 experiment."""
    stats = compute_summary(_rich_manifest_df())
    run_type = dict(stats["run_type_counts"])
    assert run_type == {"single-ended": 3, "paired-ended": 2}


def test_breakdown_donor_age_combines_value_and_units():
    """donor_age + donor_age_units are two separate manifest columns that
    should read as one fact (e.g. "30 year"), deduplicated per case
    experiment like the other donor fields."""
    stats = compute_summary(_rich_manifest_df())
    assert dict(stats["donor_age_counts"]) == {"30 year": 3}


def test_breakdown_donor_age_omits_dangling_unit_when_age_missing():
    df = _rich_manifest_df()
    df.loc[df["experiment_accession"] == "ENCSRCASE2", "donor_age"] = ""
    stats = compute_summary(df)
    # No "  year" (empty age + unit) entry - just the two experiments that
    # still have both a value and a unit.
    assert dict(stats["donor_age_counts"]) == {"30 year": 2}


def test_breakdown_perturbed_is_case_only_and_varies():
    stats = compute_summary(_rich_manifest_df())
    assert dict(stats["perturbed_counts"]) == {"False": 2, "True": 1}


def test_breakdown_control_type_is_computed_from_controls_only():
    """control_type only exists on control rows (e.g. "input library") - it
    has no case-side value, so it must come from the control experiment,
    not accidentally pick up empty strings from the case rows."""
    stats = compute_summary(_rich_manifest_df())
    assert dict(stats["control_type_counts"]) == {"input library": 1}


# ---------------------------------------------------------------------------
# print_terminal_summary
# ---------------------------------------------------------------------------


def test_print_terminal_summary_runs_without_error_default_console():
    """Exercises the ``console is None`` branch that creates its own Console."""
    print_terminal_summary(_rich_manifest_df())


def test_print_terminal_summary_runs_with_explicit_console():
    from rich.console import Console

    console = Console(file=None, force_terminal=True, width=120)
    print_terminal_summary(_rich_manifest_df(), console=console)


def test_print_terminal_summary_with_query_params():
    from rich.console import Console

    console = Console(file=None, force_terminal=True, width=120)
    print_terminal_summary(
        _rich_manifest_df(),
        console=console,
        query_params={"Accessions": "ENCSRCASE1,ENCSRCASE2", "Control presence": "any"},
    )


def test_print_terminal_summary_empty_df():
    print_terminal_summary(_empty_df())


# ---------------------------------------------------------------------------
# write_html_summary: baseline structure
# ---------------------------------------------------------------------------


def test_write_html_summary_creates_valid_file(tmp_path):
    out = tmp_path / "summary.html"
    result = write_html_summary(_rich_manifest_df(), out)
    assert result == out
    assert out.exists()
    content = out.read_text(encoding="utf-8")
    assert "<!DOCTYPE html>" in content
    assert "Manifest Summary" in content
    assert "Selection criteria" in content
    assert "Metadata index" in content
    assert "Replication" in content


def test_write_html_summary_is_dependency_free(tmp_path):
    """No external JS charting library - bioinformaticians open these on
    clusters/offline machines, so no chart is worth losing silently."""
    out = tmp_path / "summary.html"
    write_html_summary(_rich_manifest_df(), out)
    content = out.read_text(encoding="utf-8")
    assert "<script" not in content
    assert "chart.js" not in content.lower()


def test_write_html_summary_embeds_logo(tmp_path):
    out = tmp_path / "summary.html"
    write_html_summary(_rich_manifest_df(), out)
    content = out.read_text(encoding="utf-8")
    assert "data:image/png;base64," in content
    assert 'alt="ENCODEfetch"' in content


def test_logo_data_uri_falls_back_to_empty_when_asset_missing(monkeypatch):
    """If the packaged logo.png can't be read (e.g. an unusual install
    layout), loading must degrade gracefully instead of raising."""
    import encodefetch.summary as summary_mod

    monkeypatch.setattr(summary_mod, "_LOGO_DATA_URI", None)

    def _raise(*_args, **_kwargs):
        raise FileNotFoundError("logo.png not found")

    monkeypatch.setattr(summary_mod.importlib.resources, "files", _raise)
    assert summary_mod._logo_data_uri() == ""


def test_write_html_summary_falls_back_to_text_brand_without_logo(
    tmp_path, monkeypatch
):
    import encodefetch.summary as summary_mod

    monkeypatch.setattr(summary_mod, "_logo_data_uri", lambda: "")
    out = tmp_path / "summary.html"
    write_html_summary(_rich_manifest_df(), out)
    content = out.read_text(encoding="utf-8")
    assert '<div class="brand">ENCODEfetch</div>' in content
    assert "data:image/png;base64," not in content


def test_write_html_summary_escapes_xss(tmp_path):
    df = _rich_manifest_df()
    df.loc[0, "lab"] = 'Lab <script>alert("xss")</script>'
    out = tmp_path / "summary.html"
    write_html_summary(df, out)
    content = out.read_text(encoding="utf-8")
    assert "<script>alert" not in content
    assert "&lt;script&gt;alert" in content


def test_write_html_summary_empty_df(tmp_path):
    out = tmp_path / "summary.html"
    write_html_summary(_empty_df(), out)
    content = out.read_text(encoding="utf-8")
    assert "ENCODEfetch" in content
    assert "Manifest Summary" in content


# ---------------------------------------------------------------------------
# write_html_summary: Selection criteria (query_params)
# ---------------------------------------------------------------------------


def test_write_html_summary_renders_selection_criteria_table(tmp_path):
    out = tmp_path / "summary.html"
    write_html_summary(
        _rich_manifest_df(),
        out,
        query_params={
            "Accessions": "ENCSRCASE1,ENCSRCASE2",
            "File type": "fastq",
            "Status": "released",
            "Control presence": "present",
        },
    )
    content = out.read_text(encoding="utf-8")
    # Pretty-printed: file type uppercased, status/control-presence capitalized.
    assert '<th>Accessions</th><td class="mono">ENCSRCASE1,ENCSRCASE2</td>' in content
    assert "<th>File type</th><td>FASTQ</td>" in content
    assert "<th>Status</th><td>Released</td>" in content
    assert "<th>Control presence</th><td>Present</td>" in content


def test_write_html_summary_no_query_params_omits_selection_table_rows(tmp_path):
    out = tmp_path / "summary.html"
    write_html_summary(_rich_manifest_df(), out)
    content = out.read_text(encoding="utf-8")
    # The "Selection criteria" heading/table shell still renders, but empty.
    assert "Selection criteria" in content
    assert "<th>Accessions</th>" not in content


# ---------------------------------------------------------------------------
# write_html_summary: Assay panel (per-experiment, case/control split)
# ---------------------------------------------------------------------------


def test_write_html_summary_assay_panel_counts_experiments_not_files(tmp_path):
    """The regression this whole redesign step was chasing: Assay must
    report experiment counts, not manifest-row counts, and must label which
    side (case vs control) each value belongs to."""
    out = tmp_path / "summary.html"
    write_html_summary(_rich_manifest_df(), out)
    content = out.read_text(encoding="utf-8")
    # TF ChIP-seq = 2 experiments (CASE1, CASE2), NOT 3 file rows
    # (CASE1 alone has 2 rows).
    assert '<span>Case: TF ChIP-seq</span><span class="count">2</span>' in content
    assert '<span>Case: total RNA-seq</span><span class="count">1</span>' in content
    assert (
        '<span>Control: Control ChIP-seq</span><span class="count">1</span>' in content
    )


def test_write_html_summary_assay_hides_case_side_when_assay_title_queried(tmp_path):
    """When --assay-title was explicitly given, the case side of Assay is
    redundant with Selection criteria - but the control's assay is new
    information (it's always a categorically different value), so it stays."""
    out = tmp_path / "summary.html"
    write_html_summary(
        _rich_manifest_df(), out, query_params={"Assay title": "TF ChIP-seq"}
    )
    content = out.read_text(encoding="utf-8")
    assert "Case: TF ChIP-seq" not in content
    assert "Case: total RNA-seq" not in content
    assert "Control: Control ChIP-seq" in content


# ---------------------------------------------------------------------------
# write_html_summary: Per experiment / Per file grouping
# ---------------------------------------------------------------------------


def test_write_html_summary_groups_breakdowns_by_basis(tmp_path):
    out = tmp_path / "summary.html"
    write_html_summary(_rich_manifest_df(), out)
    content = out.read_text(encoding="utf-8")
    assert 'class="group-tag exp">Per experiment' in content
    assert 'class="group-tag file">Per file' in content
    # Varying experiment-level fields land under "Per experiment" as bar panels.
    exp_idx = content.index("Per experiment")
    file_idx = content.index("Per file")
    assert exp_idx < content.index("<h3>Lab</h3>") < file_idx
    assert exp_idx < content.index("<h3>Donor Sex</h3>") < file_idx
    # Varying file-level fields land under "Per file".
    assert file_idx < content.index("<h3>Run Type</h3>")
    # Perturbed varies (False:2, True:1) so it's a bar panel; Donor Age,
    # Donor Life Stage, and Control Type are constant in this fixture so
    # they show as plain facts, not panels (see the constant-fields test).
    assert exp_idx < content.index("<h3>Perturbed</h3>") < file_idx


def test_write_html_summary_donor_and_control_type_facts(tmp_path):
    """The manifest columns that were previously dropped on the floor:
    donor age (combined with its units), life stage, and the control's own
    control_type - all attached to "Per experiment" since they're constant
    facts about the sample in this fixture."""
    out = tmp_path / "summary.html"
    write_html_summary(_rich_manifest_df(), out)
    content = out.read_text(encoding="utf-8")
    exp_section = content[content.index("Per experiment") : content.index("Per file")]
    assert "<th>Donor Age</th><td>30 year</td>" in exp_section
    assert "<th>Donor Life Stage</th><td>adult</td>" in exp_section
    assert "<th>Control Type</th><td>input library</td>" in exp_section
    # Donor Ethnicity is empty in the fixture, so it shouldn't render at all.
    assert "Donor Ethnicity" not in content


def test_write_html_summary_perturbed_hidden_when_queried(tmp_path):
    out = tmp_path / "summary.html"
    write_html_summary(_rich_manifest_df(), out, query_params={"Perturbed": "false"})
    content = out.read_text(encoding="utf-8")
    # Selection criteria legitimately echoes the filter itself...
    assert "<th>Perturbed</th><td>False</td>" in content
    # ...but the now-redundant Metadata index panel is gone.
    metadata_index = content[content.index("Metadata index") :]
    assert "<h3>Perturbed</h3>" not in metadata_index


def test_write_html_summary_constant_fields_attach_to_their_own_group(tmp_path):
    """A single-valued field isn't dumped into one generic "shared" bucket -
    it's attached to whichever group (Per experiment / Per file) it
    actually belongs to, so the two bases never get remixed."""
    out = tmp_path / "summary.html"
    write_html_summary(_rich_manifest_df(), out)
    content = out.read_text(encoding="utf-8")
    assert "Shared across all results" not in content
    exp_idx = content.index("Per experiment")
    file_idx = content.index("Per file")
    # Organism/Award/Classification/Platform are constant + experiment-level.
    assert exp_idx < content.index("<th>Organism</th>") < file_idx
    assert exp_idx < content.index("<th>Platform</th>") < file_idx
    # File Format/Output Type/File Status are constant + file-level.
    assert file_idx < content.index("<th>File Format</th>")
    assert file_idx < content.index("<th>Output Type</th>")


def test_write_html_summary_skips_metadata_categories_matching_query(tmp_path):
    out = tmp_path / "summary.html"
    write_html_summary(
        _rich_manifest_df(),
        out,
        query_params={
            "Organism": "Homo sapiens",
            "File type": "fastq",
        },
    )
    content = out.read_text(encoding="utf-8")
    # Selection criteria legitimately echoes both filters as their own rows -
    # only the Metadata index section (everything after its <h2>) should
    # have dropped the now-redundant Organism/File Format categories.
    metadata_index = content[content.index("Metadata index") :]
    assert "<th>Organism</th>" not in metadata_index
    assert "<th>File Format</th>" not in metadata_index
    # Unrelated categories are unaffected.
    assert "<h3>Lab</h3>" in metadata_index


# ---------------------------------------------------------------------------
# write_html_summary: metrics row and facts table
# ---------------------------------------------------------------------------


def test_write_html_summary_metrics_row(tmp_path):
    out = tmp_path / "summary.html"
    write_html_summary(_rich_manifest_df(), out)
    content = out.read_text(encoding="utf-8")
    assert (
        '<div class="number">3</div><div class="label">Case experiments</div>'
        in content
    )
    assert (
        '<div class="number">1</div><div class="label">Control experiment</div>'
        in content
    )
    assert '<div class="number">6</div><div class="label">Case files</div>' in content
    assert (
        '<div class="number">1</div><div class="label">Control files</div>' in content
    )
    assert (
        '<div class="number">5</div><div class="label">Manifest rows</div>' in content
    )


def test_write_html_summary_facts_table(tmp_path):
    out = tmp_path / "summary.html"
    write_html_summary(_rich_manifest_df(), out)
    content = out.read_text(encoding="utf-8")
    assert "2 of 3 case experiments have a linked control" in content
    total = _format_bytes(5_800_000_000)
    case = _format_bytes(5_200_000_000)
    control = _format_bytes(600_000_000)
    assert (
        f'<span class="mono">{total}</span> total ({case} case, {control} control)'
        in content
    )
    # No em dash or en dash anywhere in the facts table.
    assert "—" not in content
    assert "–" not in content
    assert "12 Jan" not in content  # sanity: dates are actually formatted
    assert "10 Jan 2023" in content
    assert "1 Jun 2024" in content


def test_write_html_summary_no_control_experiment(tmp_path):
    out = tmp_path / "summary.html"
    write_html_summary(_single_experiment_no_control_df(), out)
    content = out.read_text(encoding="utf-8")
    assert (
        '<div class="number">0</div><div class="label">Control experiments</div>'
        in content
    )
    assert "0 of 1 case experiments have a linked control" in content
    # Only the Case replication article renders - no Control panel, and no
    # 0/0 division artifacts from the control side.
    assert "Case experiment<span" in content
    assert "Control experiment<span" not in content


# ---------------------------------------------------------------------------
# write_json_summary
# ---------------------------------------------------------------------------


def test_write_json_summary_creates_valid_json(tmp_path):
    out = tmp_path / "summary.json"
    result = write_json_summary(_rich_manifest_df(), out)
    assert result == out
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["n_case_exps"] == 3
    assert data["n_control_exps"] == 1
    assert data["n_unique_case_files"] == 6
    assert data["n_cases_with_control_link"] == 2
    assert "generated_at" in data
    assert "case_replication" in data
    assert "control_replication" in data


def test_write_json_summary_records_query_params(tmp_path):
    out = tmp_path / "summary.json"
    write_json_summary(
        _rich_manifest_df(), out, query_params={"Control presence": "any"}
    )
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["query"] == {"Control presence": "any"}


def test_write_json_summary_omits_query_key_when_not_given(tmp_path):
    out = tmp_path / "summary.json"
    write_json_summary(_rich_manifest_df(), out)
    data = json.loads(out.read_text(encoding="utf-8"))
    assert "query" not in data


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------


def test_format_bytes_scales_units():
    assert _format_bytes(500) == "500.0 B"
    assert _format_bytes(2048) == "2.0 KB"
    assert _format_bytes(5_800_000_000) == "5.4 GB"
    assert _format_bytes(1024**4 * 2) == "2.0 TB"
    assert _format_bytes(1024**5 * 3) == "3.0 PB"


def test_pretty_query_value_file_type():
    assert _pretty_query_value("File type", "any") == "Any"
    assert _pretty_query_value("File type", "fastq") == "FASTQ"
    assert _pretty_query_value("File type", "fastq,bam") == "FASTQ, BAM"


def test_pretty_query_value_capitalized_keys():
    assert _pretty_query_value("Status", "released") == "Released"
    assert _pretty_query_value("Control presence", "present") == "Present"


def test_pretty_query_value_passthrough_for_unknown_keys():
    assert _pretty_query_value("Accessions", "ENCSR644VYX") == "ENCSR644VYX"
    assert _pretty_query_value("Organism", "Homo sapiens") == "Homo sapiens"


def test_pretty_date_formats_iso_date():
    iso, display = _pretty_date("2023-01-10")
    assert iso == "2023-01-10"
    assert display == "10 Jan 2023"


def test_pretty_date_falls_back_on_unparseable_input():
    iso, display = _pretty_date("not-a-date")
    assert iso == display == "not-a-date"
