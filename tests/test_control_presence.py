"""Tests for encodefetch.postprocess.filter_control_presence.

``--control-presence`` describes only whether ENCODE metadata records a
control relationship for an experiment (``matched_control_experiments``).
It must not depend on which control *files* happen to be present after
file-type/assembly/status filtering.

Scenarios covered:
  - No linked control (``matched_control_experiments`` is empty)
  - Linked control, regardless of whether matching file rows are present
  - Shared control used by multiple cases
  - Multiple controls for a single case
  - All three modes: any / present / none
"""

import pytest
import pandas as pd

from encodefetch.postprocess import filter_control_presence


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _row(
    *,
    experiment_accession: str,
    is_control: bool,
    matched_control_experiments: str = "",
    file_accession: str = "",
) -> dict:
    """Return a minimal manifest row for use in test DataFrames."""
    return {
        "experiment_accession": experiment_accession,
        "is_control": is_control,
        "matched_control_experiments": matched_control_experiments,
        "file_accession": file_accession,
    }


def _make_df(rows: list) -> pd.DataFrame:
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def df_no_control():
    """Single case experiment with no linked control at all."""
    return _make_df([
        _row(experiment_accession="ENCSRCASE1", is_control=False,
             matched_control_experiments="", file_accession="ENCFF001"),
    ])


@pytest.fixture()
def df_linked_with_files():
    """Case linked to a control experiment that also has file rows present."""
    return _make_df([
        _row(experiment_accession="ENCSRCASE2", is_control=False,
             matched_control_experiments="ENCSRCTRL2",
             file_accession="ENCFF002"),
        _row(experiment_accession="ENCSRCTRL2", is_control=True,
             matched_control_experiments="", file_accession="ENCFFCTRL2"),
    ])


@pytest.fixture()
def df_linked_no_control_files():
    """Case has an experiment-level control link (per ENCODE metadata) but
    no rows for that control experiment survive the current file filters -
    e.g. requesting BED files for a ChIP-seq control that only has FASTQ.
    'present' must still keep this case since the link itself is what
    matters, not whether a control file row happens to be present."""
    return _make_df([
        _row(experiment_accession="ENCSR644VYX", is_control=False,
             matched_control_experiments="ENCSR513BAE",
             file_accession="ENCFF003"),
        # ENCSR513BAE has no file rows in this DataFrame (e.g. filtered
        # out by --file-type bed).
    ])


@pytest.fixture()
def df_shared_control():
    """Two case experiments sharing a single control experiment."""
    return _make_df([
        _row(experiment_accession="ENCSRCASE5A", is_control=False,
             matched_control_experiments="ENCSRCTRL5",
             file_accession="ENCFF005A"),
        _row(experiment_accession="ENCSRCASE5B", is_control=False,
             matched_control_experiments="ENCSRCTRL5",
             file_accession="ENCFF005B"),
        _row(experiment_accession="ENCSRCTRL5", is_control=True,
             matched_control_experiments="", file_accession="ENCFFCTRL5"),
    ])


@pytest.fixture()
def df_multiple_controls():
    """Single case with two distinct control experiments."""
    return _make_df([
        _row(experiment_accession="ENCSRCASE6", is_control=False,
             matched_control_experiments="ENCSRCTRL6A,ENCSRCTRL6B",
             file_accession="ENCFF006"),
        _row(experiment_accession="ENCSRCTRL6A", is_control=True,
             matched_control_experiments="", file_accession="ENCFFCTRL6A"),
        _row(experiment_accession="ENCSRCTRL6B", is_control=True,
             matched_control_experiments="", file_accession="ENCFFCTRL6B"),
    ])


@pytest.fixture()
def df_mixed():
    """Mix of: case-with-control-and-files, case-with-no-control, and
    case-with-link-but-no-surviving-control-files."""
    return _make_df([
        # Case A: control link with file rows present
        _row(experiment_accession="ENCSRCASE_A", is_control=False,
             matched_control_experiments="ENCSRCTRL_A", file_accession="ENCFFA"),
        _row(experiment_accession="ENCSRCTRL_A", is_control=True,
             matched_control_experiments="", file_accession="ENCFFCTRLA"),
        # Case B: no control
        _row(experiment_accession="ENCSRCASE_B", is_control=False,
             matched_control_experiments="", file_accession="ENCFFB"),
        # Case C: control link, but no file rows for that control present
        _row(experiment_accession="ENCSRCASE_C", is_control=False,
             matched_control_experiments="ENCSRCTRL_MISSING", file_accession="ENCFFC"),
    ])


# ---------------------------------------------------------------------------
# Mode: any (no-op)
# ---------------------------------------------------------------------------

def test_any_returns_all_rows(df_mixed):
    result = filter_control_presence(df_mixed, "any")
    assert len(result) == len(df_mixed)


def test_any_is_equal_copy(df_linked_with_files):
    result = filter_control_presence(df_linked_with_files, "any")
    pd.testing.assert_frame_equal(
        result.reset_index(drop=True),
        df_linked_with_files.reset_index(drop=True),
        check_like=True,
    )


# ---------------------------------------------------------------------------
# Mode: present
# ---------------------------------------------------------------------------

def test_present_drops_cases_with_no_control(df_mixed):
    result = filter_control_presence(df_mixed, "present")
    case_accs = set(result[~result["is_control"]]["experiment_accession"])
    assert "ENCSRCASE_B" not in case_accs


def test_present_keeps_cases_with_experiment_link(df_mixed):
    result = filter_control_presence(df_mixed, "present")
    case_accs = set(result[~result["is_control"]]["experiment_accession"])
    # Both A (link with files) and C (link, no surviving control files)
    # satisfy the 'present' criterion - it only checks the metadata link.
    assert "ENCSRCASE_A" in case_accs
    assert "ENCSRCASE_C" in case_accs


def test_present_keeps_case_even_without_matching_control_files(df_linked_no_control_files):
    """Regression test for the ENCSR644VYX/ENCSR513BAE scenario: a case
    linked to a control that has no BED files must still be kept under
    'present' when --file-type bed removed the control's rows."""
    result = filter_control_presence(df_linked_no_control_files, "present")
    case_accs = set(result[~result["is_control"]]["experiment_accession"])
    assert "ENCSR644VYX" in case_accs


def test_present_includes_control_rows_for_retained_cases(df_linked_with_files):
    result = filter_control_presence(df_linked_with_files, "present")
    ctrl_accs = set(result[result["is_control"]]["experiment_accession"])
    assert "ENCSRCTRL2" in ctrl_accs


def test_present_no_control_df_returns_empty(df_no_control):
    result = filter_control_presence(df_no_control, "present")
    assert result.empty


# ---------------------------------------------------------------------------
# Mode: none
# ---------------------------------------------------------------------------

def test_none_keeps_only_cases_without_control(df_mixed):
    result = filter_control_presence(df_mixed, "none")
    case_accs = set(result[~result["is_control"]]["experiment_accession"])
    assert case_accs == {"ENCSRCASE_B"}


def test_none_drops_cases_with_any_control_link(df_linked_with_files):
    result = filter_control_presence(df_linked_with_files, "none")
    assert result[~result["is_control"]].empty


def test_none_has_no_control_rows(df_mixed):
    """When only no-control cases are retained there are no controls to include."""
    result = filter_control_presence(df_mixed, "none")
    assert result[result["is_control"]].empty


# ---------------------------------------------------------------------------
# Shared / multiple controls
# ---------------------------------------------------------------------------

def test_shared_control_rows_present_for_both_cases(df_shared_control):
    """A single control experiment referenced by two cases must appear once."""
    result = filter_control_presence(df_shared_control, "present")
    ctrl_rows = result[result["is_control"]]
    assert set(ctrl_rows["experiment_accession"]) == {"ENCSRCTRL5"}
    # exactly one control file row (not duplicated)
    assert len(ctrl_rows) == 1


def test_multiple_controls_both_present(df_multiple_controls):
    result = filter_control_presence(df_multiple_controls, "present")
    ctrl_accs = set(result[result["is_control"]]["experiment_accession"])
    assert ctrl_accs == {"ENCSRCTRL6A", "ENCSRCTRL6B"}


def test_multiple_controls_case_retained(df_multiple_controls):
    result = filter_control_presence(df_multiple_controls, "present")
    case_accs = set(result[~result["is_control"]]["experiment_accession"])
    assert "ENCSRCASE6" in case_accs


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

def test_invalid_mode_raises_value_error(df_no_control):
    with pytest.raises(ValueError, match="Invalid control-presence mode"):
        filter_control_presence(df_no_control, "bad_mode")


def test_empty_dataframe_all_modes():
    empty = pd.DataFrame(columns=[
        "experiment_accession", "is_control",
        "matched_control_experiments", "file_accession",
    ])
    for mode in ("any", "present", "none"):
        result = filter_control_presence(empty, mode)
        assert result.empty, f"Expected empty result for mode={mode!r}"
