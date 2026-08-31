from encodefetch.core import (
    build_file_record,
    experiments_to_df,
    extract_accessions_from_paths,
    search_accessions,
)
from encodefetch.encode_client import EncodeNotFoundError


def test_extract_accessions_from_paths_handles_strings_dicts_and_dedupes():
    values = [
        "/files/ENCFF000AAA/",
        "ENCFF000BBB",
        {"accession": "ENCFF000CCC"},
        {"@id": "/files/ENCFF000AAA/"},
        {},
        "",
        None,
    ]

    assert extract_accessions_from_paths(values) == [
        "ENCFF000AAA",
        "ENCFF000BBB",
        "ENCFF000CCC",
    ]


def test_build_file_record_extracts_controlled_by_files():
    exp_json = {
        "accession": "ENCSR000CASE",
        "target": {"label": "BRD4"},
        "biosample_ontology": {"term_name": "K562"},
        "lab": {"title": "Lab"},
        "award": {"rfa": "ENCODE4"},
        "replicates": [{"library": {"biosample": {"donor": {}}}}],
    }
    file_json = {
        "accession": "ENCFF000CASE",
        "file_format": "fastq",
        "status": "released",
        "controlled_by": [
            "/files/ENCFF000CTRL1/",
            {"@id": "/files/ENCFF000CTRL2/"},
        ],
    }

    record = build_file_record(
        file_json,
        exp_json=exp_json,
        is_control=False,
        matched_controls="ENCSR000CTRL",
    )

    assert record["controlled_by_files"] == "ENCFF000CTRL1,ENCFF000CTRL2"


# ---------------------------------------------------------------------------
# Graceful handling of a nonexistent accession (regression test: this used
# to crash the whole run with a raw traceback - see EncodeNotFoundError in
# encode_client.py). A bad accession should be skipped with a warning, not
# take down accessions that DO exist alongside it.
# ---------------------------------------------------------------------------

def _fake_experiment(accession: str) -> dict:
    return {
        "accession": accession,
        "assay_title": "TF ChIP-seq",
        "target": {"label": "CTCF"},
        "biosample_ontology": {"term_name": "K562"},
        "lab": {"title": "Test Lab"},
        "award": {"rfa": "ENCODE4"},
        "replicates": [{"library": {"biosample": {"donor": {}}}}],
        "possible_controls": [],
        "files": [
            {
                "accession": "ENCFF000CASE",
                "file_format": "fastq",
                "status": "released",
            }
        ],
    }


def test_search_accessions_skips_not_found_and_keeps_valid(monkeypatch):
    def fake_fetch_experiment(acc, auth=None, embedded=True):
        if acc == "ENCSRBADACCESSION999":
            raise EncodeNotFoundError(f"/experiments/{acc}/")
        return _fake_experiment(acc)

    monkeypatch.setattr("encodefetch.core.fetch_experiment", fake_fetch_experiment)

    df, records = search_accessions(["ENCSR000GOOD", "ENCSRBADACCESSION999"])

    assert not df.empty
    assert set(df["experiment_accession"]) == {"ENCSR000GOOD"}


def test_search_accessions_all_not_found_returns_empty_df(monkeypatch):
    def always_not_found(acc, auth=None, embedded=True):
        raise EncodeNotFoundError(f"/experiments/{acc}/")

    monkeypatch.setattr("encodefetch.core.fetch_experiment", always_not_found)

    df, records = search_accessions(["ENCSRBAD1", "ENCSRBAD2"])

    assert df.empty
    assert records == []


def test_experiments_to_df_skips_case_experiment_not_found(monkeypatch):
    """A case experiment (not fetched ahead of time by search_accessions,
    e.g. one that came from a /search/ hit that later 404s on the detail
    fetch) must be skipped, not crash the whole batch."""
    def fake_fetch_experiment(acc, auth=None, embedded=True):
        if acc == "ENCSRBADACCESSION999":
            raise EncodeNotFoundError(f"/experiments/{acc}/")
        return _fake_experiment(acc)

    monkeypatch.setattr("encodefetch.core.fetch_experiment", fake_fetch_experiment)

    df, records = experiments_to_df([
        {"accession": "ENCSR000GOOD"},
        {"accession": "ENCSRBADACCESSION999"},
    ])

    assert set(df["experiment_accession"]) == {"ENCSR000GOOD"}
