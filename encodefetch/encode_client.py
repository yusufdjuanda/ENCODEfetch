from __future__ import annotations

from typing import List, Tuple, Optional, Dict, Any
import requests

ENCODE_BASE = "https://www.encodeproject.org"
HEADERS = {"accept": "application/json"}


class EncodeNotFoundError(Exception):
    """A specific ENCODE resource (e.g. an experiment accession) doesn't exist.

    Distinct from a search that legitimately matched zero results (which
    ``encode_get`` returns normally, not as an error) and from other HTTP
    failures (which still raise ``requests.HTTPError`` - server errors,
    auth failures, etc. are not "not found" and callers shouldn't silently
    skip them the way they can skip a bad accession).
    """

    def __init__(self, path: str):
        self.path = path
        super().__init__(f"Not found on ENCODE: {path}")


def encode_get(path_or_url: str,
               params: Optional[list | dict] = None,
               auth=None,
               timeout: int = 120,
               raw_query: Optional[str] = None) -> Dict[str, Any] | None:
    url = path_or_url if path_or_url.startswith("http") else (
        ENCODE_BASE.rstrip("/") + "/" + path_or_url.lstrip("/")
    )
    if params is None:
        params = []
    elif isinstance(params, dict):
        params = list(params.items())
    else:
        params = list(params)

    if not any(k == "format" for k, _ in params):
        params.append(("format", "json"))

    s = requests.Session()
    req = requests.Request("GET", url, headers=HEADERS, params=params, auth=auth)
    prepped = s.prepare_request(req)
    if raw_query:
        sep = '&' if '?' in prepped.url else '?'
        prepped.url = f"{prepped.url}{sep}{raw_query}"
    r = s.send(prepped, timeout=timeout)

    if r.status_code == 404:
        # ENCODE's /search/ endpoint returns HTTP 404 even for a perfectly
        # valid query that just matches zero experiments (confirmed: same
        # 404 + "No results found" body whether the filter is valid-but-
        # unmatched or nonsense - ENCODE doesn't distinguish). A response
        # shaped like a search result (has "@graph") is that case: return
        # it as-is so callers see an empty result, not a crash. Anything
        # else at 404 (e.g. /experiments/{accession}/ for an accession that
        # doesn't exist) is a genuine "not found".
        try:
            data = r.json()
        except ValueError:
            data = None
        if isinstance(data, dict) and "@graph" in data:
            return data
        raise EncodeNotFoundError(path_or_url)

    r.raise_for_status()
    return r.json()

def fetch_experiment(accession: str, auth=None, embedded: bool = True):
    params = {"format": "json"}
    if embedded:
        params["frame"] = "embedded"
    return encode_get(f"/experiments/{accession}/", params=params, auth=auth)

def build_params(
    assay_title: Optional[str] = None,
    target_labels: Optional[list[str]] = None,
    organism: Optional[str] = None,
    biosample: Optional[str] = None,
    status: str = "released",
    limit: str = "all",
    extra_params: Optional[dict] = None,
    perturbed: Optional[str] = None,
    series: Optional[str] = None,
) -> List[Tuple[str, str]]:
    """Build a list of (key, value) params (repeated keys preserved)."""
    p: List[Tuple[str, str]] = [("type","Experiment")]
    if assay_title:
        p.append(("assay_title", assay_title))
    if status:
        p.append(("status", status))
    if organism:
        p.append(("replicates.library.biosample.donor.organism.scientific_name", organism))
    if biosample:
        p.append(("biosample_ontology.term_name", biosample))
    if perturbed is not None:
        p.append(("perturbed", perturbed.lower()))
    if series:
        p.append(("related_series.@type", series))
    if target_labels:
        for lbl in target_labels:
            for v in str(lbl).split(","):
                v = v.strip()
                if v:
                    p.append(("target.label", v))
    if limit:
        p.append(("limit", limit))
    else:
        p.append(("limit", "all"))
    if extra_params:
        for k, v in extra_params.items():
            p.append((k, str(v)))
    p.append(("format","json"))
    return p
