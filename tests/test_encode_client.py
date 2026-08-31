"""Tests for encodefetch.encode_client.encode_get's HTTP-status handling.

ENCODE's /search/ endpoint returns HTTP 404 even for a perfectly valid
query that just matches zero experiments (confirmed against the live API:
identical 404 + "No results found" body whether the filter is valid-but-
unmatched or a nonsense field name - ENCODE doesn't distinguish). That is
a legitimate empty result, not an error, and must not raise. A 404 for an
item lookup (e.g. /experiments/{accession}/ for an accession that doesn't
exist) is a genuine "not found" and should raise EncodeNotFoundError so
callers can catch it specifically (as opposed to requests.HTTPError, which
still covers real failures like 5xx server errors).
"""

import json

import pytest
import requests

from encodefetch.encode_client import EncodeNotFoundError, encode_get


def _fake_response(status_code, json_body=None, text_body=None, url="https://www.encodeproject.org/fake/"):
    resp = requests.Response()
    resp.status_code = status_code
    resp.url = url
    if json_body is not None:
        resp._content = json.dumps(json_body).encode("utf-8")
    elif text_body is not None:
        resp._content = text_body.encode("utf-8")
    else:
        resp._content = b""
    return resp


def test_encode_get_returns_json_on_success(monkeypatch):
    payload = {"@graph": [{"accession": "ENCSR000AAA"}], "total": 1}
    monkeypatch.setattr(requests.Session, "send", lambda self, *a, **k: _fake_response(200, payload))
    assert encode_get("/search/") == payload


def test_encode_get_returns_empty_search_result_on_404_with_graph(monkeypatch):
    """A search that legitimately matched zero experiments - must not raise."""
    payload = {"@graph": [], "total": 0, "notification": "No results found"}
    monkeypatch.setattr(requests.Session, "send", lambda self, *a, **k: _fake_response(404, payload))
    result = encode_get("/search/")
    assert result == payload
    assert result["@graph"] == []


def test_encode_get_raises_not_found_on_404_without_graph(monkeypatch):
    """An item lookup (e.g. a bad accession) - a genuine "not found"."""
    payload = {
        "@type": ["HTTPNotFound", "Error"],
        "status": "error",
        "code": 404,
        "title": "Not Found",
        "description": "The resource could not be found.",
    }
    monkeypatch.setattr(requests.Session, "send", lambda self, *a, **k: _fake_response(404, payload))
    with pytest.raises(EncodeNotFoundError) as exc_info:
        encode_get("/experiments/ENCSRBADACCESSION999/")
    assert "ENCSRBADACCESSION999" in str(exc_info.value)


def test_encode_get_raises_not_found_on_404_with_non_json_body(monkeypatch):
    """Even a non-JSON 404 body (e.g. an HTML error page) must degrade to
    EncodeNotFoundError rather than crashing on the .json() call itself."""
    monkeypatch.setattr(
        requests.Session, "send", lambda self, *a, **k: _fake_response(404, text_body="<html>not found</html>")
    )
    with pytest.raises(EncodeNotFoundError):
        encode_get("/experiments/ENCSRBADACCESSION999/")


def test_encode_get_raises_http_error_on_server_error(monkeypatch):
    """A real server-side failure (not a "not found") still raises the
    standard requests.HTTPError - callers should not silently skip these
    the way they skip a bad accession."""
    monkeypatch.setattr(requests.Session, "send", lambda self, *a, **k: _fake_response(503))
    with pytest.raises(requests.exceptions.HTTPError):
        encode_get("/search/")


def test_encode_not_found_error_message_includes_path():
    err = EncodeNotFoundError("/experiments/ENCSRBADACCESSION999/")
    assert "/experiments/ENCSRBADACCESSION999/" in str(err)
