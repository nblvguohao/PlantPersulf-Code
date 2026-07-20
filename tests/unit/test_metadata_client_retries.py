from collections.abc import Callable
from typing import Any

import pytest

from plantpersulf.download import geo, iprox, pride


@pytest.mark.parametrize(
    ("module", "client", "request_method"),
    [
        (pride, pride.PrideClient(retries=2), "_request_json"),
        (geo, geo.GeoClient(retries=2), "_request_text"),
        (iprox, iprox.IproxClient(retries=2), "_request_xml"),
    ],
)
def test_metadata_clients_retry_socket_timeouts(
    monkeypatch: pytest.MonkeyPatch,
    module: Any,
    client: Any,
    request_method: str,
) -> None:
    attempts = 0

    def raise_timeout(*args: object, **kwargs: object) -> None:
        nonlocal attempts
        attempts += 1
        raise TimeoutError("policy-test timeout")

    monkeypatch.setattr(module, "urlopen", raise_timeout)
    monkeypatch.setattr(module.time, "sleep", lambda _: None)
    request: Callable[[str], object] = getattr(client, request_method)

    with pytest.raises(RuntimeError, match="request failed"):
        request("https://example.invalid/metadata")
    assert attempts == 2
