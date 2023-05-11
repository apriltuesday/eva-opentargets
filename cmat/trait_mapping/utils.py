import logging
import requests
from retry import retry
logger = logging.getLogger(__package__)


@retry(exceptions=(ConnectionError, requests.RequestException), logger=logger,
       tries=8, delay=2, backoff=1.2, jitter=(1, 3))
def json_request(url: str, payload: dict = None, method=requests.get) -> dict:
    """Makes a request of a specified type (by default GET) with the specified URL and payload, attempts to parse the
    result as a JSON string and return it as a dictionary, on failure raises an exception."""
    result = method(url, data=payload)
    result.raise_for_status()
    return result.json()
