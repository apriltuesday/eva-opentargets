import logging
import requests
logger = logging.getLogger(__package__)


def json_request(url: str, payload: dict) -> dict:
    """Makes a GET request with the specified URL and payload, attempts to parse the result as a JSON string and
    return it as a dictionary, on failure raises an exception."""
    result = requests.get(url, data=payload)
    assert result.ok
    return result.json()


def retry_helper(function, kwargs: dict, retry_count: int):
    """
    Given a function, make a `retry_count` number of attempts to call it until it returns a value without raising
    an exception, and subsequently return this value. If all attempts to run the function are unsuccessful, return None.

    :param function: Function that could need multiple attempts to return a value
    :param kwargs: Dictionary with function's keyword arguments
    :param retry_count: Number of attempts to make
    :return: Returned value of the function.
    """
    for retry_num in range(retry_count):
        try:
            return function(**kwargs)
        except Exception as e:
            logger.warning("Attempt {}: failed running function {} with kwargs {}".format(retry_num, function, kwargs))
            logger.warning(e)
    logger.warning("Error on last attempt, skipping")
    return None


def request_retry_helper(url: str, payload: dict = None, retry_count: int = 4):
    """Makes a GET request with the specified URL and payload via `json_request`, makes several attempts and handles
    the exceptions via `retry_helper`."""
    return retry_helper(json_request, {'url': url, 'payload': payload}, retry_count)
