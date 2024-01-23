from __future__ import annotations

import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Generator

import httpx
from tqdm import tqdm


def make_request(url: str) -> dict:
    with httpx.Client() as client:
        response = client.get(url)
        # Additional delay to simulate a slow request.
        time.sleep(1)

        return response.json()


def make_requests(
    urls: list[str],
) -> Generator[list[dict], None, None]:
    with tqdm(total=len(urls)) as pbar:
        with ThreadPoolExecutor(max_workers=5) as executor:
            futures = [executor.submit(make_request, url) for url in urls]
            for future in as_completed(futures):
                pbar.update(1)
                yield future.result()


def main() -> None:
    urls = [
        "https://httpbin.org/get",
        "https://httpbin.org/get?foo=bar",
        "https://httpbin.org/get?foo=baz",
        "https://httpbin.org/get?foo=qux",
        "https://httpbin.org/get?foo=quux",
    ]

    results = []
    for result in make_requests(urls):
        results.append(result)

    print(results)


if __name__ == "__main__":
    main()