#!/usr/bin/env python
"""Send a Feishu bot message — run after the GPU chain completes or fails.

Usage (called automatically by the monitor wrapper below):
    python scripts/notify_feishu.py "message text"
"""

import json
import sys
import urllib.request

WEBHOOK = (
    "https://open.feishu.cn/open-apis/bot/v2/hook/b9c6b3d4-f289-4308-bb7f-a72e7ddbf9fe"
)


def send(text: str) -> None:
    payload = json.dumps({"msg_type": "text", "content": {"text": text}}).encode()
    req = urllib.request.Request(
        WEBHOOK, data=payload, headers={"Content-Type": "application/json"}
    )
    urllib.request.urlopen(req, timeout=10)


if __name__ == "__main__":
    send(" ".join(sys.argv[1:]))
