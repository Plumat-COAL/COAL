"""Hugging Face Spaces entry: the COAL GUI + MCP server, listening on 0.0.0.0:7860.

Used by space/app.py. Locally you can try the same thing with ``python -m coal_app.space_app``.
"""
from __future__ import annotations

import os

os.environ.setdefault("GRADIO_ANALYTICS_ENABLED", "False")


def main() -> None:
    from . import ui
    port = int(os.environ.get("PORT", os.environ.get("GRADIO_SERVER_PORT", "7860")))
    ui.launch(host="0.0.0.0", port=port, mcp=True, demo_mode=False)


if __name__ == "__main__":
    main()
