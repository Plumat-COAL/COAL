"""Entry point: ``python -m coal_app`` / ``coal-gui``.

    python -m coal_app                 web GUI on http://127.0.0.1:7860 (+ MCP at /gradio_api/mcp/sse)
    python -m coal_app --demo          same, opening with the sample report
    python -m coal_app --mcp           MCP server over stdio for AI agents (no Gradio needed)
    python -m coal_app --check         print which engine features were detected and exit
"""
from __future__ import annotations

import argparse
import os
import sys


def _is_loopback(host: str) -> bool:
    import ipaddress
    if host.strip().lower() == "localhost":
        return True
    try:
        return ipaddress.ip_address(host.strip("[]")).is_loopback
    except ValueError:
        return False


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="coal-gui", description="COAL web GUI and MCP server.")
    ap.add_argument("--host", default="127.0.0.1", help="address to listen on (default 127.0.0.1, this computer only)")
    ap.add_argument("--port", type=int, default=7860)
    ap.add_argument("--allow-remote", action="store_true",
                    help="allow a --host other than this computer (DANGER: the GUI has no login, and anyone who can "
                         "reach it can make COAL read any file you can read)")
    ap.add_argument("--demo", action="store_true", help="open with the sample report loaded")
    ap.add_argument("--no-mcp", action="store_true", help="do not serve the MCP tools from the web GUI")
    ap.add_argument("--mcp", action="store_true", help="run only the MCP server, over stdio")
    ap.add_argument("--check", action="store_true", help="show detected engine features and exit")
    args = ap.parse_args(argv)
    os.environ.setdefault("GRADIO_ANALYTICS_ENABLED", "False")
    if not _is_loopback(args.host) and not args.mcp and not args.check:
        if not args.allow_remote:
            print(f"coal-gui: refusing to listen on {args.host}: the GUI has no login, and anyone who can reach it "
                  "could read any file you can read through it. Use the default 127.0.0.1, or add --allow-remote "
                  "if you really mean it.", file=sys.stderr)
            return 2
        print(f"coal-gui: warning: listening on {args.host} with --allow-remote: anyone who can reach this address "
              "can read your files through COAL.", file=sys.stderr)

    if args.check:
        from . import core
        f = core.detect_features()
        print(f"coal {f.version}: formats={','.join(f.output_formats)} report={f.report} "
              f"json={f.json} keep_going={f.keep_going}")
        print(f"examples: {', '.join(core.available_examples()) or 'none found'}")
        return 0
    if args.mcp:
        try:
            from . import mcp_stdio
            return mcp_stdio.main()
        except ImportError:
            print("coal-gui: the MCP server needs the 'mcp' package (version 2): pip install -e \".[mcp]\"",  # CHANGED(IMP-M-4): quoted
                  file=sys.stderr)
            return 1
    try:
        from . import ui
    except ImportError as e:
        print(f"coal-gui: the web GUI needs Gradio ({e}). Install it with:  pip install -e \".[gui]\"",  # CHANGED(IMP-M-4): quoted
              file=sys.stderr)
        return 1
    ui.launch(host=args.host, port=args.port, mcp=not args.no_mcp, demo_mode=args.demo)
    return 0


if __name__ == "__main__":
    sys.exit(main())
