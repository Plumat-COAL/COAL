"""Entry point: ``python -m coal_app`` / ``coal-gui``."""
import sys


def main(argv=None) -> int:
    print("coal_app: the GUI is not built yet", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
