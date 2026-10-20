"""Run: python -m supermagx case.json result.json"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .case import run_case


def main() -> None:
    parser = argparse.ArgumentParser(description="SuperMagX phase-one calculation")
    parser.add_argument("case", type=Path, help="input case JSON")
    parser.add_argument("output", type=Path, help="output result JSON")
    args = parser.parse_args()
    case = json.loads(args.case.read_text(encoding="utf-8"))
    result = run_case(case)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2),
                           encoding="utf-8")
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
