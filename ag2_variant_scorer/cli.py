"""Command-line interface for AlphaGenome 2 variant scoring."""

import argparse
from pathlib import Path

from . import score, summarize


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="ag2-variant-scorer", description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    score_parser = commands.add_parser("score", help="Score variants or preview planned requests")
    score.add_score_arguments(score_parser)
    score.add_execution_arguments(score_parser)
    score_parser.set_defaults(resume=None)

    resume_parser = commands.add_parser("resume", help="Continue a saved scoring run")
    resume_parser.add_argument(
        "--run-dir", dest="resume", type=Path, required=True, help="Saved scoring run directory",
    )
    score.add_execution_arguments(resume_parser)

    summary_parser = commands.add_parser("summarize", help="Summarize gene effects from saved scores")
    summarize.add_arguments(summary_parser)

    args = parser.parse_args(argv)
    if args.command == "summarize":
        return summarize.run_cli(args)
    return score.run_cli(args, score_parser if args.command == "score" else resume_parser)
