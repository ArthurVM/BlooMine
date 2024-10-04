#! /usr/bin/env python3
"""
Command line interface for BlooMine.
"""
from .parser import parser_bloomine


def main():
    args = parser_bloomine.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser_bloomine.print_help()


if __name__ == "__main__":
    main()
