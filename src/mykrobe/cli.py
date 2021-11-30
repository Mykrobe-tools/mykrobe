#! /usr/bin/env python
from __future__ import print_function
from mykrobe.parser import parser


def main():
    args = parser.parse_args()
    if hasattr(args, "func"):
        args.func(parser, args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
