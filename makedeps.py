#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
Script for generating Makefile dependencies.
"""

__all__ = [
    ]

# ---------------------------------------------------------

from path import *

# ---------------------------------------------------------

def print_dependencies(n):
    r"""
    Print the Makefile dependencies for size n.
    """
    print "PATHS-{n} :=".format(n=n)
    for path in iter_path(n):
        path_string = ''.join(map(str, path))
        print "PATHS-{n} += {path_string}".format(n=n, path_string=path_string)
    print "OUTFILES += $(PATHS-{n}:%=output/csf-%.py)".format(n=n)
    print "OUTFILES += $(PATHS-{n}:%=output/hess-%.py)".format(n=n)
    print "$(PATHS-{n}:%=output/csf-%.py): var/csf-size-{n}".format(n=n)
    print "\ttouch $@"

# ---------------------------------------------------------

def doctest():
    import doctest
    doctest.testmod(verbose=False)

# ---------------------------------------------------------

def argparse():
    import argparse
    parser = argparse.ArgumentParser(
        description='Generate Makefile dependencies for Dyck paths of size n.',
        )
    parser.add_argument(
        'n',
        type=int,
        choices=range(1, 11),
        metavar='n',
        help='The size of Dyck path to consider (from 1 to 10).',
        )
    args = parser.parse_args()
    return args.n

# ---------------------------------------------------------
if __name__ == '__main__':
    doctest()
    n = argparse()
    print_dependencies(n)
# ---------------------------------------------------------

