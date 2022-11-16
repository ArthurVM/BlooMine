""" General utility funtions for Afanc
"""

import shutil
import re
import time
import sys
from os import path, listdir
from time import localtime, strftime


def isFile(filename):
    """ Checks if a path is an existing file """

    if not path.isfile(filename):
        vprint("MAIN", f"No file found at {filename}", "prRed", )
        exit(3)
    else:
        return path.abspath(path.realpath(path.expanduser(filename)))


def isDir(dirname):
    """ Checks if a path is an existing directory """

    if not path.isdir(dirname):
        vprint("MAIN", f"No file found at {dirname}", "prRed")
        exit(3)
    else:
        return path.abspath(path.realpath(path.expanduser(dirname)))


###################
# PRINT UTILITIES #
###################

colfuncs = {}
colfunc = lambda f: colfuncs.setdefault(f.__name__, f)

@colfunc
def prRed(sp):
    return f"\033[91m{sp}\033[00m"

@colfunc
def prGreen(sp):
    return f"\033[92m{sp}\033[00m"

@colfunc
def prYellow(sp):
    return f"\033[93m{sp}\033[00m"


def vprint(subprocess, info_text, colour, f=sys.stdout, end="\n"):
    """ controls process output
    """

    time = strftime("%H:%M:%S", localtime())

    print(
        f"\n{time} {colfuncs[colour](subprocess)} :: {info_text}",
        end=end,
        file=f,
        flush=True,
    )
