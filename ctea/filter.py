#!/usr/bin/env python

import sys
import os
import argparse
from pathlib import Path
import re
from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd

# Set options
threads = 2
genome_build = "hg38"
rtea_script = "rtea_functions.R"
refdir = None

# Initialize RTEA environment
RTEA = {}
exec(Path(rtea_script).read_text(), {}, RTEA)

# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("-c", "--cteafile", help="ctea output file path")
parser.add_argument("-o", "--outfile", help="output directory")
parser.add_argument("-t", "--threads", type=int, default=threads, help="number of threads")
parser.add_argument("--genome_build", default=genome_build, help="reference genome build")
parser.add_argument("--rtea_script", default=rtea_script, help="rtea Rscript file")
parser.add_argument("--refdir", help="directory containing reference files")
args = parser.parse_args()

# Update options based on command line arguments
threads = args.threads
genome_build = args.genome_build
rtea_script = args.rtea_script
refdir = args.refdir or os.path.join(os.path.dirname(args.rtea_script), "ref", args.genome_build)

# Set options
options = {
    "mc.cores": threads,
    "genome_build": genome_build,
    "refdir": refdir
}
RTEA.update(options)

# Read ctea file
ctea = pd.read_table(args.cteafile)

# Perform filtering steps
ctea = RTEA['filterUnlocalized.ctea'](ctea)
ctea = RTEA['filterSimpleRepeat.ctea'](ctea)
ctea = RTEA['repeatPositon.ctea'](ctea)
ctea = RTEA['filterSimpleSite'](ctea)
ctea = RTEA['filterNoClip.ctea'](ctea, threads=threads)
ctea = RTEA['TEcoordinate'](ctea, threads=threads)

# Update class column
ctea.loc[ctea['isPolyA'], 'class'] = 'PolyA'

# Write the result to the outfile
ctea.to_csv(args.outfile, sep='\t', quotechar='', quoting=pd.QUOTE_NONE, na_rep='NA', index=False)
