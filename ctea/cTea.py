#!/usr/bin/env python

import argparse
import pandas as pd
from data.table import fread, fwrite
from magrittr import setDTthreads

parser = argparse.ArgumentParser(description="Process ctea output file.")
parser.add_argument("-c", "--cteafile", help="ctea output file path", required=True)
parser.add_argument("-o", "--outfile", help="output directory", required=True)
parser.add_argument("-t", "--threads", type=int, default=2, help="number of threads")
parser.add_argument("--genome_build", default="hg38", help="reference genome build")
parser.add_argument("--rtea_script", default="rtea_functions.R", help="rtea Rscript file")
parser.add_argument("--refdir", help="directory containing reference files")

args = parser.parse_args()

if args.refdir is None:
    args.refdir = "ref/{}".format(args.genome_build)

# Set threads
setDTthreads(args.threads)

# Load RTEA
RTEA = {}
exec(open(args.rtea_script).read(), {}, RTEA)

# Load ctea data
ctea = fread(args.cteafile)

# Data processing steps
ctea = RTEA["readctea"](ctea)
ctea = RTEA["filterUnlocalized.ctea"](ctea)
ctea = RTEA["filterSimpleRepeat.ctea"](ctea)
ctea = RTEA["repeatPositon.ctea"](ctea)
ctea = RTEA["filterSimpleSite"](ctea)
ctea = RTEA["filterNoClip.ctea"](ctea, threads=args.threads)
ctea = RTEA["TEcoordinate"](ctea, threads=args.threads)

# Modify DataFrame in Python style
ctea.loc[ctea['isPolyA'] == True, 'class'] = 'PolyA'
ctea = ctea[ctea['isPolyA'] | (ctea['TEscore'] > 0)]

# Save the result
ctea.to_csv(args.outfile, sep="\t", index=False, quoting=None, na_rep="NA")
