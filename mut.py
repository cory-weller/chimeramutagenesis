#!/usr/bin/env python


# PARSE ARGUMENTS
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta_1_filename', type=str,
    help="File name for first (i.e. upstream) sequence.")
parser.add_argument('fasta_2_filename', type=str,
    help="File name for second (i.e. downstream) sequence.")
parser.add_argument('codon_table', type=str,
    help='''File name for tab-delimited 4-column codon table.
        Should include no header, with columns representing
        triplet, amino-acid code, relative frequency,
        fraction per 1000, and count. See README for help
        finding and generating this file.''')
parser.add_argument("-f", "--filter-freq",
    type=float,
    help='''Minimum required fraction for including a codon.
    Triplets with (amino-acid specific) frequency
    below this provided value will be excluded.
    Default: 0.10 (i.e. include all triplets with 10%%
    or greater represention for its amino acid).''')

args = parser.parse_args()

print(args.filter_freq)


# IMPORT MODULES
import pandas as pd
import numpy as np
import sys
import functions as f
from functions import fasta

# RUN CODE

## 
codons = pd.read_csv(   "codonTableCleaned.txt", 
             sep = "\t",
             header = None,
            )

codons.columns = ["triplet", "AA", "freq", "fkp", "count"]


putativeCodonsTable = codons[codons.freq >= args.filter_freq]
putativeCodonsTable['triplet'] = putativeCodonsTable['triplet'].str.replace("U","T")


assert putativeCodonsTable['AA'].nunique() == 21, (
    "Codon frequency filter of %s eliminates some amino acids. Try a lower filter." % (args.filter_freq)
)

allAminoAcids = list(putativeCodonsTable['AA'].unique())

# retrieve rows with maximum frequency, grouped by amino acid
abundantCodons = putativeCodonsTable.loc[putativeCodonsTable.groupby('AA')['count'].idxmax()]


f.sampleCodon(putativeCodonsTable, "F")



exit