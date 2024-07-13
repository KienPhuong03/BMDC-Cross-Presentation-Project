### IMPORTS
import os
import sys
import json
import re
import math
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import seaborn as sns

###Tasks: The filtered PSM are essentially PSMs that passed the quality filtering. 
# From this pool, we would like to know the fraction of peptides (note that you need to concatenate PSMs) that has a "Y" or "F" or "N" in the sequence,
# and from this subset, what is the breakdown of #peptides where all the presented "Y" or "F" or "N" are heavy labeled in all peptides,
# #peptides that are partially labeled (some PSMs of this peptide has heavy label but some don't), and #peptides that are not heavy labeled at all.

def run_labelling_analysis(filename: str):
    """
    Input: .csv file containing all the PSMs that passed the quality filtering.
    Output: a tuple containing
    - Dictionary storing number of fully labeled, partially labeled and not labeled unique peptides
    - Total number of unique peptides.
    Notation: 
    1. fully labeled means that all Y, F, N residues in a peptide is labeled in all spectra measurements.
    2. not labeled means either the peptide has no modification or irrelevant modifications.
    3. partially labeled is the rest.
    """
    peptide_data = pd.read_csv(filename)
    unique_peptide_data = peptide_data.groupby("Sequence")["Modifications"].apply(tuple).reset_index()
    num_peptides = unique_peptide_data.shape[0]
    labeled_aa = {"Y", "F", "N"}
    output = {
        "Fully Labeled": 0,
        "Partially Labeled": 0,
        "Not Labeled": 0,
        }
    
    #Modifications is now a tuple of all concatenated modifications on all measurements of a unique peptide.
    for index, row in unique_peptide_data.iterrows():
        peptide_seq = row["Sequence"]
        modifications = row["Modifications"]
        num_spectra = len(modifications)
        num_yfn_residues = sum([aa in labeled_aa for aa in peptide_seq])
        count = 0
        for mod in modifications:
            if mod is np.nan:
                continue
            # if the number of labeled residue in a spectra = number of Y,F,N residues in that peptide then such spectra is "fully labeled"
            num_labels_residue = sum(["Label" in token for token in mod.split(sep = ";")])
            if num_labels_residue == num_yfn_residues:
                count += 1
        # A peptide is fully labeled if all spectras are fully labeled. If none are labeled, then it's not labeled.
        # Otherwise the peptide is partially labeled.
        if count == 0:
            output["Not Labeled"] += 1
        elif count == num_spectra:
            output["Fully Labeled"] += 1
        else:
            output["Partially Labeled"] += 1
    return (output, num_peptides)

print(f"Result of CT2A sample is {run_labelling_analysis("SILAC_MHC_CT2A_filtered.csv")}")
print(f"Result of GL261 sample is {run_labelling_analysis("SILAC_MHC_GL261_filtered.csv")}")