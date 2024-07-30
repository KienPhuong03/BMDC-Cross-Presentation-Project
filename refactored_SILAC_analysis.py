## IMPORTS
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

def peptide_label_helper(peptides_df : pd.DataFrame, labeled_residues : list):
    """
    Helper function
    Input: two arguments
    - residues: a list of amino acids that dictates what kind of labeled amino acid we are accounting for
    - peptides_df: dataframe of unique peptides
    Output: total number of peptides that has the residues specified in variable condition,
    and a dictionary specifying how many peptides are fully labeled, partially labeled, and
    not labeled.
    """
    #build regex engine for parsing modifications
    pattern = rf'[{"".join(labeled_residues)}]\d{{1,2}}\(Label.*\)$'
    regex = re.compile(pattern)

    #initiate result
    fully_labeled = 0
    partially_labeled = 0
    not_labeled = 0
    total_count_peptides = 0
    verifier = []


    for index, row in peptides_df.iterrows():
        
        peptide_seq = row["Sequence"]
        modifications = row["Modifications"]
        num_spectra = len(modifications)
        num_residues_in_question = sum([aa in labeled_residues for aa in peptide_seq])
        if num_residues_in_question != 0:
            total_count_peptides += 1
        else:
            continue
        count = 0
        
        for mod in modifications:
            if mod is np.nan:
                continue
            # if the number of labeled residue in a spectra = number of Y,F,N residues in that peptide then such spectra is "fully labeled"
            num_labeled_residue = sum([bool(regex.match(token.strip())) for token in mod.split(sep = ";")])

            verifier.append((num_residues_in_question, num_labeled_residue))

            if num_labeled_residue == num_residues_in_question:
                count += 1
        # A peptide is fully labeled if all spectras are fully labeled. If none are labeled, then it's not labeled.
        # Otherwise the peptide is partially labeled.
        if count == 0:
            not_labeled += 1
        elif count == num_spectra:
            fully_labeled += 1
        else:
            partially_labeled += 1
    
    return (total_count_peptides, 
        {
            "Fully labeled": fully_labeled,
            "Partially labeled": partially_labeled,
            "Not labeled": not_labeled,
        },
    )

def run_labeling_analysis(filename: str, labeled_aa_list, separate = False):
    """
    Notation: 
    1. fully labeled means that all Y, F, N residues in a peptide is labeled in all spectra measurements.
    2. not labeled means either the peptide has no modification or irrelevant modifications.
    3. partially labeled is the rest.
    Input: .csv file containing all the PSMs that passed the quality filtering.
    Output: return a list of tuples of the form (num of unique peptides, labeling breakdown dictionary)
    , in which dictionary will report number of fully labeled, partially labeled, and not labeled peptides.
    If separate parameter is true, each dicionary reflect data for each labeled amino acid separately.
    """
    #Read and Consolidate Peptide Data
    peptide_data = pd.read_csv(filename)
    peptide_data = peptide_data.groupby("Sequence")["Modifications"].apply(tuple).reset_index()
    num_peptides = peptide_data.shape[0]
    
    if separate: # Analyze each labeled amino acid separately
        outputs = []
        for aa in labeled_aa_list:
            outputs.append((aa, peptide_label_helper(peptide_data, [aa])))
        return outputs

    else: # Process all residues in aggregate (total peptides = peptides with at least one of the labeled residues)
        total_count, output = peptide_label_helper(peptide_data, labeled_aa_list)
        return [(total_count, output)]

