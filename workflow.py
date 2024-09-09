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
import seaborn
import logging
logging.basicConfig(level=logging.INFO)
from refactored_SILAC_analysis import run_labeling_analysis


#Helper function to create output file:
def write_data_to_file(output_data, sample_names, file_path):
    print(f"Length of output data: {len(output_data)}, and number of sample names: {len(sample_names)}")
    with open(file_path, 'w') as file:
        for idx, sample in enumerate(sample_names):
            file.write(f"For {sample}, we have:\n")
            for label, (total, details) in output_data[idx]:
                file.write(f"There are total of {total} peptides with at least 1 {label}, and the result is \n")
                file.write(f" {details}\n")
            file.write("______\n")

def main():
    parser = argparse.ArgumentParser(description = "Write peptide labeling data to a specified output file.")
    parser.add_argument('--path', type = str, required = True, help = 'The path to the data files')
    parser.add_argument('--data', nargs = '+', required = True, help = 'List of data file names')
    parser.add_argument('--sample_names', nargs='+', required=True, help = 'List of sample names')
    parser.add_argument('--label', nargs = '+', required = True, help = 'List of labels')
    parser.add_argument('--sep', action ='store_true', help = 'Separator flag')
    parser.add_argument('--no-sep', dest = 'sep', action = 'store_false', help = 'Separator flag')
    parser.set_defaults(sep = True)
    parser.add_argument('--o', type=str, required = True, help = 'The output file path where data will be written')
    
    args = parser.parse_args()
    #set working directory to file path provided
    os.chdir(args.path)
    all_output_data = []
    for filename in args.data:
        logging.info(f"Starting to process file: {filename}")
        try:
            output = run_labeling_analysis(filename, args.label, args.sep)
            all_output_data.append(output)
            logging.info(f"Successfully processed file: {filename}")
        except Exception as e:
            logging.error(f"Error processing file {filename}: {e}")
   
    write_data_to_file(all_output_data, args.sample_names, args.o)
    print("Done!")

if __name__ == "__main__":
    main()
    