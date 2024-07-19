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

from refactored_SILAC_analysis import run_labeling_analysis


#Helper function to create output file:
def write_data_to_file(output_data, sample_names, file_path):
    with open(file_path, 'w') as file:
        for idx, sample in enumerate(sample_names):
            file.write(f"For {sample}, we have:\n")
            for label, (total, details) in output_data[idx]:
                file.write(f"There are total of {total} peptides with at least 1 labeled {label}, and the result is \n")
                file.write(f" {details}\n")
            file.write("______\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str, required=True)
    parser.add_argument('--data', type=list, required=True)
    parser.add_argument("--sample_names", type = list, required = True)
    parser.add_argument("--label", type = list, required = True)
    parser.add_argument("--sep", type = bool, required = True)
    parser.add_argument('--o', type=str,required=True)
    
    args = parser.parse_args()
    #set working directory to file path provided
    os.chdir = args.path
    for filename in args.data:
        output = run_labeling_analysis(filename, args.label, args.sep)
        write_data_to_file(output, args.sample_names, args.o)

if __name__ == "__main__":
    main()