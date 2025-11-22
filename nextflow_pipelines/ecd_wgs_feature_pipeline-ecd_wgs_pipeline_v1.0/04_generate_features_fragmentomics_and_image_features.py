import argparse
import pandas as pd 
import os
from tqdm import tqdm
from feature_class import *

##### Feature class definition
def main():
    parser = argparse.ArgumentParser(description='Generate an image matrix.')
    parser.add_argument('--input', type=str, required=True, help='Path to the pre-processed frag.tsv files from 01 and 02')
    parser.add_argument('--output', type=str, required=False, help='Path to save output feature csv files')
    
    args = parser.parse_args()
    input_tsv = args.input
    outputdir = args.output
    os.system(f"mkdir -p {outputdir}")
    
    output_obj = WGS_GW_Image_features(
        input_tsv = input_tsv,
        outputdir = outputdir)
    
    output_obj.generate_flen_feature(save_feature = True)
    output_obj.generate_em_feature(save_feature = True)
    output_obj.generate_nuc_feature(save_feature = True)
    output_obj.generate_all_image_features()
    
if __name__ == '__main__':
    main()