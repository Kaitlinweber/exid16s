import argparse
import pathlib
import os 
import subprocess
from extract_16S_Barrnap import * 


def get_sample_name(file_path):
    '''Extracts sample name (stem of file name) from the file path.
    '''
    file_path = pathlib.Path(file_path)
    sample_ID = str(file_path.stem)
    return sample_ID

def extract_16S_Barrnap(list_of_fasta_files, barrnap_out):
    ''' Runs barrnap and extracts the 16S rDNa sequence from Barrnap results
    ''' 
    for fasta_file in list_of_fasta_files:
        sample_name = get_sample_name(fasta_file)
        barrnap_dir = "barrnap_result/"
        for parent_dir in barrnap_out:
            path = os.path.join(parent_dir, barrnap_dir)
            try:
                os.mkdir(path)
            except FileExistsError:
                pass
        subprocess.call(["barrnap", fasta_file, "--outseq", path + f"{sample_name}.fasta"])
        
        
        
        sequence_16S_dir = "FASTA_16S_sequence/"





if __name__ == '__main__':
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument('-i', '--input', type=pathlib.Path, 
                        default=[], nargs='+', help='Filepath to assembled input files')
    argument_parser.add_argument('-o', '--output', type=pathlib.Path, 
                        default=[], nargs='+', help='Filepath to output directory')
    argument_parser.add_argument('-db', '--database', type=pathlib.Path, 
                        default=[], nargs='+', help='Path to folder with databases which need to be combined')
    argument_parser.add_argument('-e', '--email', type=str, help='Enter email from NCBI account')
    args = argument_parser.parse_args()
    extract_16S = extract_16S_Barrnap(list_of_fasta_files=args.input, barrnap_out=args.output)



    