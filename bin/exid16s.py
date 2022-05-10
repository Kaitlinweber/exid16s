import argparse
import pathlib
import os 
import subprocess
from extract16S_Barrnap import * 


def get_sample_name(file_path):
    '''Extracts sample name (stem of file name) from the file path.
    '''
    file_path = pathlib.Path(file_path)
    sample_ID = str(file_path.stem)
    return sample_ID

def run_barrnap(list_of_fasta_files, barrnap_out):
    ''' Runs barrnap with fasta files to extract rDNA
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
        #subprocess.call(["barrnap", fasta_file, "--outseq", path + f"{sample_name}.fasta"])
        return parent_dir

def extract_16S_barrnap(path):
    ''' Extracts the 16S sequence from the barrnap results 
    ''' 
    list_barrnap = []
    for root, dirs, all_files in os.walk(path):
        for data_files in all_files:
            files = os.path.join(root, data_files)
            list_barrnap.append(files)

    record_list = extraction_16Ssequence(list_barrnap)

    sequence_16S_dir = "FASTA_16S_sequence"
    output_path = os.path.join(path, sequence_16S_dir)

    create_output_file(record_list, list_barrnap, output_path)
    return output_path


def run_Kraken2(parent_dir, path_16s, database):
    ''' Uses 16S sequence from barrnap to run Kraken2
    '''
    # Get fasta files with 16S sequence 
    for root, dirs, all_files in os.walk(path_16s):
        for data_files in all_files:
            files = os.path.join(root, data_files)
    
    #database is posixpath 
  
    kraken_folder = "Kraken2_kreports"
    kraken_output_folder = os.path.join(parent_dir, kraken_folder)
    try:
        os.mkdir(kraken_output_folder)
    except FileExistsError:
        pass
    print(kraken_output_folder)



    

    




if __name__ == '__main__':
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument('-i', '--input', type=pathlib.Path, 
                        default=[], nargs='+', help='All input files, this would be the assembled WGS data')
    argument_parser.add_argument('-o', '--output', type=pathlib.Path, 
                        default=[], nargs='+', help='Filepath to output directory')
    argument_parser.add_argument('-db', '--database', type=pathlib.Path, 
                        default=[], nargs='+', help='Path to folder with database')
    argument_parser.add_argument('-e', '--email', type=str, help='Enter email from NCBI account')
    args = argument_parser.parse_args()
    barrnap_result= run_barrnap(list_of_fasta_files=args.input, barrnap_out=args.output)
    sequence_16s = extract_16S_barrnap(path=barrnap_result)
    run_Kraken2(parent_dir=barrnap_result, path_16s=sequence_16s, database=args.database)
    

