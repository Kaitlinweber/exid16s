import argparse
import pathlib
import os 
import glob
import subprocess
from .extract16S_Barrnap import extraction_16Ssequence, create_output_file
from .make_summary_kreport import read_kreport


def get_sample_name(file_path):
    '''Extracts sample name (stem of file name) from the file path.
    '''
    file_path = pathlib.Path(file_path)
    sample_ID = str(file_path.stem)
    return sample_ID


def run_barrnap(list_of_fasta_files, barrnap_out):
    ''' Runs barrnap with fasta files to extract rDNA
    '''  
    for parent_dir in barrnap_out:
            try: 
                os.mkdir(parent_dir)
            except FileExistsError:
                pass

    for fasta_file in list_of_fasta_files:
        sample_name = get_sample_name(fasta_file)
        barrnap_dir = "barrnap_result/"
        path = os.path.join(parent_dir, barrnap_dir)
        try:
            os.mkdir(path)
        except FileExistsError:
            pass
        subprocess.call(["barrnap", fasta_file, "--outseq", path + f"{sample_name}.fasta"], check=True)
    return parent_dir

def extract_16S_barrnap(path):
    ''' Extracts the 16S sequence from the barrnap results 
    ''' 
    list_barrnap = []
    barrnap_result = str(path) + '/barrnap_result'
    for root, dirs, all_files in os.walk(barrnap_result):
        for data_files in all_files:
            files = os.path.join(root, data_files)
            list_barrnap.append(files)
    list_all_samples = extraction_16Ssequence(list_barrnap)

    sequence_16S_dir = "FASTA_16S_sequence"
    output_path = os.path.join(path, sequence_16S_dir)

    create_output_file(list_all_samples, list_barrnap, output_path)
    return output_path


def run_Kraken2(parent_dir, path_16s, database_file_path):
    ''' Uses 16S sequence from barrnap to run Kraken2
    '''
    kraken_folder = "Kraken2_kreports/"
    kraken_output_folder = os.path.join(parent_dir, kraken_folder)
    try:
        os.mkdir(kraken_output_folder)
    except FileExistsError:
        pass
    
    for database_files in database_file_path:
        database = database_files
    
    for root, dirs, all_files in os.walk(path_16s):
        for data_files in all_files:
            fasta_file_path = os.path.join(root, data_files)
            sample_name = get_sample_name(fasta_file_path)
            subprocess.call(["kraken2", "-db", database, fasta_file_path, "--report", kraken_output_folder + f"{sample_name}.kraken2_kreport"], check=True)
    
    kreport_file_path = os.listdir(kraken_output_folder)
    kreports = glob.glob(kraken_output_folder + '*.kraken2_kreport')
    summary_file = read_kreport(kreports)
    summary_file.to_csv(kraken_output_folder + 'summary_file_kreport.csv', index=False)

