from Bio import SeqIO
import argparse
import pathlib
import os
 

def get_sample_name(file_path):
    '''Extracts sample name (stem of file name) from the file path.
    '''
    file_path = pathlib.Path(file_path)
    sample_ID = str(file_path.stem)
    return sample_ID

def get_multiple_sample_names(list_fasta_files):
    '''Makes list of names from Barrnap's output FASTA files
    '''
    sample_names = []
    for file_path in list_fasta_files:
        sample_ID = get_sample_name(file_path)
        sample_names.append(sample_ID)
    return sample_names

def extraction_16Ssequence(list_of_fasta_files):
    ''' Reads barrnaps FASTA file with rRNA sequences and
    extracts only the 16S rRNA sequence. 
    ''' 

    record_list = []
    for fasta_file in list_of_fasta_files:
        for record in SeqIO.parse(fasta_file, "fasta"):
            extract_record = '16S_rRNA'
            if record.id.startswith(extract_record):
                record_list.append(record)
    return record_list

def create_output_file(record_list, list_of_fasta_files, output_path):
    ''' Creates FASTA file with the extracted 16S rRNA sequence 
    ''' 
 
    try:
        os.mkdir(output_path)
    except FileExistsError:
        pass
    
    sample_number = get_multiple_sample_names(list_of_fasta_files)
    
    for file_name, record in zip(sample_number, record_list):
        with open(os.path.join(output_path, f"{file_name}.fasta"), "w") as output_handle:
            SeqIO.write(record, output_handle, 'fasta')
            



if __name__ == '__main__':
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument('-i', '--input', type=pathlib.Path, 
                        default=[], nargs='+', help='Path to Barrnap result files')
    argument_parser.add_argument('-o', '--output', type=pathlib.Path, 
                        default=[], nargs='+', help='Filepath to output directory')
    args = argument_parser.parse_args()
    record_list = extraction_16Ssequence(list_of_fasta_files=args.input)
    create_output_file(record_list, list_of_fasta_files=args.input, output_path=args.output)
    