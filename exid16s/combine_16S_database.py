from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
import argparse
import pathlib
import os


def get_database_files(file_path):
    ''' Iterates over each directory, collecting all the database files
    '''
    list_of_folders = []
    list_of_files = []
    for path in file_path:
        directory_path = path.as_posix()
        list_of_folders.append(directory_path)
        for folder in list_of_folders:
            for root, dirs, all_files in os.walk(folder):
                for data_files in all_files:
                    files = os.path.join(root, data_files)
                    list_of_files.append(files)
    return list_of_files

def extract_from_database_files(files):
    ''' Extracts each scientific name, taxonomy ID and sequence 
    ''' 
    scientific_name_list_organism = []
    scientific_name_list_definition = []
    sequence_list = []
    #always_print = False
    organism_dict = {}
    for database_file in files:
        with open(database_file, 'rb') as data:
            always_print = False
            contents = data.readlines()
            
            merged_contents_list = [] 
            for line in contents:
                line = line.decode("iso-8859-1")
                merged_contents_list.append(line)
            
            #dit is nu een string
            new_contents = ' '.join(merged_contents_list)

            if 'ORGANISM' in new_contents:
                for line in contents:
                    line = line.decode("iso-8859-1")
                    if 'ORGANISM' in line:
                        scientific_name_organism = line.replace('ORGANISM', '').replace('    ', '')
                        #print(scientific_name_organism)
                        taxonomy_ID_organism = get_taxonomy_ID(scientific_name_organism)
                        #print (taxonomy_ID_organism)
                    if 'ORIGIN' in line:
                        always_print = True
                    if always_print:
                        origin = line.replace('ORIGIN', '').replace('//', '')
                        sequence = ''.join((x for x in origin if not x.isdigit()))
                        sequence = ''.join(sequence.split())

            else:
                for line in contents:
                    line = line.decode("iso-8859-1")
                    if 'DEFINITION' in line:
                        if len(line) > 14:
                            elements = line.split(' ')
                            #new_elements = list(filter(None, elements))
                            genus = elements[2]
                            species = elements[3] # TODO: species/subspecies bespreken 
                            scientific_name_definition = genus + ' ' + species
                            print(scientific_name_definition)
                            taxonomy_ID_definition = get_taxonomy_ID(scientific_name_definition)
                            print(taxonomy_ID_definition)

                    if 'ORIGIN' in line:
                        always_print = True
                    if always_print:
                        origin = line.replace('ORIGIN', '').replace('//', '')
                        sequence = ''.join((x for x in origin if not x.isdigit()))
                        sequence = ''.join(sequence.split())
                        sequence = ''.join(sequence.strip('/n'))
                        print(sequence)

                    
def get_taxonomy_ID(scientific_name):
    ''' Search for the taxonomy ID based on the scientific name 
    ''' 
    Entrez.email = args.email 
    handle = Entrez.esearch(db="Taxonomy", term=scientific_name) 
    record = Entrez.read(handle)
    taxonomy_ID = record["IdList"][0]
    return taxonomy_ID


if __name__ == '__main__':
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument('-db', '--database', type=pathlib.Path, 
                        default=[], nargs='+', help='Path to folder with databases which need to be combined')
    argument_parser.add_argument('-e', '--email', type=str, help='Enter email from NCBI account')
    args = argument_parser.parse_args()
    extract_information = get_database_files(file_path=args.database)
    extract_from_database_files(files=extract_information)