from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
import argparse
import pathlib
import os
import itertools


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
    sequence_list_organism = []
    sequence_list_definition = []
    test_list = []
    #always_print = False
    #organism_dict = {}
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
                        print(scientific_name_organism)
                        scientific_name_list_organism.append(scientific_name_organism)
                        #print("scientific name list", scientific_name_list_organism)
            
            
            if "ORIGIN" in new_contents:
        
                always_print = True
            if always_print:

                origin = line.replace('ORIGIN', '').replace('//', '')
                sequence = ''.join((x for x in origin if not x.isdigit()))
                print(sequence)
                    # if 'ORIGIN' in line:
                    #     always_print = True
                    # if always_print:
                    
                    #     test_2_list = []
                    #     origin = line.replace('ORIGIN', '').replace('//', '')
                    #     sequence = ''.join((x for x in origin if not x.isdigit()))
                    #     sequence = ''.join(sequence.split())
                    #     # test_list.append(sequence)
                    #     # sequence = ''.join(test_list)
                        
                    #     test_list.append(sequence)
                    #     #sequence_list_organism = itertools.chain(*test_list)
                    #     for e in test_list:

                    #         sequence_list_organism.append(test_list)
                    #         print(sequence_list_organism)
                        
                        

                    
                        #sequence = sequence + sequence
                        #sequence = ''.join(sequence)
                        #print("this is sequence", sequence)
                        # sequence_list_organism.append(sequence)
                        # print(sequence_list_organism)
                        #SPLIT TUSSEN SEQUENCE WEGHALEN EN SEQUENCE LIST MOET UIT DE LOOP 


            else:
                for line in contents:
                    line = line.decode("iso-8859-1")
                    if 'DEFINITION' in line:
                        if len(line) > 14:
                            elements = line.split(' ')
                            #new_elements = list(filter(None, elements))
                            genus = elements[2]
                            species = elements[3]
                            scientific_name_definition = genus + ' ' + species
                            scientific_name_list_definition.append(scientific_name_definition)
                            #print(scientific_name_list_definition)
                            # taxonomy_ID_definition = get_taxonomy_ID(scientific_name_definition)
                            # print(taxonomy_ID_definition)

                    if 'ORIGIN' in line:
                        always_print = True
                    if always_print:
                        origin = line.replace('ORIGIN', '').replace('//', '')
                        sequence = ''.join((x for x in origin if not x.isdigit()))
                        sequence = ''.join(sequence.split())
                        sequence = ''.join(sequence.strip('/n'))
                        sequence_list_definition.append(sequence)
                        #print(sequence)
                        #print(sequence_list_definition)
            

            #add taxonomic ID 
            # tax_id_string = '|kraken:taxid|'
            # header = '>sequence' + str(count+1) + tax_id_string + taxonomy_ID_organism + "\t" + scientific_name_organism
            # print(header)
                    
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