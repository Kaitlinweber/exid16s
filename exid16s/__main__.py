import argparse
import pathlib
from exid16s import exid16s

def main():
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument('-i', '--input', type=pathlib.Path, 
                        default=[], nargs='+', help='All input files, this would be the assembled WGS data')
    argument_parser.add_argument('-o', '--output', type=pathlib.Path, 
                        default=[], nargs='+', help='Filepath to output directory')
    argument_parser.add_argument('-db', '--database', type=pathlib.Path, 
                        default=[], nargs='+', help='Path to folder with database')
    argument_parser.add_argument('-e', '--email', type=str, help='Enter email from NCBI account')
    args = argument_parser.parse_args()
    barrnap_result= exid16s.run_barrnap(list_of_fasta_files=args.input, barrnap_out=args.output)
    sequence_16s = exid16s.extract_16S_barrnap(path=barrnap_result)
    exid16s.run_Kraken2(parent_dir=barrnap_result, path_16s=sequence_16s, database_file_path=args.database)



if __name__ == '__main__':
    main()