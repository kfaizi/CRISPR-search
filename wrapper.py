#!/Users/kianfaizi/miniconda3/bin/ python

import os
import sys
import subprocess
import re
import csv
import argparse
from sys import argv



# for a given NCBI nucl blastdb with seqids parsed, automates tblastn search >
# creates a master_csv of results > extracts hits > cat hits > return .fa

# cat *.fsa_nt > new_file.fsa

# def cat_and_move(input_loc, output_loc):
#     # need to check naming conventions
#     output_name = (output_loc.split("/"))[-1]
#     subprocess.run(["cat",
#                    "*.fa",
#                    ">",
#                    f"{output_name}"
#                    ])

def make_db(input_loc):
    """Concatenated fasta > parsed nucl blastdb."""
    subprocess.run(["makeblastdb",
                   "-in",
                   f"{input_loc}",
                   "-parse_seqids",
                   "-dbtype",
                   "nucl",
                   "-out",
                   f"{output_loc}"],
                   cwd="/",
                   check=True)

# makeblastdb -in foo_cat.fsa -parse_seqids -dbtype nucl [-out bar]

def search_db(db, output_loc, threads):
    """blastdb + query > hits in ASN.1"""
    subprocess.run(["tblastn",
                   "-db",
                   f"{db}",
                   "-query",
                   "/Users/kianfaizi/Desktop/Cas13d_proteins.fasta",
                   "-out",
                   f"{output_loc}",
                   "-evalue",
                   "10",
                   "-num_threads",
                   f"{threads}",
                   "-outfmt",  # need to also make output dir if not existing already?
                   "11",
                   "-show_gis"],  # necessary?
                   cwd="/",
                   check=True)


def make_master_csv(input_loc, output_loc):
    """ASN.1 > csv with headers"""
    subprocess.run(["blast_formatter",
                   "-archive",
                   f"{input_loc}",
                   "-outfmt"
                   "10",  # csv
                   "sseq",  # aligned part of subject (but, is it complete???); the 'new orth seq'
                   # length of this new ortholog
                   # compare to length of its contig
                   "qseqid",  # known cas13d that matched ('closest known')
                   "ppos",  # % similarity (% positives)
                   "pident",  # % identity
                   "sseqid",  # ID of contig containing hit
                   "qseq",  # aligned part of query (but, is it complete???); the 'contig seq'
                   "sacc",  # accession of contig. how to get bioproj/sample?
                   "sscinames",  # scientific name for contig; N/A?
                   "staxids",  # see above. need taxid db?
                   "bitscore",
                   "evalue",
                   "-out",
                   f"{output_loc}"],
                   cwd="/",
                   check=True)

# descriptive header
# need to extract full seq or do sseq/qseq do the job? find out! then obtain lengths.
# biosample/bioproject/genus/annotations/taxonomy information?

# stdout=subprocess.PIPE

def extract_hits(database_name, 'entry batch thing?'):
    """Extract contigs by referencing sseqids to blastdb"""
    subprocess.run(["blastdbcmd",
                   "-db",
                   f"{database_name}",
                   "-dbtype",
                   "nucl",
                   "-entry_batch",
                   f"{what to call this?}"],
                   cwd="/",
                   check=True)
    return True # how to redirect output? send to .csv!


def find_CRISPRs(): # need to install cf_v1 on andromeda
    """Search extracted contigs for DRs and retain output"""
    return True


def csv_updater(): # can I define a function to which I can pass info
# in order to easily update the csv? Or just use module csv
    return True



def wrapper(concatenated_raw_fasta):  # full path
    try:
        blast_database = make_db(concatenated_raw_fasta)  # track exit codes
        print("Making blastdb...")
        return ####
    except subprocess.CalledProcessError as e:
        sys.exit(e)

# or use if-statement logical structure with breaks for errors?

# def read_csv():
#     with open(args.input, encoding='utf8') as f:  # gets the accession IDs
#         ID_finder = csv.reader(f)
#         output_list = []  # to be list of IDs
#         # next(ID_finder, None)  # skip the headers (use if present)
#         try:
#             for row in ID_finder:
#                 output_list.append(row[0])  # make list of IDs
#             print("Working...", sum(1 for line in f), "IDs to parse")  # count IDs
#             print("Parsed...")
#             return output_list
#         except csv.Error as e:
#             sys.exit('file %s, line %d: %s' % (args.input, reader.line_num, e))


# search_db("/Users/kianfaizi/Desktop/S_search/S_cat.fsa", "/Users/kianfaizi/Desktop/results_S_search", "8")

# parser = argparse.ArgumentParser(description='This script accepts a .csv file,\
# gets accession numbers, and writes them to a .txt file.')

# parser.add_argument('-i', '--input', help='Input .csv filename', required=True)
# parser.add_argument('-o', '--output', help='Output .txt filename', required=False)
# args = parser.parse_args()



