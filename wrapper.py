#!/Users/kianfaizi/miniconda3/bin/python
# remove shebang before export

# this script is customized (hard-coded?) to match my andromeda setup and so
# requires some care before use.

# in a/home directory, create blastdb/, genomes/, and output/.
# crisprfinder, and the query file, etc should be in output/.
# genomes including cat fastas should be in genomes/.
# indexed dbs only should be in blastdb/.
# the script will work from that assumption forwards.

import argparse
import subprocess
import sys
import csv
import re
import os


greeter = argparse.ArgumentParser(
    description='''For a given concatenated FASTA file, this script generates an
    indexed nucleotide BLASTDB, searches for similar seqs to a given query with
    TBLASTN, and searches those hits for DRs with CRISPRFinder. Results are
    tabulated in a CSV file.''')

greeter.add_argument("--input_fasta", help="Full path to the concatenated \
                    input .fasta file.", required=True)
greeter.add_argument("--db", help="Desired path to the blastdb to be created.",
                     required=True)
greeter.add_argument("--query", help="Full path to file with known protein \
                     sequences for BLAST search.", required=True)
greeter.add_argument("--script", help="Full path to CRISPRFinder script.",
                     required=True)
greeter.add_argument("--threads", help="Number of threads for TBLASTN.",
                     required=True)
greeter.add_argument("--archive", help="Desired path to the ASN.1 archive \
                     to be created from the results of TBLASTN.", required=True)
hello = greeter.parse_args()

# cat *.fsa_nt > new_file.fsa

def make_db(fasta_path, new_database_path):
    """Concatenated fasta > parsed blastdb."""
    subprocess.run(["makeblastdb",
               "-in",
               f"{fasta_path}",  # in ~/genomes/
               "-parse_seqids",
               "-dbtype",
               "nucl",
               "-out",
               f"{new_database_path}"],
               cwd="/",
               check=True)
    return new_database_path

# makeblastdb -in foo_cat.fsa -parse_seqids -dbtype nucl [-out bar]

def search_db(database_path, query_path, new_archive_path, threads):
    """blastdb + query > hits in ASN.1"""
    subprocess.run(["tblastn",
                   "-db",
                   f"{database_path}",  # feed arg from make_db()
                   "-query",
                   f"{query_path}",  # known protein reference
                   "-out",
                   f"{new_archive_path}",  # name, to differentiate results
                   "-evalue",
                   "10",
                   "-num_threads",
                   f"{threads}",
                   "-outfmt",  # need to also make output dir if not existing already?
                   "11"],  # ASN.1
                   cwd="/",
                   check=True)
    return new_archive_path


def make_master_csv(archive_path, new_csv_path):
    """ASN.1 > csv with headers"""
    subprocess.run(["blast_formatter",
                   "-archive",
                   f"{archive_path}",  # path to ASN.1 from search_db()
                   "-outfmt"
                   "10",  # csv
                   "sseq",  # aligned part of subject; the 'new orth seq'
                   # length of this new ortholog
                   # compare to length of its contig
                   "qseqid",  # known cas13d that matched ('closest known')
                   "ppos",  # % similarity (% positives)
                   "pident",  # % identity
                   "sseqid",  # ID of contig containing hit
                   "qseq",  # aligned part of query; from the known cas13d
                   "sacc",  # accession of contig. Can I get bioproj/sample?
                   "sscinames",  # scientific name for contig?
                   "staxids",  # see above. need taxid db?
                   "bitscore",
                   "evalue",
                   "-out",
                   f"{new_csv_path}"],  # should probably use same dir? need standardization
                   cwd="/",
                   check=True)

# descriptive header
# biosample/bioproject/genus/annotations/taxonomy information?
# stdout=subprocess.PIPE
# Extract sseqids for the next step:::

def extract_hits(database_path, IDs_path, new_contigs_path):
    """Extract contigs by referencing sseqids to blastdb"""
    subprocess.run(["blastdbcmd",
                   "-db",
                   f"{database_path}",  # parent
                   "-dbtype",
                   "nucl",
                   "-entry_batch",
                   f"{IDs_path}",  # list of sseqids to extract
                   "-out",
                   f"{contigs_path}"],  #
                   cwd="/",
                   check=True)

# put the FULl contig in the csv too.
# best to concatenate these hits for CF? Or pipe in one-at-a-time?

# must first move these final csvs to a new directory within the cf_v1 directory
# to avoid the blanking bug.

def find_CRISPRs(crisprfinder_path, contigs_path, new_CRISPRs_path):
    """Search extracted contigs for DRs; retain output"""
    subprocess.run(["perl",
               f"{crisprfinder_path}",
               f"{contigs_path}",
               f"{new_CRISPRs_path}"],
               cwd="/",
               check=True)
    return True


def csv_updater(csv_path, CRISPRs_path):
    return True


def wrapper():
    try:  # making the blastdb
        db_name = (hello.db.split("/"))[-1]
        print(f"Making blastdb {db_name} ...")
        blast_database = make_db(hello.input_fasta, hello.db)  # captures new database path
        print("Success!")
        return db_name, blast_database
    except subprocess.CalledProcessError as e:
        print(f"Error making blastdb {db_name}")
        sys.exit(e)
    try:  # blasting the blastdb
        query_name = (hello.query.split("/"))[-1]
        archive_name = (hello.archive.split("/"))[-1]
        print(f"Blasting proteins in the file {query_name} against database {db_name}, \
              using {hello.threads} threads...")
        blast_results = search_db(blast_database, hello.query, hello.archive, hello.threads)  # Better way to specify archive destination? Convention.
        print("Success!")
        print(f"Blast results written to {archive_name}.")
        return blast_results
    except subprocess.CalledProcessError as e:
        print(f"Error blasting {query_name} against blastdb {db_name}")
        sys.exit(e)
    try:  # creating a .csv from the results


