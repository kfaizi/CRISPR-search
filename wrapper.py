#!/Users/kianfaizi/miniconda3/bin/python
# remove shebang before export

# this script is customized (hardcoded?) to match my andromeda setup and so
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
from pathlib import Path


greeter = argparse.ArgumentParser(
    description='''For a given concatenated FASTA file, this script generates an
    indexed nucleotide BLASTDB, searches for similar seqs to a given query with
    TBLASTN, and searches those hits for DRs with CRISPRFinder. Results are
    tabulated in a CSV file.''')

greeter.add_argument("--input_fasta", help="Full path to the concatenated \
                    input .fasta file.", required=True)
greeter.add_argument("--db", help="Full path to the blastdb to be created.",
                     required=True)
greeter.add_argument("--query", help="Full path to file with known protein \
                     sequences for BLAST search.", required=True)
greeter.add_argument("--script", help="Full path to CRISPRFinder script.",
                     required=True)
greeter.add_argument("--threads", help="Number of threads for TBLASTN.",
                     required=True)
greeter.add_argument("--output", help="Full path to output directory.",
                     required=True)

hello = greeter.parse_args()

# cat *.fsa_nt > new_file.fsa


def make_db(path_to_fasta, path_to_new_database):
    """Concatenated fasta > parsed blastdb."""
    subprocess.run(["makeblastdb",
                    "-in",
                    f"{path_to_fasta}",  # in ~/genomes/
                    "-parse_seqids",
                    "-dbtype",
                    "nucl",
                    "-out",
                    f"{path_to_new_database}"],
                   cwd="/",
                   check=True)

# makeblastdb -in foo_cat.fsa -parse_seqids -dbtype nucl [-out bar]


def search_db(path_to_database, path_to_query, path_to_new_archive, threads):
    """blastdb + query > hits in ASN.1 format."""
    subprocess.run(["tblastn",
                    "-db",
                    f"{path_to_database}",  # feed arg from make_db()
                    "-query",
                    f"{path_to_query}",  # known protein reference
                    "-out",
                    f"{path_to_new_archive}",  # name, to differentiate results
                    "-evalue",
                    "10",  # could go bigger...
                    "-num_threads",
                    f"{threads}",
                    "-outfmt",  # need to also make output dir if not existing already?
                    "11"],  # ASN.1
                   cwd="/",
                   check=True)


def make_csv(path_to_archive, path_to_new_csv):
    """ASN.1 > csv with headers."""
    subprocess.run(["blast_formatter",
                    "-archive",
                    f"{path_to_archive}",  # path to ASN.1 from search_db()
                    "-outfmt"
                    "10",  # csv
                    "sseqid",  # contig ID
                    "sseq",  # aligned part of subject; the 'new orth seq'
                    # length of sseq
                    # length of the whole source contig
                    "qseqid",  # known cas13d that matched ('closest known')
                    "ppos",  # % similarity (% positives)
                    "pident",  # % identity
                    "qseq",  # aligned part of query; from the known cas13d
                    "sacc",  # accession of contig. Can I get bioproj/sample?
                    "sscinames",  # scientific name for contig?
                    "staxids",  # see above. need taxid db?
                    "bitscore",
                    "evalue",
                    "-out",
                    f"{path_to_new_csv}"],
                   cwd="/",
                   check=True)

# descriptive header. Is there one already?
# biosample/bioproject/genus/annotations/taxonomy information?

def read_csv(path_to_csv):  # gets accession IDs
    with open(path_to_csv, encoding='utf8') as f:
        ID_finder = csv.reader(f)
        output_list = []  # to be list of IDs
        next(ID_finder, None)  # skip the headers. Does this work w/ blast csv?
        try:
            for row in ID_finder:
                output_list.append(row[0])  # make list of IDs
            print("Working...", sum(1 for line in f), "IDs to parse")
            print("Done...")
            return output_list
        except csv.Error as e:
            sys.exit(f"file {path_to_csv}, line {reader.line_num}: {e}")


def list_csv(path_to_csv, path_to_new_IDs):  # writes to .txt (creates if DNE)
    with open(path_to_new_IDs, 'w+') as f:
        output_list = read_csv()
        print(f"Writing IDs to text file {path_to_new_IDs}...")
        for ID in output_list:
            f.write("{ID}\n")
        print("Done...")


def extract_contigs(path_to_database, path_to_IDs, path_to_new_contigs):
    """Extract contigs by referencing sseqids to blastdb."""
    subprocess.run(["blastdbcmd",
                    "-db",
                    f"{path_to_database}",  # parent
                    "-dbtype",
                    "nucl",
                    "-entry_batch",
                    f"{path_to_IDs}",  # list of sseqids to extract
                    "-out",
                    f"{path_to_new_contigs}"],
                   cwd="/",
                   check=True)

# put the FULL contig in the csv too.
# best to concatenate these hits for cf? Or pipe in one-at-a-time?

# must first move these final csvs to a new directory within the cf directory
# to avoid the blanking bug.


def find_CRISPRs(path_to_crisprfinder, path_to_contigs, path_to_new_results):
    """Search extracted contigs for DRs; retain output."""
    start = Path(f"{path_to_crisprfinder}")
    start = start.resolve()
    basedir = start.parent  # to make cwd same as cf's. see below
    p = subprocess.run(["perl",
                        f"{path_to_crisprfinder}",
                        f"{path_to_contigs}",
                        f"{path_to_new_results}"],
                       cwd=f"{basedir}",  # important, prevents bug
                       capture_output=True,
                       text=True,
                       check=True)
    print("Results: \n" + "-"*7 + "\n" + p.stdout)


def update_csv(path_to_csv, path_to_results):




# using the initial input paths across all blocks, rather than feeding
# returned output sequentially.

# exceptions: sys.exit vs print. When something breaks, does the script stop?
# or still attempt the next try block (bad)?
# review logic. "Exception as e" OK?


def wrapper():

    db_name = str((hello.db.split("/"))[-1])  # archive named in kind
    query_name = str((hello.query.split("/"))[-1])
    archive_path = Path(hello.output, db_name)
    csv_path = Path(archive_path).with_suffix(".csv")

    txt_name = db_name + "_acc"  # text file of ACCession numbers
    txt_path = Path(hello.output, txt_name).with_suffix(".txt")

    contigs_name = db_name + "_hits"  # extracted contigs
    contigs_path = Path(hello.output, "in", contigs_name).with_suffix(".fa")

    contigs_path = contigs_path.resolve()
    (contigs_path.parent).mkdir(parents=True,
                                exist_ok=True)  # make "in" directory, if DNE

    results_name = db_name + "_results"  # CRISPRFinder output
    results_path = Path(hello.output, results_name)

    try:  # making the blastdb
        print(f"Making blastdb {db_name}...")
        make_db(hello.input_fasta, hello.db)
        print(f"Success! Blast database created at {hello.db}")
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error making blastdb {db_name}: {e}")

    try:  # blasting the blastdb
        print(f"Blasting proteins in {query_name} against {db_name} with \
              {hello.threads} threads...")
        search_db(hello.db, hello.query, archive_path, hello.threads)
        print(f"Success! Blast results written to {archive_path}")
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error blasting {query_name} against blastdb {db_name}: {e}")

    try:  # creating a .csv from the results
        print("Creating a .csv file from the blast results...")
        make_csv(archive_path, csv_path)
        print(f"Success! New .csv file created at {csv_path}")
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error creating .csv file from the archive \
                 at {archive_path}: {e}")

    try:  # make a list of accession IDs. Need to verify this step (csv format)
        print("Extracting accession numbers from the .csv file...")
        list_csv(csv_path, txt_path)
        print(f"Success! Accession list created at {txt_path}")
    except Exception as e:
        sys.exit(f"Error creating accession list at {txt_path}: {e}")

    try:  # extract the corresponding contigs. TODO: ensure results are cat
        print(f"Reading accession list at {txt_path} and extracting \
              corresponding contigs from {hello.db}...")
        extract_contigs(hello.db, txt_path, contigs_path)
        print(f"Success! Contigs extracted to {contigs_path}")
    except Exception as e:
        sys.exit(f"Error extracting contigs to {contigs_path}: {e}")

    try:  # search 'em with crisprfinder!
        print(f"Almost there. Searching the sequences in {contigs_name} \
              for DRs...")
        find_CRISPRs(hello.script, contigs_path, results_path)
        print(f"Success! CRISPRFinder output created at {results_path}")
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error searching for DRs: {e}")

# finally, update the csv
