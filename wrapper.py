"""A Python (3) wrapper for identifying Cas protein orthologs.

Searches prokaryotic/metagenomic DNA sequence data for
new orthologs of known Cas proteins. Strings together
analysis and output from several existing tools.

Accepts a single FASTA file.
Generates an indexed nucleotide BLASTDB, searches it for ORFs matching a
protein query using TBLASTN, and then checks hits for the presence of DRs
with CRISPRFinder. The results are tabulated in a CSV file.

This script is not well-optimized... yet ;)
"""
__author__ = "Kian Faizi"
__copyright__ = "Copyright 2019, Kian Faizi"
#  __license__ = null
# __version__ = "1.0"
__maintainer__ = "Kian Faizi"
__email__ = "kfaizi@ucsd.edu"
__status__ = "Prototype"

import argparse
import subprocess
import sys
import csv
# import pandas as pd

from pathlib import Path

# in working directory with CF + query, create blastdb/ and genomes/.
# genomes including cat fastas should be in genomes/.
# indexed dbs should be in blastdb/.

greeter = argparse.ArgumentParser(
    description='''For a given concatenated FASTA file, this script generates an
    indexed nucleotide BLASTDB, searches for similar seqs to a given query with
    TBLASTN, and searches those hits for DRs with CRISPRFinder. Results are
    tabulated in a CSV file.''')

greeter.add_argument("-i", "--input_name", help="Name/identifier of the input file.", required=True)

greeter.add_argument("-w", "--working_dir", help="Path to working directory.", required=True)

greeter.add_argument("-b", "--blastdb_dir", help="Path to blastdb/.", required=True)

greeter.add_argument("-g", "--genomes_dir", help="Path to genomes/.", required=True)

greeter.add_argument("-q", "--query", help="Name of protein query file.", required=True)

greeter.add_argument("-s", "--script", help="Name of CRISPRFinder script file.", required=True)

greeter.add_argument("-t", "--threads", help="Number of threads for TBLASTN.", required=True)

## for interactivity, use: ##
# hello = greeter.parse_args()

## for hardcoded testing, use: ##
hello = greeter.parse_args(['-i', 'A',
                            '-w', '/Users/kianfaizi/dev/',  # use pathlib to avoid / errors
                            '-b', '/Users/kianfaizi/dev/blastdb/',  # use pathlib to avoid / errors
                            '-g', '/Users/kianfaizi/dev/genomes/',  # use pathlib to avoid / errors
                            '-q', 'Cas13d_proteins.fa',
                            '-s', 'cf_v2.pl',
                            '-t', '4'])

# make all pathlike args full paths (resolve "~/")
working_path = Path(hello.working_dir).expanduser()
genomes_dir = Path(hello.genomes_dir).expanduser()
blastdb_dir = Path(hello.blastdb_dir).expanduser()

# Constants
name = hello.input_name

blastdb_path = Path(blastdb_dir, name)
genomes_path = Path(genomes_dir, name)
query_path = sorted(working_path.glob(hello.query))[0]
og_script_path = sorted(working_path.glob(hello.script))[0]

output_path = Path(working_path, "output", name)  # for each name
output_path.mkdir(parents=True, exist_ok=True)  # make dirs incl parent
script_path = Path(output_path, hello.script)

genome_name = name + "_cat"
genome_path = Path(genomes_dir, genome_name).with_suffix(".fa")

archive_name = name + "_archive"
archive_path = Path(output_path, archive_name)

csv_name = name + "_csv"
csv_path = Path(output_path, csv_name).with_suffix(".csv")

accessions_name = name + "_acc"
accessions_path = Path(output_path, accessions_name).with_suffix(".txt")

contigs_name = name + "_hits"
contigs_path = Path(output_path, "in", contigs_name).with_suffix(".fa").resolve()
(contigs_path.parent).mkdir(parents=True, exist_ok=True)

results_name = name + "_crisprs"
results_path = Path(output_path, results_name)

# sections (columns) for the .csv file:
sections = "10 sseqid slen sstart send sseq length sframe qseqid ppos pident qseq bitscore evalue saccver"
sections_list = sections.split()
sections_list.pop(0)  # remove filetype specifier (10)

# Currently:
# sseqid, the ID for the contig/subject (the fasta header)
# slen, the total length of contig/subject ???
# sstart, start of alignment in contig
# send, end of alignment in contig
# sseq, the aligned part of the contig sequence
# length, the length of the aligned part
# sframe, the frame of subject
# qseqid, the name of the match (query/CasRx, Dj, etc)
# ppos, percent positive/sim
# pident, percent identity
# qseq, the aligned part of the match
# bitscore
# evalue
# saccver, accession ID.version

# Need:
# DR? y/p/n
# DR count
# Consensus DR?
# DR 1... 2... 3...
# A given contig may contain multiple arrays, need to be able to handle this


def move_script(path_to_script, path_to_copy):
    """To avoid CF errors, make a copy of it in each output directory."""
    subprocess.run(["cp",
                    f"{path_to_script}",
                    f"{path_to_copy}"],
                   cwd="/",
                   check=True)


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
    """Blastdb + query > hits in ASN.1 format."""
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
                    "-outfmt",  # need to make dir if not existing already?
                    "11"],  # ASN.1
                   cwd="/",
                   check=True)


def make_csv(path_to_archive, path_to_new_csv):
    """ASN.1 > csv with headers."""
    subprocess.run(["blast_formatter",
                    "-archive",
                    f"{path_to_archive}",  # path to ASN.1 from search_db()
                    "-outfmt",
                    sections,
                    "-out",
                    f"{path_to_new_csv}"],
                   cwd="/",
                   check=True)


def read_csv(path_to_csv):
    """Get accession IDs from csv."""
    with open(path_to_csv, newline='') as f:
        ID_finder = csv.reader(f)
        output_list = []  # to be list of IDs
        # next(ID_finder, None)  # skip the headers. Do not use for raw blast csv as it has none
        try:
            linecount = sum(1 for line in ID_finder)  # this iterates through the csv...
            print("Working. There are", linecount, "contig IDs to parse...")
            f.seek(0)  # ...so we need to move back to start
            ID_finder = csv.reader(f)  # reset generator
            for line in ID_finder:
                output_list.append(line[0])  # make list of IDs
            return output_list
        except csv.Error as e:
            sys.exit(f"Error gathering accession data for {path_to_csv}, line {ID_finder.line_num}: {e}")
        finally:
            if len(output_list) != linecount:
                print("Error: initial and final line counts do not match! Exiting")
                sys.exit()


def list_csv(path_to_csv, path_to_new_IDs):
    """Write accession IDs to .txt, creating one if it doesn't exist."""
    with open(path_to_new_IDs, 'w+', newline='') as f:
        try:
            output_list = read_csv(path_to_csv)
            print(f"Writing IDs to text file {path_to_new_IDs}...")
            for ID in output_list:
                f.write(f"{ID}\n")
            # print("Done!")
        except csv.Error as e:
            sys.exit(f"Error writing accession data for {path_to_csv}: {e}")


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


def find_CRISPRs(path_to_crisprfinder, path_to_contigs, path_to_new_gff):
    """Search extracted contigs for DRs; retain output."""
    # make cwd same as CF's to prevent bug
    p = subprocess.run(["perl",
                        f"{path_to_crisprfinder}",
                        f"{path_to_contigs}",
                        f"{path_to_new_gff}"],
                       cwd=output_path,  # prevent bug
                       capture_output=True,
                       text=True,
                       check=True)
    print("\nResults: \n" + "-"*7 + "\n" + p.stdout)


def parseGffToDict(path_to_gff):  # ie results_path
    """Take gff filename, return parsed dict. Adapted from PDH."""
    # dict will be populated with crispr array info from gff
    gff_file = sorted(path_to_gff.glob("*.gff"))[0]  # find actual gff file
    results_dict = {'crisprs_list': []}
    with open(gff_file) as file:
        data = file.read()  # entire gff file
    data_split = data.split('\n\n')  # split into found crisprs
    data_split = list(filter(None, data_split))  # remove empty elements
    num_crisprs = len(list(data_split))
    results_dict.update({
        'num_crisprs': num_crisprs,
        # 'assembly_accession': assembly_accession
    })
    for section in data_split:  # for each found crispr section
        # split into lines and remove top '##gff-version 3'
        section_split = section.split('\n')
        section_split = [line for line in section_split if '##' not in line]
        for line_idx, line in enumerate(section_split):  # for each line in section
            line_split = line.split('\t')
            if line_idx == 0:  # first line is summary of found crispr
                # from last element in line (attributes), get values
                attr_split = line_split[-1].split(';')

                dr_consensus_seq = attr_split[0][len('DR='):]
                dr_consensus_length = attr_split[1][len('DR_length='):]
                crispr_id = attr_split[-1][len('ID='):]  # ie 'parent'

                # construct item to add to dictionary
                dict_item = {
                  'crispr_id': crispr_id,
                  'genomic_accession': line_split[0],
                  'crispr_type': line_split[2],
                  'start': line_split[3],
                  'end': line_split[4],
                  'dr_consensus_seq': dr_consensus_seq,
                  'dr_consensus_length': dr_consensus_length,
                  'elements': []
                }

            else:  # other lines are found crispr elements (dr or spacer)
                attr_split = line_split[-1].split(';')

                element_type = line_split[2]
                element_start = line_split[3]
                element_end = line_split[4]

                if element_type == 'CRISPRdr':
                    dict_item['elements'].append({
                      'element_type': element_type,
                      'start': element_start,
                      'end': element_end
                    })

                elif element_type == 'CRISPRspacer':
                    spacer_seq = attr_split[0][len('sequence='):]
                    spacer_name = attr_split[1][len('name='):]
                    dict_item['elements'].append({
                      'element_type': element_type,
                      'start': element_start,
                      'end': element_end,
                      'spacer_seq': spacer_seq,
                      'spacer_name': spacer_name
                    })

        # add to crisprs list in dictionary
        results_dict['crisprs_list'].append(dict_item)

    return gff_file, results_dict


# def update_csv(parsed_dict, path_to_csv, path_to_new_csv):
#     """Add headers, and append CF data."""
#     pass


def wrapper():
    """Put it all together."""
    try:  # copy CF to each result dir to avoid error
        move_script(og_script_path, script_path)
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error initializing file structure: {e}")

    # future: add a concatenation step? cat *.fsa_nt > new_file.fsa

    try:  # making the blastdb
        print(f"\nMaking blastdb {name}...")
        make_db(genome_path, blastdb_path)
        print(f"Success! Blast database created at {blastdb_path}")
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error making blastdb {name}: {e}")

    try:  # blasting the blastdb
        print(f"\nBlasting proteins in {query_path} against {name} with {hello.threads} threads...")
        search_db(blastdb_path, query_path, archive_path, hello.threads)
        print(f"Success! Blast results written to {archive_path}")
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error blasting {query_path} against blastdb {name}: {e}")

    try:  # creating a .csv from the results
        print("\nCreating .csv file from the blast results...")
        make_csv(archive_path, csv_path)
        print(f"Success! New .csv file created at {csv_path}")
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error creating .csv file from the archive at {archive_path}: {e}")

    try:  # make a list of accession IDs
        print("\nExtracting accession ID numbers from the .csv file...")
        list_csv(csv_path, accessions_path)
        print(f"Success! Accession list created at {accessions_path}")
    except Exception as e:
        sys.exit(f"Error creating accession list at {accessions_path}: {e}")

    try:  # extract the corresponding contigs
        print(f"\nReading accession list at {accessions_path} and extracting corresponding contigs from {blastdb_path}...")
        extract_contigs(blastdb_path, accessions_path, contigs_path)
        print(f"Success! Contigs extracted to {contigs_path}")
    except Exception as e:
        sys.exit(f"Error extracting contigs to {contigs_path}: {e}")

    try:  # search 'em with crisprfinder!
        print(f"\nSearching sequences in {contigs_path} for DRs...")
        find_CRISPRs(script_path, contigs_path, results_path)
        print(f"Success! CRISPRFinder output created at {results_path}")
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error searching for DRs: {e}")

    try:  # parse the gff output
        print(f"\nParsing the CRISPRFinder output...")
        gff_file, results_dict = parseGffToDict(results_path)
        print(f"Success! Parsed gff {gff_file} to dict")
        return results_dict
    except Exception as e:
        sys.exit(f"Error parsing results at {results_path}: {e}")

## work ##

    try:  # update the csv
        print(f"\nUpdating the .csv at {csv_path}...")
        print("Updating headers...")
        # function call for header remake
        print("Adding CRISPRFinder results...")
        # function call for update
        print(f"Success! Final output file located at {csv_path}")
    except Exception as e:
        sys.exit(f"Error updating .csv at {csv_path}: {e}")
