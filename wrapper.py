#DB DIFF

"""A Python (3.7) wrapper for identifying Cas protein orthologs.

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
import pandas as pd
# import logging
# import logging.config

from pathlib import Path
from itertools import groupby
from Bio import SeqIO

# in working directory with CF + query, create blastdb/ and genomes/.
# genomes including cat fastas should be in genomes/.
# indexed dbs should be in blastdb/.

greeter = argparse.ArgumentParser(
    description='''For a given concatenated FASTA file, this script generates an
    indexed nucleotide BLASTDB, searches for similar seqs to a given query with
    TBLASTN, and searches those hits for DRs with CRISPRFinder. Results are
    tabulated in a CSV file.''')

greeter.add_argument("-in", "--input_name", help="Name/identifier of the input file.", required=True)

greeter.add_argument("-wdir", "--working_dir", help="Path to working directory.", required=True)

greeter.add_argument("-bdir", "--blastdb_dir", help="Path to blastdb/.", required=True)

greeter.add_argument("-gdir", "--genomes_dir", help="Path to genomes/.", required=True)

greeter.add_argument("-q", "--query", help="Name of protein query file.", required=True)

greeter.add_argument("-s", "--script", help="Name of CRISPRFinder script file.", required=True)

greeter.add_argument("-t", "--threads", help="Number of threads for TBLASTN.", required=True)

######################## for interactivity, use: ###########################
# hello = greeter.parse_args()

######################## for hardcoded testing, use: ########################
# hello = greeter.parse_args(['-i', 'A',
#                             '-w', '/Users/kianfaizi/dev/',  # use pathlib to avoid / errors
#                             '-b', '/Users/kianfaizi/dev/blastdb/',  # use pathlib to avoid / errors
#                             '-g', '/Users/kianfaizi/dev/genomes/',  # use pathlib to avoid / errors
#                             '-q', 'Cas13d_proteins.fa',
#                             '-s', 'cf_v2.pl',
#                             '-t', '4'])

hello = greeter.parse_args(['-i', 'V',
                            '-w', '~/search/',  # use pathlib to avoid / errors
                            '-b', '~/search/blastdb/',  # use pathlib to avoid / errors
                            '-g', '~/search/genomes/',
                            '-q', 'Cas13d_proteins.fa',
                            '-s', 'cf_v2.pl',
                            '-t', '4'])

# configure logging
# logging.config.dictConfig('logging.yaml')
# logger = logging.getLogger('statusLog')

NCBI_datestamp = "jul25"  # track date of last database update (NCBI)

# make all pathlike args full paths (resolve "~/")
working_path = Path(hello.working_dir).expanduser()
genomes_dir = Path(hello.genomes_dir).expanduser()
blastdb_dir = Path(hello.blastdb_dir).expanduser()

# Constants
name = hello.input_name + "_" + NCBI_datestamp

blastdb_path = Path(blastdb_dir, hello.input_name, name)
# genomes_path = Path(genomes_dir, name)
query_path = sorted(working_path.glob(hello.query))[0]
og_script_path = sorted(working_path.glob(hello.script))[0]

output_path = Path(working_path, "output", name)  # for each name
output_path.mkdir(parents=True, exist_ok=True)  # make dirs incl parent
script_path = Path(output_path, hello.script)

genome_path = Path(genomes_dir, hello.input_name, name).with_suffix(".fa")

archive_name = name + "_archive"
archive_path = Path(output_path, archive_name)

csv_name = name + "_csv"
csv_path = Path(output_path, csv_name).with_suffix(".csv")

accessions_name = name + "_acc"
accessions_path = Path(output_path, accessions_name).with_suffix(".txt")

contigs_name = name + "_hits"
contigs_path = Path(output_path, "in", contigs_name).with_suffix(".fa").resolve()
(contigs_path.parent).mkdir(parents=True, exist_ok=True)

dedupe_name = name + "_dedupe"
dedupe_path = Path(contigs_path.parent, dedupe_name).with_suffix(".fa")

results_name = name + "_crisprs"
results_path = Path(output_path, results_name)

summary_csv_name = name + "_summary_csv"
summary_csv_path = Path(output_path, summary_csv_name).with_suffix(".csv")

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


def move_script(path_to_script, path_to_copy):
    """To avoid CF errors, make a copy of it in each output directory."""
    subprocess.run(["cp",
                    f"{path_to_script}",
                    f"{path_to_copy}"],
                   cwd="/",
                   check=True)

# def unzip_fasta(path_to_zipped):
#   """(G)unzip compressed fasta files."""
#   # from '~/search/raw/foo/ to ~/search/genomes/foo/'


# def cat_and_trim(path_to_fastas):
#   """Concatenate unzipped fastas, delete source."""
#   # in ~/search/genomes/foo/, delete all of type "wgs.*", keep the datestamp

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


def dedupe_contigs(path_to_contigs, path_to_new_contigs):
    """Dedupe DNA fasta by sequence; from https://www.biostars.org/p/3003/."""
    all_seqs = set()
    with open(path_to_contigs, 'r') as infile:
        with open(path_to_new_contigs, 'w+') as outfile:
            head = None
            for h, lines in groupby(infile, lambda x: x.startswith('>')):
                if h:
                    head = next(lines)
                else:
                    seq = ''.join(lines)
                    if seq not in all_seqs:
                        all_seqs.add(seq)
                        outfile.write(f"{head}{seq}")


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


def update_csv(results_dict, path_to_blast_csv, path_to_new_csv):
    """Add headers, and append CF data."""
    # make blast csv into dataframe
    df_blast = pd.read_csv(path_to_blast_csv, dtype=object, header=None, names=sections_list, index_col=None)

    # make parsed crisprfinder results into dataframe
    d_crispr = results_dict['crisprs_list']
    df_crispr = pd.DataFrame(d_crispr, dtype=object, columns=['genomic_accession',
                                                              'crispr_id',
                                                              'crispr_type',
                                                              'start',
                                                              'end',
                                                              'dr_consensus_seq',
                                                              'dr_consensus_length',
                                                              'elements'], index=None)

    # make master dataframe. creates duplicates as needed to
    # handle multiple hits from blast/crisprfinder/both
    combo = pd.merge(df_blast, df_crispr, how='left', left_on='saccver', right_on='genomic_accession')
    combo = combo.drop(columns='genomic_accession')  # redundant; use crispr_id

    # make new summary dataframe
    summary = combo[['sseq',
                     'qseqid',
                     'ppos',
                     'pident',
                     'crispr_type',
                     'crispr_id',
                     'dr_consensus_seq',
                     'dr_consensus_length',
                     'bitscore',
                     'evalue',
                     'slen',
                     'sstart',
                     'send',
                     'start',
                     'end',
                     'elements']].copy()

    # create columns w/ whether DR found (y/n) and if yes, how many
    dr_status = []
    dr_sum = []

    for eltup in summary[['dr_consensus_seq', 'elements']].itertuples():
        dr_count = 0
        dr_exists = 'no'

        # if dr_consensus_seq is there, then DR was found
        if isinstance(eltup.dr_consensus_seq, str):
            dr_exists = 'yes'

        dr_status.append(dr_exists)

        try:
            for dictx in eltup.elements:
                for v in dictx.values():
                    if v == "CRISPRdr":
                        dr_count += 1
        except Exception:
            pass

        dr_sum.append(dr_count)

    # add columns to summary
    summary.loc[:, 'DR_found'] = dr_status
    summary.loc[:, 'num_DRs'] = dr_sum

    # reorder columns for readability
    summary = summary[['sseq',
                       'qseqid',
                       'ppos',
                       'pident',
                       'bitscore',
                       'evalue',
                       'crispr_type',
                       'crispr_id',
                       'DR_found',
                       'num_DRs',
                       'dr_consensus_seq',
                       'dr_consensus_length',
                       'slen',
                       'sstart',
                       'send',
                       'start',
                       'end',
                       'elements']]

    # rename columns more intuitively
    summary = summary.rename(columns={'sseq': 'protein_sequence',
                                      'slen': 'contig_length',
                                      'qseqid': 'ortholog_match',
                                      'dr_consensus_seq': 'DR_consensus_seq',
                                      'dr_consensus_length': 'DR_consensus_length',
                                      'sstart': 'protein_start',
                                      'send': 'protein_end',
                                      'start': 'array_start',
                                      'end': 'array_end',
                                      'elements': 'array'})

    # parse contigs for subsequent extraction, without overloading memory
    record_dict = SeqIO.index(str(dedupe_path), 'fasta')  # w/o str, breaks

    # extract up to +/- 20kb flanking sequence around crispr locus
    extracted_heads = []
    extracted_paths = []

    for row in summary.itertuples():
        extracted_head = None
        extracted_path = None

        try:
            acc = row.crispr_id.split('_')[0]
            crispr_id = row.crispr_id
            start = int(row.array_start)
            end = int(row.array_end)
            length = int(row.contig_length)

            if acc in record_dict:  # then contig has crispr
                extracted_start = max(start-20000, 1)
                extracted_end = min(end+20000, length)
                extracted_len = extracted_end - extracted_start + 1
                extracted_seq = (record_dict[acc])[extracted_start-1: extracted_end]

                # print('For', acc, 'start at', extracted_start, 'and end at', extracted_end, '. The seq is', length, 'long, of which we\'re taking', extracted_len)
                extracted_head = str(record_dict[acc].description)
                # extracted_body = str(record_dict[acc].seq)
                extracted_name = crispr_id + "_extracted.fa"
                extracted_path = Path(results_path, extracted_name)

                # write extracted sequence to a new fasta file...
                # with open(extracted_path, 'w') as outfile:
                #     SeqIO.write(extracted_seq, outfile, 'fasta')

                # print(f"Wrote extracted CRISPR locus to {extracted_path}")

        except Exception:
            pass

        extracted_heads.append(extracted_head)
        extracted_paths.append(extracted_path)

    # ...and add the header and filename to the summary dataframe
    summary.loc[:, 'extracted_header'] = extracted_heads
    summary.loc[:, 'extracted_filename'] = extracted_paths

    # write summary dataframe to csv
    with open(summary_csv_path, 'w') as f:
        summary.to_csv(f, mode='w', index=False, header=True)


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

    try:  # deduplicate the contigs fasta
        print(f"\nDeduplicating BLAST hits at {contigs_path}...")
        dedupe_contigs(contigs_path, dedupe_path)
        print(f"Success! Contigs deduplicated and moved to {dedupe_path}")
    except Exception as e:
        sys.exit(f"Error deduplicating contigs at {contigs_path}: {e}")

    try:  # search 'em with crisprfinder
        print(f"\nSearching sequences in {dedupe_path} for DRs...")
        find_CRISPRs(script_path, dedupe_path, results_path)
        print(f"Success! CRISPRFinder output created at {results_path}")
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error searching for DRs: {e}")

    try:  # parse the gff output
        print(f"\nParsing the CRISPRFinder output...")
        gff_file, results_dict = parseGffToDict(results_path)
        print(f"Success! Parsed gff {gff_file} to dict")
    except Exception as e:
        sys.exit(f"Error parsing results at {results_path}: {e}")

    try:  # update the csv
        print(f"\nUpdating the .csv at {csv_path} with new info...")
        update_csv(results_dict, csv_path, summary_csv_path)
        print(f"Success! Final output file located at {summary_csv_path}")
    except Exception as e:
        sys.exit(f"Error updating .csv at {csv_path}: {e}")

# results_dict = {'crisprs_list': [{'crispr_id': 'AAAA02021996.1_21523_21912', 'genomic_accession': 'AAAA02021996.1', 'crispr_type': 'CRISPR', 'start': '21523', 'end': '21912', 'dr_consensus_seq': 'TTCTGTATATTTTCGGACTTGTCCGAAA', 'dr_consensus_length': '28', 'elements': [{'element_type': 'CRISPRdr', 'start': '21523', 'end': '21550'}, {'element_type': 'CRISPRspacer', 'start': '21551', 'end': '21589', 'spacer_seq': 'CTCATTTCCGGACTTTCCGAAAACACATAGAACCAGATT', 'spacer_name': 'spacer_21551_39'}, {'element_type': 'CRISPRdr', 'start': '21590', 'end': '21617'}, {'element_type': 'CRISPRspacer', 'start': '21618', 'end': '21661', 'spacer_seq': 'GTGATTTTCGGACTTTCCGAAAAAGACTGCGAAGGCAAAAGTGG', 'spacer_name': 'spacer_21618_44'}, {'element_type': 'CRISPRdr', 'start': '21662', 'end': '21689'}, {'element_type': 'CRISPRspacer', 'start': '21690', 'end': '21734', 'spacer_seq': 'CTGATTTTCGGATTTTCCGAAAATCATCAGTAGAGTCAATTTCGC', 'spacer_name': 'spacer_21690_45'}, {'element_type': 'CRISPRdr', 'start': '21735', 'end': '21762'}, {'element_type': 'CRISPRspacer', 'start': '21763', 'end': '21807', 'spacer_seq': 'GTTATTTTCGGACTTTTCCGAGAACATCCAGAAGGATGTTGTTGG', 'spacer_name': 'spacer_21763_45'}, {'element_type': 'CRISPRdr', 'start': '21808', 'end': '21835'}, {'element_type': 'CRISPRspacer', 'start': '21836', 'end': '21884', 'spacer_seq': 'AGTTTTTGGAATATCCGAAAAATCCTTCGTTGACTTTGCTGTTGGTGCT', 'spacer_name': 'spacer_21836_49'}, {'element_type': 'CRISPRdr', 'start': '21885', 'end': '21912'}]},{'crispr_id': 'AAAA02027088.1_723_794', 'genomic_accession': 'AAAA02027088.1', 'crispr_type': 'PossibleCRISPR', 'start': '723', 'end': '794', 'dr_consensus_seq': 'TATCAATTTGACTCTAATTTTTT', 'dr_consensus_length': '23', 'elements': [{'element_type': 'CRISPRdr', 'start': '723', 'end': '745'}, {'element_type': 'CRISPRspacer', 'start': '746', 'end': '771', 'spacer_seq': 'TTAGTTTACATATTTTAATATCTAAC', 'spacer_name': 'spacer_746_26'}, {'element_type': 'CRISPRdr', 'start': '772', 'end': '794'}]}, {'crispr_id': 'AAAA02027088.1_13954_14066', 'genomic_accession': 'AAAA02027088.1', 'crispr_type': 'PossibleCRISPR', 'start': '13954', 'end': '14066', 'dr_consensus_seq': 'AAGCAAATAAAGCACGCTCTCAAATTT', 'dr_consensus_length': '27', 'elements': [{'element_type': 'CRISPRdr', 'start': '13954', 'end': '13980'}, {'element_type': 'CRISPRspacer', 'start': '13981', 'end': '14039', 'spacer_seq': 'CAGCAATTGATGACAAACTATTTGAGATATTATCTTTTGGTAAGGTTTAAATTCAGTCT', 'spacer_name': 'spacer_13981_59'}, {'element_type': 'CRISPRdr', 'start': '14040', 'end': '14066'}]}, {'crispr_id': 'AAAA02022349.1_23281_23389', 'genomic_accession': 'AAAA02022349.1', 'crispr_type': 'PossibleCRISPR', 'start': '23281', 'end': '23389', 'dr_consensus_seq': 'TGATGTGATGGAAAGTTGAAAGTTTGGA', 'dr_consensus_length': '28', 'elements': [{'element_type': 'CRISPRdr', 'start': '23281', 'end': '23308'}, {'element_type': 'CRISPRspacer', 'start': '23309', 'end': '23361', 'spacer_seq': 'AAAAAAACTTTGGAACTAAATAGGGCCTGTGTAGGAAAGTTTTGGATGTGATA', 'spacer_name': 'spacer_23309_53'}, {'element_type': 'CRISPRdr', 'start': '23362', 'end': '23389'}]}, {'crispr_id': 'AAAA02006613.1_10366_10470', 'genomic_accession': 'AAAA02006613.1', 'crispr_type': 'PossibleCRISPR', 'start': '10366', 'end': '10470', 'dr_consensus_seq': 'ATGAAATCCATTTTTACCAAACCA', 'dr_consensus_length': '24', 'elements': [{'element_type': 'CRISPRdr', 'start': '10366', 'end': '10389'}, {'element_type': 'CRISPRspacer', 'start': '10390', 'end': '10446', 'spacer_seq': 'TTTTAATTATTACCATGAAATCCCTGCGGTACGGAGCTCGATTCTCTCTGGGTTGCC', 'spacer_name': 'spacer_10390_57'}, {'element_type': 'CRISPRdr', 'start': '10447', 'end': '10470'}]}, {'crispr_id': 'AAAA02026716.1_8899_8996', 'genomic_accession': 'AAAA02026716.1', 'crispr_type': 'PossibleCRISPR', 'start': '8899', 'end': '8996', 'dr_consensus_seq': 'CGCCGAGCCGCCGCCGTAGCCGT', 'dr_consensus_length': '23', 'elements': [{'element_type': 'CRISPRdr', 'start': '8899', 'end': '8921'}, {'element_type': 'CRISPRspacer', 'start': '8922', 'end': '8973', 'spacer_seq': 'GCCGCCCGCCGACGCCATCGACGCCGCTACGGAGCAGGGGCCGTAGCTGCTA', 'spacer_name': 'spacer_8922_52'}, {'element_type': 'CRISPRdr', 'start': '8974', 'end': '8996'}]}, {'crispr_id': 'AAAA02005981.1_22110_22187', 'genomic_accession': 'AAAA02005981.1', 'crispr_type': 'PossibleCRISPR', 'start': '22110', 'end': '22187', 'dr_consensus_seq': 'TTACTACTATAGTACAAATGCATTT', 'dr_consensus_length': '25', 'elements': [{'element_type': 'CRISPRdr', 'start': '22110', 'end': '22134'}, {'element_type': 'CRISPRspacer', 'start': '22135', 'end': '22162', 'spacer_seq': 'ATCTATTTAACTTTATGCATTGTTCTCT', 'spacer_name': 'spacer_22135_28'}, {'element_type': 'CRISPRdr', 'start': '22163', 'end': '22187'}]}, {'crispr_id': 'AAAA02005981.1_27204_27318', 'genomic_accession': 'AAAA02005981.1', 'crispr_type': 'PossibleCRISPR', 'start': '27204', 'end': '27318', 'dr_consensus_seq': 'GGAGTCGGACTTTTTTATTTTTT', 'dr_consensus_length': '23', 'elements': [{'element_type': 'CRISPRdr', 'start': '27204', 'end': '27226'}, {'element_type': 'CRISPRspacer', 'start': '27227', 'end': '27244', 'spacer_seq': 'TAAAAGGGAATGGATATA', 'spacer_name': 'spacer_27227_18'}, {'element_type': 'CRISPRdr', 'start': '27245', 'end': '27267'}, {'element_type': 'CRISPRspacer', 'start': '27268', 'end': '27295', 'spacer_seq': 'ATTTTTCGCTATTGTTTTTGTTTTGACG', 'spacer_name': 'spacer_27268_28'}, {'element_type': 'CRISPRdr', 'start': '27296', 'end': '27318'}]}, {'crispr_id': 'AAAA02026886.1_5066_5145', 'genomic_accession': 'AAAA02026886.1', 'crispr_type': 'PossibleCRISPR', 'start': '5066', 'end': '5145', 'dr_consensus_seq': 'TATTGCTTTTATTATTATGTCTATAT', 'dr_consensus_length': '26', 'elements': [{'element_type': 'CRISPRdr', 'start': '5066', 'end': '5091'}, {'element_type': 'CRISPRspacer', 'start': '5092', 'end': '5119', 'spacer_seq': 'ATAAATAGGAGACACTATGCAGTGTAAA', 'spacer_name': 'spacer_5092_28'}, {'element_type': 'CRISPRdr', 'start': '5120', 'end': '5145'}]}, {'crispr_id': 'AAAA02016417.1_14766_14847', 'genomic_accession': 'AAAA02016417.1', 'crispr_type': 'PossibleCRISPR', 'start': '14766', 'end': '14847', 'dr_consensus_seq': 'TACTCCTTCGTATTTAGCCCCGGTT', 'dr_consensus_length': '25', 'elements': [{'element_type': 'CRISPRdr', 'start': '14766', 'end': '14790'}, {'element_type': 'CRISPRspacer', 'start': '14791', 'end': '14822', 'spacer_seq': 'GGACGGCCATTACAGGTTGTCATAAGGAACTC', 'spacer_name': 'spacer_14791_32'}, {'element_type': 'CRISPRdr', 'start': '14823', 'end': '14847'}]}, {'crispr_id': 'AAAA02018790.1_16953_17048', 'genomic_accession': 'AAAA02018790.1', 'crispr_type': 'PossibleCRISPR', 'start': '16953', 'end': '17048', 'dr_consensus_seq': 'CAAGGCCACTAGATCTAATGTTG', 'dr_consensus_length': '23', 'elements': [{'element_type': 'CRISPRdr', 'start': '16953', 'end': '16975'}, {'element_type': 'CRISPRspacer', 'start': '16976', 'end': '17025', 'spacer_seq': 'GCGAGCTCGTGGCCATGGGATCCTGCATCGCGATGGCTGAAAGAGAGGAC', 'spacer_name': 'spacer_16976_50'}, {'element_type': 'CRISPRdr', 'start': '17026', 'end': '17048'}]}, {'crispr_id': 'AAAA02014824.1_9290_9397', 'genomic_accession': 'AAAA02014824.1', 'crispr_type': 'PossibleCRISPR', 'start': '9290', 'end': '9397', 'dr_consensus_seq': 'CGTCCTGCTGCTTCCGCCGCATCC', 'dr_consensus_length': '24', 'elements': [{'element_type': 'CRISPRdr', 'start': '9290', 'end': '9313'}, {'element_type': 'CRISPRspacer', 'start': '9314', 'end': '9373', 'spacer_seq': 'ACATCCTTCGCCGCCGCCGCCGCCGCCTGCGCCGATCCCCCCATCCATCCCGCGAGCGAG', 'spacer_name': 'spacer_9314_60'}, {'element_type': 'CRISPRdr', 'start': '9374', 'end': '9397'}]}, {'crispr_id': 'AAAA02007038.1_5522_5633', 'genomic_accession': 'AAAA02007038.1', 'crispr_type': 'PossibleCRISPR', 'start': '5522', 'end': '5633', 'dr_consensus_seq': 'CTGTAACAATTTATTATGTTATATCA', 'dr_consensus_length': '26', 'elements': [{'element_type': 'CRISPRdr', 'start': '5522', 'end': '5547'}, {'element_type': 'CRISPRspacer', 'start': '5548', 'end': '5607', 'spacer_seq': 'GTTTATTGCAAAAAGATCCTAATGTTCAGGTGAGAAAGATATAACCAATTCTGATGTGCT', 'spacer_name': 'spacer_5548_60'}, {'element_type': 'CRISPRdr', 'start': '5608', 'end': '5633'}]}, {'crispr_id': 'AAAA02008053.1_59432_59521', 'genomic_accession': 'AAAA02008053.1', 'crispr_type': 'PossibleCRISPR', 'start': '59432', 'end': '59521', 'dr_consensus_seq': 'CTGGAAACTGCAATGTGGTGCTG', 'dr_consensus_length': '23', 'elements': [{'element_type': 'CRISPRdr', 'start': '59432', 'end': '59454'}, {'element_type': 'CRISPRspacer', 'start': '59455', 'end': '59498', 'spacer_seq': 'TCAACTGTGGCTTCTGTCTCCATAATTACAGTGCTTGTTTCTCC', 'spacer_name': 'spacer_59455_44'}, {'element_type': 'CRISPRdr', 'start': '59499', 'end': '59521'}]}], 'num_crisprs': 14}


# if __name__ == "__main__":
#     wrapper()
gff_file, results_dict = parseGffToDict(results_path)
update_csv(results_dict, csv_path, summary_csv_path)
