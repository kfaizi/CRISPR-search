"""A Python (3.7) wrapper for identifying Cas protein orthologs.

Searches prokaryotic/metagenomic DNA sequence data for
new orthologs of known Cas proteins. Strings together
analysis and output from several existing tools.

Accepts a single FASTA file.
Generates an indexed nucleotide BLASTDB, searches it for ORFs matching a
protein query using TBLASTN, and then checks hits for the presence of DRs
with CRISPRFinder. The results are tabulated in a CSV file.
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
import logging
import texter
import signal

from pathlib import Path
from itertools import groupby
from Bio import SeqIO

# in working directory with CF + query, create blastdb/ and genomes/.
# genomes including cat fastas should be in genomes/.
# indexed dbs should be in blastdb/.

greeter = argparse.ArgumentParser(
    description='''For a given directory of FASTA files, this script generates an
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
hello = greeter.parse_args()
# python wcopy.py -in 'AAF' -wdir '~/search/' -bdir '~/search/blastdb/' -gdir '~/search/genomes/' -q 'Cas13d_proteins.fa' -s 'cf_v2.pl' -t '6'

######################## for hardcoded testing, use: ########################
# hello = greeter.parse_args(['-i', 'A',
#                             '-w', '~/search/',  # use pathlib to avoid / errors
#                             '-b', '~/search/blastdb/',  # use pathlib to avoid / errors
#                             '-g', '~/search/genomes/',
#                             '-q', 'Cas13d_proteins.fa',
#                             '-s', 'cf_v2.pl',
#                             '-t', '8'])

# catch sigints
signal.signal(signal.SIGINT, signal.default_int_handler)

# track date of last database update (NCBI)
NCBI_datestamp = "jul25"

# make all pathlike args full paths (resolve "~/")
working_path = Path(hello.working_dir).expanduser()
genomes_dir = Path(hello.genomes_dir).expanduser()
blastdb_dir = Path(hello.blastdb_dir).expanduser()

# Constants
name = hello.input_name + "_" + NCBI_datestamp
logname = "log_" + name + ".log"

blastdb_path = Path(blastdb_dir, hello.input_name, name)
genomes_path = Path(genomes_dir, hello.input_name)
query_path = sorted(working_path.glob(hello.query))[0]
og_script_path = sorted(working_path.glob(hello.script))[0]

gz_extractor_path = '/mnt/md0/kfaizi/cellophane/gz_extractor.sh'  # TEMP
zip_path = Path(working_path, 'raw', hello.input_name)

output_path = Path(working_path, "output", name)  # for each name
output_path.mkdir(parents=True, exist_ok=True)  # make dirs incl parent
script_path = Path(output_path, hello.script)

genome_path = Path(genomes_path, name).with_suffix(".fa")
cleaned_path = Path(genomes_path, name + "_cleaned").with_suffix(".fa")
removed_path = Path(genomes_path, name + "_removed").with_suffix(".fa")

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

# configure logging
logger = logging.getLogger()  # use root logger
logger.setLevel(logging.DEBUG)

sf = logging.Formatter(fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

file = logging.FileHandler(filename=logname, mode='a')
file.setLevel(logging.DEBUG)
file.setFormatter(sf)
logger.addHandler(file)

console = logging.StreamHandler(stream='ext://sys.stdout')
console.setLevel(logging.ERROR)
console.setFormatter(sf)
logger.addHandler(console)

# define functions


def move_script(path_to_script, path_to_copy):
    """To avoid CF errors, make a copy of it in each output directory."""
    subprocess.run(["cp",
                    f"{path_to_script}",
                    f"{path_to_copy}"],
                   cwd="/",
                   check=True)


def unzip_fasta(path_to_gz_extractor, path_to_zipped, path_to_dest):
    """(G)unzip compressed fasta files."""
    subprocess.run(["bash",
                    f"{path_to_gz_extractor}",
                    f"{path_to_zipped}",
                    f"{path_to_dest}"],
                   cwd="/",
                   check=True)


def cat_and_cut(path_to_fastas, path_to_cat):
    """Concatenate unzipped fastas, delete source."""
    cat_list = []

    for fa in sorted(path_to_fastas.glob("*.fsa_nt")):  # genomes_path
        cat_list.append(fa)

    combine = (['cat'] + cat_list)

    with open(path_to_cat, "w") as outfile:  # genome_path
        try:  # create master file...
            numfiles = str(len(cat_list))
            logger.info(f"Working. There are {numfiles} files to merge...")
            subprocess.run(combine,
                           cwd="/",
                           stdout=outfile,
                           check=True)
            logger.info("Done.")
        except subprocess.CalledProcessError as e:
            logger.critical(f"Error combining files: {e}")
            texter.send_text(f"Failed {name}")
            sys.exit()

    remove = (['rm'] + cat_list)

    try:  # ...then delete the individual files
        logger.info("Working. Removing the individual files...")
        subprocess.run(remove,
                       cwd='/',
                       check=True)
        logger.info("Done.")
    except subprocess.CalledProcessError as e:
        logger.critical(f"Error removing redundant individual files: {e}")
        texter.send_text(f"Failed {name}")
        sys.exit()


def dedupe_fasta(path_to_cat, path_to_new_cat, path_to_removed):
    """Dedupe fasta by sequence id."""
    cleaned_accs = {}
    removed_accs = {}
    skipcount = 0

    record = SeqIO.index(str(path_to_cat), 'fasta')
    logger.info(f"There are {len(record)} sequences to dedupe...")

    for idx in record:
        id_split = idx.split('|')

        if len(id_split) == 1:  # then it's an accession
            acc = idx
        elif 'gb' in id_split:
            accindex = id_split.index('gb') + 1
            acc = id_split[accindex]
        elif 'dbj' in id_split:
            accindex = id_split.index('dbj') + 1
            acc = id_split[accindex]
        else:
            logger.warning(f"Warning! Sequence id format unrecognized: {idx}")

        try:
            if acc not in cleaned_accs:
                cleaned_accs[acc] = idx
            else:  # then it's a duplicate; keep track
                removed_accs[acc] = idx
        except Exception:
            skipcount += 1
            pass

    with open(path_to_new_cat, 'w+') as cleanfile:
        for ids in cleaned_accs.values():
            SeqIO.write(record[ids], cleanfile, 'fasta')  # new deduped fasta

        if skipcount >= 1:
            logger.warning(f"Warning! Skipped {skipcount} sequences.")

    with open(path_to_removed, 'w+') as remfile:
        for ids in removed_accs.values():
            SeqIO.write(record[ids], remfile, 'fasta')  # holds deleted seqs

    return len(cleaned_accs), len(removed_accs)

    # try:  # deletes og fasta to make space. will validate b4 using this
    #     logger.info(f"Deleting original file {path_to_cat}")
    #     subprocess.run(["rm",
    #                    f"{path_to_cat}"],
    #                    cwd='/',
    #                    check=True)
    # except subprocess.CalledProcessError as e:
    #     logger.critical(f"Error deleting original file: {e}")
    #     texter.send_text(f"Failed {name}")
    #     sys.exit()


def make_db(path_to_fasta, path_to_new_database):
    """Concatenated fasta > parsed blastdb."""
    p = subprocess.run(["makeblastdb",
                        "-in",
                        f"{path_to_fasta}",  # in ~/genomes/
                        "-parse_seqids",
                        "-max_file_sz",
                        "4GB",
                        "-dbtype",
                        "nucl",
                        "-out",
                        f"{path_to_new_database}"],
                       cwd="/",
                       text=True,
                       capture_output=True,
                       check=True)
    results = str(p.stdout)
    logger.info(results)


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
        try:
            linecount = sum(1 for line in ID_finder)  # this iterates through the csv...
            logger.info(f"Working. There are {linecount} contig IDs to parse...")

            if linecount == 0:  # then blast found no results
                logger.info("Blast found no hits; finishing up.")
                sys.exit()

            else:
                f.seek(0)  # ...so we need to move back to start
                ID_finder = csv.reader(f)  # reset generator
                for line in ID_finder:
                    output_list.append(line[0])  # make list of IDs
                return output_list

        except csv.Error as e:
            logger.critical(f"Error gathering accession data for {path_to_csv}, line {ID_finder.line_num}: {e}")
            texter.send_text(f"Failed {name}")
            sys.exit()

        finally:
            if len(output_list) != linecount:
                logger.critical("Error: initial and final line counts do not match! Exiting")
                texter.send_text(f"Failed {name}")
                sys.exit()


def list_csv(path_to_csv, path_to_new_IDs):
    """Write accession IDs to .txt, creating one if it doesn't exist."""
    with open(path_to_new_IDs, 'w+', newline='') as f:
        try:
            output_list = read_csv(path_to_csv)
            logger.info(f"Writing IDs to text file {path_to_new_IDs}...")

            for ID in output_list:
                f.write(f"{ID}\n")
            logger.info("Done!")

        except csv.Error as e:
            logger.critical(f"Error parsing csv for {path_to_csv}: {e}")
            texter.send_text(f"Failed {name}")
            sys.exit()


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
    """Dedupe fasta by sequence; from https://www.biostars.org/p/3003/."""
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
    results = str(p.stdout)
    logger.info(results)


def parseGffToDict(path_to_gff):  # ie results_path
    """Take gff filename, return parsed dict. Adapted from PDH."""
    # dict will be populated with crispr array info from gff
    gff_file = sorted(path_to_gff.glob("*.gff"))[0]  # find actual gff file
    results_dict = {'crisprs_list': []}
    with open(gff_file) as file:
        data = file.read()  # entire gff file

    data_split = data.split('\n\n')  # split into found crisprs
    data_split = list(filter(None, data_split))  # remove empty elements

    if data_split == ['##gff-version 3\n']:
        logger.info("CRISPRFinder found no hits.")
        results_dict['crisprs_list'] = None
        results_dict['num_crisprs'] = 0

    else:
        num_crisprs = len(data_split)
        results_dict['num_crisprs'] = num_crisprs

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
    df_blast['acc2'] = df_blast['saccver']  # make duplicate accession column for later

    # send parsed crisprfinder results into dataframe
    if results_dict['num_crisprs'] == 0:
        # do stuff
        pass

    else:
        d_crispr = results_dict['crisprs_list']
        df_crispr = pd.DataFrame(d_crispr,
                                 dtype=object,
                                 columns=['genomic_accession',
                                          'crispr_id',
                                          'crispr_type',
                                          'start',
                                          'end',
                                          'dr_consensus_seq',
                                          'dr_consensus_length',
                                          'elements'],
                                 index=None)

        # make master dataframe. creates duplicates as needed to
        # handle multiple hits from blast/crisprfinder/both
        combo = pd.merge(df_blast, df_crispr, how='left', left_on='saccver', right_on='genomic_accession')
        combo = combo.drop(columns='genomic_accession')  # redundant; use crispr_id

        # make new summary dataframe
        summary = combo[['sseq',
                         'acc2',
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
                         'length',
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
                           'acc2',
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
                           'length',
                           'start',
                           'end',
                           'elements']]

        # rename columns more intuitively
        summary = summary.rename(columns={'sseq': 'protein_sequence',  # only aligned part
                                          'acc2': 'accession_number',
                                          'slen': 'contig_length',
                                          'qseqid': 'ortholog_match',
                                          'dr_consensus_seq': 'DR_consensus_seq',
                                          'dr_consensus_length': 'DR_consensus_length',
                                          'sstart': 'protein_start',
                                          'send': 'protein_end',
                                          'length': 'protein_length',
                                          'start': 'array_start',
                                          'end': 'array_end',
                                          'elements': 'array'})

        # parse contigs for subsequent extraction, without overloading memory
        record_dict = SeqIO.index(str(dedupe_path), 'fasta')  # can't handle pure paths

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

                    # # print('For', acc, 'start at', extracted_start, 'and end at', extracted_end, '. The seq is', length, 'long, of which we\'re taking', extracted_len)
                    extracted_head = str(record_dict[acc].description)
                    extracted_name = crispr_id + "_extracted.fa"
                    extracted_path = Path(results_path, extracted_name)

                    # write extracted sequence to a new fasta file...
                    with open(extracted_path, 'w') as outfile:
                        SeqIO.write(extracted_seq, outfile, 'fasta')

                    logger.info(f"Wrote extracted CRISPR locus to {extracted_path}")

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
    try:  # initialize
        texter.send_text(f"Starting search on {name}")
        logger.info("Starting...")
    except Exception as e:
        logger.error(f"Something went wrong while starting: {e}")
        texter.send_text(f"Startup {name} failed :(")
        sys.exit()

    try:  # copy CF to each result dir to avoid error
        move_script(og_script_path, script_path)
    except subprocess.CalledProcessError as e:
        logger.critical(f"Error initializing file structure: {e}")
        texter.send_text(f"Failed {name}")
        sys.exit()

    try:  # unzip fastas
        logger.info(f"Extracting zipped fastas {hello.input_name}...")
        unzip_fasta(gz_extractor_path, zip_path, genomes_path)
        logger.info(f"Success! Unzipped fastas to {genomes_path}")
    except subprocess.CalledProcessError as e:
        logger.critical(f"Error extracting zipped fastas: {e}")
        texter.send_text(f"Failed {name}")
        sys.exit()

    try:  # combine fastas and make space
        logger.info(f"Concatenating fastas at {genomes_path} and cleaning up...")
        cat_and_cut(genomes_path, genome_path)
        logger.info(f"Success! New file created at {genome_path}")
    except subprocess.CalledProcessError as e:
        logger.critical(f"Error combining and/or deleting fastas: {e}")
        texter.send_text(f"Failed {name}")
        sys.exit()

    try:  # deduplicating the raw data
        logger.info(f"Deduplicating sequences in {genome_path}...")
        num_saved, num_removed = dedupe_fasta(genome_path, cleaned_path, removed_path)
        logger.info(f"Success! Cleaned fasta created at {cleaned_path}. There were {num_saved} unique sequences with {num_removed} duplicates, which are saved at {removed_path}")
        texter.send_text(f"deduped {name}")
    except Exception as e:
        logger.critical(f"Error deduping: {e}")
        texter.send_text(f"Failed {name}")
        sys.exit()

    try:  # making the blastdb
        logger.info(f"Making blastdb {name}...")
        make_db(cleaned_path, blastdb_path)
        logger.info(f"Success! Blast database created at {blastdb_path}")
        texter.send_text(f"made db {name}")
    except subprocess.CalledProcessError as e:
        logger.critical(f"Error making blastdb: {e}")
        texter.send_text(f"Failed {name}")
        sys.exit()

    try:  # blasting the blastdb
        logger.info(f"Blasting proteins in {query_path} against {name} with {hello.threads} threads...")
        search_db(blastdb_path, query_path, archive_path, hello.threads)
        logger.info(f"Success! Blast results written to {archive_path}")
        texter.send_text(f"blasted {name}")
    except subprocess.CalledProcessError as e:
        logger.critical(f"Error blasting {query_path} against blastdb {name}: {e}")
        texter.send_text(f"Failed {name}")
        sys.exit()

    try:  # creating a .csv from the results
        logger.info("Creating .csv file from the blast results...")
        make_csv(archive_path, csv_path)
        logger.info(f"Success! New .csv file created at {csv_path}")
    except subprocess.CalledProcessError as e:
        logger.critical(f"Error creating .csv file from the archive at {archive_path}: {e}")
        texter.send_text(f"Failed {name}")
        sys.exit()

    try:  # make a list of accession IDs
        logger.info("Extracting accession ID numbers from the .csv file...")
        list_csv(csv_path, accessions_path)
        logger.info(f"Success! Accession list created at {accessions_path}")
    except Exception as e:
        logger.critical(f"Error creating accession list at {accessions_path}: {e}")
        texter.send_text(f"Failed {name}")
        sys.exit()

    try:  # extract the corresponding contigs
        logger.info(f"Reading accession list at {accessions_path} and extracting corresponding contigs from {blastdb_path}...")
        extract_contigs(blastdb_path, accessions_path, contigs_path)
        logger.info(f"Success! Contigs extracted to {contigs_path}")
    except Exception as e:
        logger.critical(f"Error extracting contigs to {contigs_path}: {e}")
        texter.send_text(f"Failed {name}")
        sys.exit()

    try:  # deduplicate the contigs fasta
        logger.info(f"Deduplicating BLAST hits at {contigs_path}...")
        dedupe_contigs(contigs_path, dedupe_path)
        logger.info(f"Success! Contigs deduplicated and moved to {dedupe_path}")
    except Exception as e:
        logger.critical(f"Error deduplicating contigs at {contigs_path}: {e}")
        texter.send_text(f"Failed {name}")
        sys.exit()

    try:  # search 'em with crisprfinder
        logger.info(f"Searching sequences in {dedupe_path} for DRs...")
        find_CRISPRs(script_path, dedupe_path, results_path)
        logger.info(f"Success! CRISPRFinder output created at {results_path}")
    except subprocess.CalledProcessError as e:
        logger.critical(f"Error searching for DRs: {e}")
        texter.send_text(f"Failed {name}")
        sys.exit()

    try:  # parse the gff output
        logger.info("Parsing the CRISPRFinder output...")
        gff_file, results_dict = parseGffToDict(results_path)
        logger.info(f"Success! Parsed gff {gff_file} to dict")
    except Exception as e:
        logger.critical(f"Error parsing results at {results_path}: {e}")
        texter.send_text(f"Failed {name}")
        sys.exit()

    try:  # update the csv
        logger.info(f"Updating the .csv at {csv_path} with new info...")
        update_csv(results_dict, csv_path, summary_csv_path)
        logger.info(f"Success! Final output file located at {summary_csv_path}")
        texter.send_text(f"All done {name}")
    except Exception as e:
        logger.critical(f"Error updating .csv at {csv_path}: {e}")
        texter.send_text(f"Failed {name}")
        sys.exit()


if __name__ == "__main__":
    try:
        wrapper()
    except KeyboardInterrupt:
        logger.error("Keyboard interrupt")
        sys.exit()




