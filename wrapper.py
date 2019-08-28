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
# python wcopy.py -in 'C' -wdir '~/search/' -bdir '~/search/blastdb/' -gdir '~/search/genomes/' -q 'Cas13d_proteins.fa' -s 'cf_v2.pl' -t '6'

######################## for hardcoded testing, use: ########################
# hello = greeter.parse_args(['-i', 'C',
#                             '-w', '~/search/',  # use pathlib to avoid / errors
#                             '-b', '~/search/blastdb/',  # use pathlib to avoid / errors
#                             '-g', '~/search/genomes/',
#                             '-q', 'Cas13d_proteins.fa',
#                             '-s', 'cf_v2.pl',
#                             '-t', '4'])

# hello = greeter.parse_args(['-i', 'AAA',
#                             '-w', '~/meethere/search/',  # use pathlib to avoid / errors
#                             '-b', '~/meethere/search/blastdb/',  # use pathlib to avoid / errors
#                             '-g', '~/meethere/search/genomes/',
#                             '-q', 'Cas13d_proteins.fa',
#                             '-s', 'cf_v2.pl',
#                             '-t', '2'])

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
int_cleaned_path = Path(genomes_path, name + "_int_cleaned").with_suffix(".fa")
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


def move_script(script, copy):
    """To avoid CF errors, make a copy of it in each output directory."""
    subprocess.run(["cp",
                    f"{script}",
                    f"{copy}"],
                   cwd="/",
                   check=True)


def unzip_fasta(gz_extractor, zipped, unzipped):
    """(G)unzip compressed fasta files."""
    subprocess.run(["bash",
                    f"{gz_extractor}",
                    f"{zipped}",
                    f"{unzipped}"],
                   cwd="/",
                   check=True)


def cat_and_cut(lone, grouped):
    """Concatenate unzipped fastas, delete source."""
    cat_list = []

    for fa in sorted(lone.glob("*.fsa_nt")):  # genomes_path
        cat_list.append(fa)

    combine = (['cat'] + cat_list)

    with open(grouped, "w") as outfile:  # genome_path
        try:  # create master file...
            numfiles = str(len(cat_list))
            logger.info(f"Working. There are {numfiles} files to merge...")
            subprocess.run(combine,
                           cwd="/",
                           stdout=outfile,
                           check=True)
            logger.info("Done.")
        except subprocess.CalledProcessError as e:
            logger.critical(f"Error combining files: {e}", exc_info=True)
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
        logger.critical(f"Error removing individual files: {e}", exc_info=True)
        texter.send_text(f"Failed {name}")
        sys.exit()


def dedupe_raw(cat, int_cleaned, cleaned, removed):
    """Dedupe fasta by sequence id."""
    cleaned_accs = {}
    removed_accs = {}
    skipcount = 0

    # looks for exact header duplicates and removes them. Needed for NCBI-Q
    # adapted from https://www.biostars.org/p/3003/
    seqids = set()
    with open(cat, 'r') as infile:
        with open(int_cleaned, 'w+') as outfile:
            head = None
            for h, lines in groupby(infile, lambda x: x.startswith('>')):
                if h is True:
                    head = next(lines)
                else:
                    seq = ''.join(lines)
                if head not in seqids:
                    seqids.add(head)
                    outfile.write(f"{head}{seq}")
                elif head in seqids:
                    logger.info(f"Exact duplicate removed: {head}")

    # now this should deal with cases of same accession, different header
    record = SeqIO.index(str(int_cleaned), 'fasta')
    logger.info(f"There are {len(record)} sequences to dedupe...")

    for idx in record:
        id_split = idx.split('|')

        if len(id_split) == 1:  # then it's an accession
            acc = idx
        elif len(id_split) == 3:  # form "dbj|AAAAAnnn|"
            acc = id_split[1]
        elif len(id_split) == 5:  # form "gb|123|dbj|AAAAAnnn|"
            acc = id_split[3]
        else:
            logger.warning(f"Warning! Sequence id format unrecognized: {idx}")
            acc = None

        if acc is None:
            skipcount += 1
        elif acc not in cleaned_accs:
            cleaned_accs[acc] = idx
        elif acc in cleaned_accs:  # then it's a duplicate; keep track
            removed_accs[acc] = idx

    with open(cleaned, 'w+') as cleanfile:
        for ids in cleaned_accs.values():
            SeqIO.write(record[ids], cleanfile, 'fasta')  # new deduped fasta

        if skipcount >= 1:
            logger.warning(f"Warning! Skipped {skipcount} sequences.")

    with open(removed, 'w+') as remfile:
        for ids in removed_accs.values():
            SeqIO.write(record[ids], remfile, 'fasta')  # holds deleted seqs

    return len(cleaned_accs), len(removed_accs)


def make_db(fasta, db):
    """Concatenated fasta > parsed blastdb."""
    p = subprocess.run(["makeblastdb",
                        "-in",
                        f"{fasta}",  # in ~/genomes/
                        "-parse_seqids",
                        "-max_file_sz",
                        "4GB",
                        "-dbtype",
                        "nucl",
                        "-out",
                        f"{db}"],
                       cwd="/",
                       text=True,
                       capture_output=True,
                       check=True)
    results = str(p.stdout)
    logger.info(results)


def search_db(db, query, archive, threads):
    """Blastdb + query > hits in ASN.1 format."""
    subprocess.run(["tblastn",
                    "-db",
                    f"{db}",
                    "-query",
                    f"{query}",  # known protein reference
                    "-out",
                    f"{archive}",
                    "-evalue",
                    "10",  # could go bigger...
                    "-num_threads",
                    f"{threads}",
                    "-outfmt",
                    "11"],  # ASN.1
                   cwd="/",
                   check=True)


def make_csv(archive, csvfile):
    """ASN.1 > csv with headers."""
    subprocess.run(["blast_formatter",
                    "-archive",
                    f"{archive}",
                    "-outfmt",
                    sections,
                    "-out",
                    f"{csvfile}"],
                   cwd="/",
                   check=True)


def read_csv(csvfile):
    """Get accession IDs from csv."""
    with open(csvfile, newline='') as f:
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
                    if line[0] not in output_list:  # deduplicates blast hits by accession
                        output_list.append(line[0])  # make list of IDs
                return output_list

        except csv.Error as e:
            logger.critical(f"Error gathering accession data for {csvfile}, line {ID_finder.line_num}: {e}", exc_info=True)
            texter.send_text(f"Failed {name}")
            sys.exit()


def list_csv(csvfile, accessions):
    """Write accession IDs to .txt, creating one if it doesn't exist."""
    with open(accessions, 'w+', newline='') as f:
        try:
            output_list = read_csv(csvfile)
            logger.info(f"Writing IDs to text file {accessions}...")

            for ID in output_list:
                f.write(f"{ID}\n")
            logger.info("Done!")

        except csv.Error as e:
            logger.critical(f"Error parsing csv for {csvfile}: {e}", exc_info=True)
            texter.send_text(f"Failed {name}")
            sys.exit()


def extract_contigs(db, accessions, contigs):
    """Extract contigs by referencing sseqids to blastdb."""
    subprocess.run(["blastdbcmd",
                    "-db",
                    f"{db}",  # parent
                    "-dbtype",
                    "nucl",
                    "-entry_batch",
                    f"{accessions}",  # list of sseqids to extract
                    "-out",
                    f"{contigs}"],
                   cwd="/",
                   check=True)


def find_CRISPRs(cf, contigs, gff):
    """Search extracted contigs for DRs; retain output."""
    p = subprocess.run(["perl",
                        f"{cf}",
                        f"{contigs}",
                        f"{gff}"],
                       cwd=output_path,  # prevents CRISPRFinder from erring
                       capture_output=True,
                       text=True,
                       check=True)
    results = str(p.stdout)
    logger.info(results)


def parseGffToDict(gff):  # ie results_path
    """Take gff filename, return parsed dict. Adapted from PDH."""
    # dict will be populated with crispr array info from gff
    gff_file = sorted(gff.glob("*.gff"))[0]  # find actual gff file
    results_dict = {'crisprs_list': []}
    with open(gff_file) as file:
        data = file.read()  # entire gff file

    data_split = data.split('\n\n')  # split into found crisprs
    data_split = list(filter(None, data_split))  # remove empty elements

    if data_split == ['##gff-version 3\n']:  # then no results
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


def update_csv(results_dict, blast_csv, summary_csv):
    """Add headers, and append CF data."""
    # make blast csv into dataframe
    df_blast = pd.read_csv(blast_csv, dtype=object, header=None, names=sections_list, index_col=None)
    df_blast['acc2'] = df_blast['saccver']  # preserve extra accession column for after merge

    # send parsed crisprfinder results into dataframe
    if results_dict['num_crisprs'] == 0:
        new_cols = ['sseq',
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
                    'sframe',
                    'length',
                    'start',
                    'end',
                    'elements']

        summary = df_blast.reindex(columns=new_cols, copy=True)

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
                         'sframe',
                         'length',
                         'start',
                         'end',
                         'elements']].copy()

    # create columns w/ whether DR found (y/n) and if yes, how many
    dr_status = []
    dr_sum = []
    array_lens = []

    for eltup in summary[['dr_consensus_seq', 'elements', 'start', 'end']].itertuples():
        dr_count = 0
        dr_exists = 'no'
        array_len = 0

        # if dr_consensus_seq is there, then DR was found
        if isinstance(eltup.dr_consensus_seq, str):
            dr_exists = 'yes'

            array_len = abs(int(eltup.end) - int(eltup.start)) + 1

            for dictx in eltup.elements:
                for v in dictx.values():
                    if v == "CRISPRdr":
                        dr_count += 1

        dr_status.append(dr_exists)
        dr_sum.append(dr_count)
        array_lens.append(array_len)

    # add columns to summary
    summary.loc[:, 'DR_found'] = dr_status
    summary.loc[:, 'num_DRs'] = dr_sum
    summary.loc[:, 'array_length'] = array_lens

    # rename columns more intuitively
    summary = summary.rename(columns={'sseq': 'protein_sequence',  # only aligned part
                                      'acc2': 'accession_number',
                                      'qseqid': 'ortholog_match',
                                      'dr_consensus_seq': 'DR_consensus_seq',
                                      'dr_consensus_length': 'DR_consensus_length',
                                      'slen': 'contig_length',
                                      'sstart': 'protein_start',
                                      'send': 'protein_end',
                                      'length': 'protein_length',  # in amino acids
                                      'sframe': 'protein_frame',
                                      'start': 'array_start',
                                      'end': 'array_end',
                                      'elements': 'array'})

    # parse contigs for subsequent extraction, without overloading memory
    record_dict = SeqIO.index(str(contigs_path), 'fasta')  # can't handle pure paths
    # extract up to +/- 20kb flanking sequence around crispr locus
    extracted_heads = []
    extracted_paths = []
    distances = []

    for row in summary.itertuples():
        extracted_head = None
        extracted_path = None
        dist = None  # protein-array distance (nt)
        acc = row.accession_number

        if row.DR_found == "yes":  # then extract flanking sequence

            crispr_id = row.crispr_id

            extracted_start = max(int(row.array_start)-20000, 1)
            extracted_end = min(int(row.array_end)+20000, int(row.contig_length))
            extracted_seq = (record_dict[acc])[extracted_start-1: extracted_end]

            extracted_head = str(record_dict[acc].description)
            extracted_name = crispr_id + "_extracted.fa"
            extracted_path = Path(results_path, extracted_name)

            # write extracted sequence to a new fasta file...
            with open(extracted_path, 'w') as outfile:
                SeqIO.write(extracted_seq, outfile, 'fasta')
            logger.info(f"Wrote extracted CRISPR locus to {extracted_path}")

            max_prot = max(int(row.protein_start), int(row.protein_end))
            max_array = max(int(row.array_start), int(row.array_end))

            # clearly, assuming no overlap
            # if the array is downstream of the protein,
            if max_array > max_prot:
                dist = max_array - max_prot - int(row.array_length)
            # else if the protein is downstream of the array,
            elif max_prot > max_array:
                dist = max_prot - max_array - int(row.protein_length)

        else:  # then just extract header
            extracted_head = str(record_dict[acc].description)

        extracted_heads.append(extracted_head)
        extracted_paths.append(extracted_path)
        distances.append(dist)

    # add header/filename info to the summary dataframe
    summary.loc[:, 'extracted_header'] = extracted_heads
    summary.loc[:, 'extracted_filename'] = extracted_paths
    summary.loc[:, 'distance_protein_and_array'] = distances

    # reorder columns for readability
    summary = summary[['protein_sequence',
                       'accession_number',
                       'extracted_header',
                       'ortholog_match',
                       'ppos',
                       'pident',
                       'evalue',
                       'bitscore',
                       'crispr_type',
                       'crispr_id',
                       'DR_found',
                       'num_DRs',
                       'DR_consensus_seq',
                       'DR_consensus_length',
                       'contig_length',
                       'protein_start',
                       'protein_end',
                       'protein_length',
                       'protein_frame',
                       'array_start',
                       'array_end',
                       'array_length',
                       'distance_protein_and_array',
                       'array',
                       'extracted_filename']]

    # drop the redundant accessions (they're in the headers)
    summary = summary.drop(columns='accession_number')

    # write summary dataframe to csv
    with open(summary_csv, 'w') as f:
        summary.to_csv(f, mode='w', index=False, header=True)


def wrapper():
    """Put it all together."""
    try:  # initialize
        texter.send_text(f"Starting search on {name}")
        logger.info("Starting...")
    except Exception as e:
        logger.error(f"Something went wrong while starting: {e}", exc_info=True)
        texter.send_text(f"Startup {name} failed :(")
        sys.exit()

    try:  # copy CF to each result dir to avoid error
        move_script(og_script_path, script_path)
    except subprocess.CalledProcessError as e:
        logger.critical(f"Error initializing file structure: {e}", exc_info=True)
        texter.send_text(f"Failed copy {name}")
        sys.exit()

    try:  # unzip fastas
        logger.info(f"Extracting zipped fastas {hello.input_name}...")
        unzip_fasta(gz_extractor_path, zip_path, genomes_path)
        logger.info(f"Success! Unzipped fastas to {genomes_path}")
    except subprocess.CalledProcessError as e:
        logger.critical(f"Error extracting zipped fastas: {e}", exc_info=True)
        texter.send_text(f"Failed extract {name}")
        sys.exit()

    try:  # combine fastas and make space
        logger.info(f"Concatenating fastas at {genomes_path} and cleaning up...")
        cat_and_cut(genomes_path, genome_path)
        logger.info(f"Success! New file created at {genome_path}")
    except subprocess.CalledProcessError as e:
        logger.critical(f"Error combining and/or deleting fastas: {e}", exc_info=True)
        texter.send_text(f"Failed concat {name}")
        sys.exit()

    try:  # deduplicating the raw data
        logger.info(f"Deduplicating sequences in {genome_path}...")
        num_saved, num_removed = dedupe_raw(genome_path, int_cleaned_path, cleaned_path, removed_path)
        logger.info(f"Success! Cleaned fasta created at {cleaned_path}. There were {num_saved} unique sequences with {num_removed} duplicates, moved to {removed_path}")
    except Exception as e:
        logger.critical(f"Error deduping: {e}", exc_info=True)
        texter.send_text(f"Failed dedupe {name}")
        sys.exit()

    try:  # making the blastdb
        logger.info(f"Making blastdb {name}...")
        make_db(cleaned_path, blastdb_path)
        logger.info(f"Success! Blast database created at {blastdb_path}")
    except subprocess.CalledProcessError as e:
        logger.critical(f"Error making blastdb: {e}", exc_info=True)
        texter.send_text(f"Failed blastdb {name}")
        sys.exit()

    try:  # blasting the blastdb
        logger.info(f"Blasting proteins in {query_path} against {name} with {hello.threads} threads...")
        search_db(blastdb_path, query_path, archive_path, hello.threads)
        logger.info(f"Success! Blast results written to {archive_path}")
    except subprocess.CalledProcessError as e:
        logger.critical(f"Error blasting {query_path} against blastdb {name}: {e}", exc_info=True)
        texter.send_text(f"Failed blast {name}")
        sys.exit()

    try:  # creating a .csv from the results
        logger.info("Creating .csv file from the blast results...")
        make_csv(archive_path, csv_path)
        logger.info(f"Success! New .csv file created at {csv_path}")
    except subprocess.CalledProcessError as e:
        logger.critical(f"Error creating .csv file from the archive at {archive_path}: {e}", exc_info=True)
        texter.send_text(f"Failed csv {name}")
        sys.exit()

    try:  # make a list of accession IDs
        logger.info("Extracting accession ID numbers from the .csv file...")
        list_csv(csv_path, accessions_path)
        logger.info(f"Success! Accession list created at {accessions_path}")
    except Exception as e:
        logger.critical(f"Error creating accession list at {accessions_path}: {e}", exc_info=True)
        texter.send_text(f"Failed acc {name}")
        sys.exit()

    try:  # extract the corresponding contigs
        logger.info(f"Reading accession list at {accessions_path} and extracting corresponding contigs from {blastdb_path}...")
        extract_contigs(blastdb_path, accessions_path, contigs_path)
        logger.info(f"Success! Contigs extracted to {contigs_path}")
    except Exception as e:
        logger.critical(f"Error extracting contigs to {contigs_path}: {e}", exc_info=True)
        texter.send_text(f"Failed acc 2 {name}")
        sys.exit()

    try:  # search 'em with crisprfinder
        logger.info(f"Searching sequences in {contigs_path} for DRs...")
        find_CRISPRs(script_path, contigs_path, results_path)
        logger.info(f"Success! CRISPRFinder output created at {results_path}")
    except subprocess.CalledProcessError as e:
        logger.critical(f"Error searching for DRs: {e}", exc_info=True)
        texter.send_text(f"Failed cf {name}")
        sys.exit()

    try:  # parse the gff output
        logger.info("Parsing the CRISPRFinder output...")
        gff_file, results_dict = parseGffToDict(results_path)
        logger.info(f"Success! Parsed gff {gff_file} to dict")
    except Exception as e:
        logger.critical(f"Error parsing results at {results_path}: {e}", exc_info=True)
        texter.send_text(f"Failed parse {name}")
        sys.exit()

    try:  # update the csv
        logger.info(f"Updating the .csv at {csv_path} with new info...")
        update_csv(results_dict, csv_path, summary_csv_path)
        logger.info(f"Success! Final output file located at {summary_csv_path}")
        texter.send_text(f"Done {name}")
    except Exception as e:
        logger.critical(f"Error updating .csv at {csv_path}: {e}", exc_info=True)
        texter.send_text(f"Failed update {name}")
        sys.exit()


if __name__ == "__main__":
    try:
        wrapper()
    except KeyboardInterrupt:
        logger.error("Keyboard interrupt")
        sys.exit()
