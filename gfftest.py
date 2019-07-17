"""Docstring."""
# pdh
# june 17, 11pm. v1
# this python script parses gff files and creates 4 files:
# [1] <assembly_accession>__<genomic_accession>_<crispr_start>_<crispr_end>__extracted__feature_table_new.txt
# [2] <assembly_accession>__<genomic_accession>_<crispr_start>_<crispr_end>__extracted.fasta
# [3] <assembly_accession>__<genomic_accession>_<crispr_start>_<crispr_end>__extracted__<n>_proteins.fasta
# [4] <assembly_accession>__<genomic_accession>_<crispr_start>_<crispr_end>__extracted__<n>_proteins_nt.fasta

import glob

from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
from Bio import Entrez
Entrez.email = 'pdhsu@salk.edu'  # required for esearch/efetch anti-abuse

# CONSTANTS
gff_dir = Path("/Users/kianfaizi/devtemp/test")  # dir of gff output from crisprfinder
genome_dir = Path("/Users/kianfaizi/devtemp/test")  # dir of *_genomic.fna and *_feature_table.txt
extracted_dir = Path("/Users/kianfaizi/devtemp/test")


def parseGffToDict(filename):
    """Take gff filename, return parsed dict."""
    # dict will be populated with crispr array info from gff
    dict = {'crisprs_list': []}
    with open(Path(gff_dir, filename)) as file:
        data = file.read()  # entire gff file
    data_split = data.split('\n\n')  # split into found crisprs
    data_split = list(filter(None, data_split))  # remove empty elements
    num_crisprs = len(list(data_split))
    dict.update({
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
        dict['crisprs_list'].append(dict_item)

    print(f'Parsed gff {filename} to dict')


# # extract sequence of 40kb locus around found crispr to file
# def extractSequence(crispr, record_dict, assembly_accession):

#   # grab SeqRecord for genomic accession that this found crispr is in
#   record = record_dict[crispr['genomic_accession']]

#   crispr_id = crispr['crispr_id']
#   start = int(crispr['start'])
#   end = int(crispr['end'])
#   seq = record.seq[start-1: end]
#   record_length = len(record.seq)

#   # extract up to 20kb upstream and downstream from crispr
#   extracted_start = max(start-20000, 1)  # go upstream 20kb or start at 1
#   extracted_end = min(end+20000, record_length)  # go downstream 20kb or stop at record_length
#   extracted_length = extracted_end - extracted_start + 1
#   extracted_seq = record.seq[extracted_start-1: extracted_end]
#   print(f'Extracted {extracted_start} to {extracted_end}')

#   # write extracted fasta to file
#   extracted_id = '{}__{}__extracted'.format(assembly_accession, crispr_id)
#   extracted_record = SeqRecord(
#     extracted_seq,
#     id = extracted_id,
#     name = extracted_id,
#     description = """EXTRACTED LOCUS, assembly accession: {}, genomic accession: {}, crispr id: {}, extracted_start: {}, extracted_end: {}, extracted_length: {}""".format(
#         assembly_accession,
#         crispr['genomic_accession'],
#         crispr_id,
#         extracted_start,
#         extracted_end,
#         extracted_length
#       )
#     )

#   # write to file
#   extracted_path = '{}{}.fasta'.format(extracted_dir, extracted_id)
#   with open(extracted_path, 'w') as outf:
#     SeqIO.write(extracted_record, outf, 'fasta')

#   print(f'Wrote extracted crispr locus {crispr_id} to {extracted_path}')
#   return extracted_id, extracted_start, extracted_end


# take crispr dict, write new feature table
def createNewFeatureTable(crispr, record_dict, extracted_id, extracted_start, extracted_end, assembly_accession):

  record = record_dict[crispr['genomic_accession']]

  crispr_id = crispr['crispr_id']
  start = int(crispr['start'])
  end = int(crispr['end'])
  seq = record.seq[start-1: end]
  length = end - start + 1

  # open original feature table

  ft_filename = assembly_accession + '_feature_table.txt'
  try:
    with open(genome_dir + ft_filename) as file:
      ft = file.read()
  except EnvironmentError as err:
    print(f'ERROR opening feature table: {err}')
    return None, None

  # with open(genome_dir + ft_filename) as file:
  #   ft = file.read()

  ft_split = ft.split('\n')  # split feature table into rows
  ft_split = filter(None, ft_split)  # remove empty elements

  newft = []  # will hold new feature table
  protein_seqrecords = []
  protein_nt_seqrecords = []

  header_orig = ft_split[0]  # save original header
  ft_split = ft_split[1:]  # remove header from original ft, easier to process

  for line in ft_split:  # for row in feature table
    line_split = line.split('\t')

    feature_ga = line_split[6]
    feature_start = int(line_split[7])
    feature_end = int(line_split[8])

    if feature_ga == crispr['genomic_accession'] and \
       feature_start >= extracted_start and \
       feature_end <= extracted_end:  # keep rows in range and add to newft

      line_split.insert(9, '')  # add sequence column after start/end
      feature_refseq = line_split[12]  # non-redundant_refseq column, e.g. WP_001205278.1
      if len(feature_refseq) > 0:
        fetched_header, fetched_seq = fetchProteinSeq(feature_refseq)
        line_split[9] = fetched_seq  # write seq into seq column
        line_split[14] += '; ' + fetched_header  # write header into name column

        # create protein seqrecord to return
        protein_id = '{}__{}__{}'.format(assembly_accession, crispr_id, feature_refseq)
        protein_seqrecord = SeqRecord(
          Seq(fetched_seq, SingleLetterAlphabet()),
          id=protein_id,
          name=protein_id,
          description="""efetch header: {}, non-redundant_refseq: {}, assembly: {}, genomic_accession: {}, nt start: {}, nt end: {}, nt length: {}, product aa length: {}""".format(
              fetched_header[1:],
              feature_refseq,
              assembly_accession,
              feature_ga,
              feature_start,
              feature_end,
              line_split[18],
              line_split[19]
            )
        )
        protein_seqrecords.append(protein_seqrecord)

        # create protein nt seqrecord to return
        protein_nt_id = '{}__{}__{}'.format(assembly_accession, crispr_id, feature_refseq)
        protein_nt_seqrecord = SeqRecord(
          record.seq[feature_start-1: feature_end],
          id=protein_nt_id,
          name=protein_nt_id,
          description="""**nt sequence** efetch header: {}, non-redundant_refseq: {}, assembly: {}, genomic_accession: {}, nt start: {}, nt end: {}, nt length: {}""".format(
              fetched_header[1:],
              feature_refseq,
              assembly_accession,
              feature_ga,
              feature_start,
              feature_end,
              line_split[18]
            )
        )
        protein_nt_seqrecords.append(protein_nt_seqrecord)


      # add row to new feature table
      joined = '\t'.join(line_split)
      newft.append(joined)

  # add crispr summary row to newft
  crispr_summary = [
    'CRISPR array', crispr['crispr_type'], '', '', '', '', crispr['genomic_accession'], crispr['start'], crispr['end'], str(seq), '', '', '', '', crispr['crispr_id'], '', '', '', str(length), '', crispr['dr_consensus_seq']
  ]
  crispr_summary_joined = '\t'.join(crispr_summary)
  newft.append(crispr_summary_joined)

  # add drs and spacer rows to newft
  for el in crispr['elements']:
    el_start = int(el['start'])
    el_end = int(el['end'])
    el_seq = record.seq[el_start-1: el_end]  # actual DR seq, not consensus
    el_length = el_end - el_start + 1

    el_row = [
      el['element_type'], crispr['crispr_type'], '', '', '', '', crispr['genomic_accession'], el['start'], el['end'], str(el_seq), '', '', '', '', el['element_type'], '', '', '', str(el_length), '', ''
    ]
    el_row_joined = '\t'.join(el_row)
    newft.append(el_row_joined)

  # add header to top of newft
  header = header_orig.split('\t')
  header.insert(9, 'sequence')
  header_joined = '\t'.join(header)
  newft.insert(0, header_joined)

  # write new feature table to file
  newft_path = '{}{}__feature_table_new.txt'.format(extracted_dir, extracted_id)
  with open(newft_path, 'w') as outf:
    newft_joined = '\n'.join(newft)
    outf.write(newft_joined)

  print(f'Wrote new feature table to file {newft_path}')

  return protein_seqrecords, protein_nt_seqrecords


# # fetch protein from Entrez via non-redundant_refseq
# def fetchProteinSeq(refSeqId):
#   esearch = Entrez.esearch(db="protein", term=refSeqId)
#   esearch_read = Entrez.read(esearch)
#   uidList = esearch_read['IdList']

#   if len(uidList) > 0:
#     uid = uidList[0]
#     efetch_read = Entrez.efetch(
#       db = "protein",
#       id = uid,
#       retmode = "text",
#       rettype = "fasta"
#     ).read()

#     efetch_header = efetch_read.split('\n',1)[0]  # header is before first \n
#     efetch_seq = efetch_read.split('\n',1)[1].replace('\n','')  # sequence is after, clean up \n within
#     print(f'Fetched protein {refSeqId} from Entrez...')
#     return efetch_header, efetch_seq

#   else:
#     print(f"Error: couldn't fetch {refSeqId} from Entrez...")
#     return "couldn't fetch from Entrez", "couldn't fetch from Entrez"


# # write feature protein sequences to file
# def createProteinFile(seqrecords, extracted_id):
#   if not seqrecords: return
#   path = '{}{}__{}_proteins.fasta'.format(extracted_dir, extracted_id, len(seqrecords))
#   with open(path, 'w') as outf:
#     SeqIO.write(seqrecords, outf, 'fasta')
#   print(f'Wrote extracted region protein sequences to file {path}')

# # write feature protein nt sequences to file
# def createProteinNTFile(seqrecords, extracted_id):
#   if not seqrecords: return
#   path = '{}{}__{}_proteins_nt.fasta'.format(extracted_dir, extracted_id, len(seqrecords))
#   with open(path, 'w') as outf:
#     SeqIO.write(seqrecords, outf, 'fasta')
#   print(f'Wrote extracted region protein NT sequences to file {path}')


def main():

   for f in glob.glob(gff_dir + '*.gff'):
    # gff_filename = 'GCF_000006845-1_ASM684v1_genomic.gff'
    gff_filename = f.split('/')[-1]
    assembly_accession = gff_filename[:-len('_genomic.gff')]
    fna_filename = assembly_accession + '_genomic.fna'


    dict = parseGffToDict(gff_filename, assembly_accession)
    for crispr in dict['crisprs_list']:
      with open(genome_dir + fna_filename, 'rU') as file:
        record_dict = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))

      ex_id, ex_start, ex_end = extractSequence(crispr, record_dict, assembly_accession)
      protein_seqrecords, protein_nt_seqrecords = createNewFeatureTable(crispr, record_dict, ex_id, ex_start, ex_end, assembly_accession)
      createProteinFile(protein_seqrecords, ex_id)
      createProteinNTFile(protein_nt_seqrecords, ex_id)
      print('---')

# if __name__ == "__main__":
#     main()

# ############################################################################
# import re

# def update_csv(gff):
#     dr_test = re.compile('^\>')
#     keyword_test = re.compile(keyword)

#     with open(gff, 'r') as f:
#         line = f.readline()
#         while(line):
#             if header_test.match(line) and keyword_test.search(line):
#                 headers.append(line.rstrip())  # makes list of headers
#                 line = f.readline()
#                 reads.append(line.rstrip())  # makes list of corresponding reads
#             line = f.readline()  # exit loop


# gffInfoFields = ["seqid",
#                  "source",  # always CRISPRfinder
#                  "type",
#                  "start",
#                  "end",
#                  "score",  # .
#                  "strand",
#                  "phase",  # .
#                  "attributes"]

