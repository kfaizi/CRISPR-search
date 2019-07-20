# modify for pathlib later
import pandas as pd
from Bio import SeqIO

parsed_dict = {'crisprs_list': [{'crispr_id': 'AAAA02021996.1_21523_21912', 'genomic_accession': 'AAAA02021996.1', 'crispr_type': 'CRISPR', 'start': '21523', 'end': '21912', 'dr_consensus_seq': 'TTCTGTATATTTTCGGACTTGTCCGAAA', 'dr_consensus_length': '28', 'elements': [{'element_type': 'CRISPRdr', 'start': '21523', 'end': '21550'}, {'element_type': 'CRISPRspacer', 'start': '21551', 'end': '21589', 'spacer_seq': 'CTCATTTCCGGACTTTCCGAAAACACATAGAACCAGATT', 'spacer_name': 'spacer_21551_39'}, {'element_type': 'CRISPRdr', 'start': '21590', 'end': '21617'}, {'element_type': 'CRISPRspacer', 'start': '21618', 'end': '21661', 'spacer_seq': 'GTGATTTTCGGACTTTCCGAAAAAGACTGCGAAGGCAAAAGTGG', 'spacer_name': 'spacer_21618_44'}, {'element_type': 'CRISPRdr', 'start': '21662', 'end': '21689'}, {'element_type': 'CRISPRspacer', 'start': '21690', 'end': '21734', 'spacer_seq': 'CTGATTTTCGGATTTTCCGAAAATCATCAGTAGAGTCAATTTCGC', 'spacer_name': 'spacer_21690_45'}, {'element_type': 'CRISPRdr', 'start': '21735', 'end': '21762'}, {'element_type': 'CRISPRspacer', 'start': '21763', 'end': '21807', 'spacer_seq': 'GTTATTTTCGGACTTTTCCGAGAACATCCAGAAGGATGTTGTTGG', 'spacer_name': 'spacer_21763_45'}, {'element_type': 'CRISPRdr', 'start': '21808', 'end': '21835'}, {'element_type': 'CRISPRspacer', 'start': '21836', 'end': '21884', 'spacer_seq': 'AGTTTTTGGAATATCCGAAAAATCCTTCGTTGACTTTGCTGTTGGTGCT', 'spacer_name': 'spacer_21836_49'}, {'element_type': 'CRISPRdr', 'start': '21885', 'end': '21912'}]}, {'crispr_id': 'AAAA02027088.1_723_794', 'genomic_accession': 'AAAA02027088.1', 'crispr_type': 'PossibleCRISPR', 'start': '723', 'end': '794', 'dr_consensus_seq': 'TATCAATTTGACTCTAATTTTTT', 'dr_consensus_length': '23', 'elements': [{'element_type': 'CRISPRdr', 'start': '723', 'end': '745'}, {'element_type': 'CRISPRspacer', 'start': '746', 'end': '771', 'spacer_seq': 'TTAGTTTACATATTTTAATATCTAAC', 'spacer_name': 'spacer_746_26'}, {'element_type': 'CRISPRdr', 'start': '772', 'end': '794'}]}, {'crispr_id': 'AAAA02027088.1_13954_14066', 'genomic_accession': 'AAAA02027088.1', 'crispr_type': 'PossibleCRISPR', 'start': '13954', 'end': '14066', 'dr_consensus_seq': 'AAGCAAATAAAGCACGCTCTCAAATTT', 'dr_consensus_length': '27', 'elements': [{'element_type': 'CRISPRdr', 'start': '13954', 'end': '13980'}, {'element_type': 'CRISPRspacer', 'start': '13981', 'end': '14039', 'spacer_seq': 'CAGCAATTGATGACAAACTATTTGAGATATTATCTTTTGGTAAGGTTTAAATTCAGTCT', 'spacer_name': 'spacer_13981_59'}, {'element_type': 'CRISPRdr', 'start': '14040', 'end': '14066'}]}, {'crispr_id': 'AAAA02022349.1_23281_23389', 'genomic_accession': 'AAAA02022349.1', 'crispr_type': 'PossibleCRISPR', 'start': '23281', 'end': '23389', 'dr_consensus_seq': 'TGATGTGATGGAAAGTTGAAAGTTTGGA', 'dr_consensus_length': '28', 'elements': [{'element_type': 'CRISPRdr', 'start': '23281', 'end': '23308'}, {'element_type': 'CRISPRspacer', 'start': '23309', 'end': '23361', 'spacer_seq': 'AAAAAAACTTTGGAACTAAATAGGGCCTGTGTAGGAAAGTTTTGGATGTGATA', 'spacer_name': 'spacer_23309_53'}, {'element_type': 'CRISPRdr', 'start': '23362', 'end': '23389'}]}, {'crispr_id': 'AAAA02006613.1_10366_10470', 'genomic_accession': 'AAAA02006613.1', 'crispr_type': 'PossibleCRISPR', 'start': '10366', 'end': '10470', 'dr_consensus_seq': 'ATGAAATCCATTTTTACCAAACCA', 'dr_consensus_length': '24', 'elements': [{'element_type': 'CRISPRdr', 'start': '10366', 'end': '10389'}, {'element_type': 'CRISPRspacer', 'start': '10390', 'end': '10446', 'spacer_seq': 'TTTTAATTATTACCATGAAATCCCTGCGGTACGGAGCTCGATTCTCTCTGGGTTGCC', 'spacer_name': 'spacer_10390_57'}, {'element_type': 'CRISPRdr', 'start': '10447', 'end': '10470'}]}, {'crispr_id': 'AAAA02026716.1_8899_8996', 'genomic_accession': 'AAAA02026716.1', 'crispr_type': 'PossibleCRISPR', 'start': '8899', 'end': '8996', 'dr_consensus_seq': 'CGCCGAGCCGCCGCCGTAGCCGT', 'dr_consensus_length': '23', 'elements': [{'element_type': 'CRISPRdr', 'start': '8899', 'end': '8921'}, {'element_type': 'CRISPRspacer', 'start': '8922', 'end': '8973', 'spacer_seq': 'GCCGCCCGCCGACGCCATCGACGCCGCTACGGAGCAGGGGCCGTAGCTGCTA', 'spacer_name': 'spacer_8922_52'}, {'element_type': 'CRISPRdr', 'start': '8974', 'end': '8996'}]}, {'crispr_id': 'AAAA02005981.1_22110_22187', 'genomic_accession': 'AAAA02005981.1', 'crispr_type': 'PossibleCRISPR', 'start': '22110', 'end': '22187', 'dr_consensus_seq': 'TTACTACTATAGTACAAATGCATTT', 'dr_consensus_length': '25', 'elements': [{'element_type': 'CRISPRdr', 'start': '22110', 'end': '22134'}, {'element_type': 'CRISPRspacer', 'start': '22135', 'end': '22162', 'spacer_seq': 'ATCTATTTAACTTTATGCATTGTTCTCT', 'spacer_name': 'spacer_22135_28'}, {'element_type': 'CRISPRdr', 'start': '22163', 'end': '22187'}]}, {'crispr_id': 'AAAA02005981.1_27204_27318', 'genomic_accession': 'AAAA02005981.1', 'crispr_type': 'PossibleCRISPR', 'start': '27204', 'end': '27318', 'dr_consensus_seq': 'GGAGTCGGACTTTTTTATTTTTT', 'dr_consensus_length': '23', 'elements': [{'element_type': 'CRISPRdr', 'start': '27204', 'end': '27226'}, {'element_type': 'CRISPRspacer', 'start': '27227', 'end': '27244', 'spacer_seq': 'TAAAAGGGAATGGATATA', 'spacer_name': 'spacer_27227_18'}, {'element_type': 'CRISPRdr', 'start': '27245', 'end': '27267'}, {'element_type': 'CRISPRspacer', 'start': '27268', 'end': '27295', 'spacer_seq': 'ATTTTTCGCTATTGTTTTTGTTTTGACG', 'spacer_name': 'spacer_27268_28'}, {'element_type': 'CRISPRdr', 'start': '27296', 'end': '27318'}]}, {'crispr_id': 'AAAA02026886.1_5066_5145', 'genomic_accession': 'AAAA02026886.1', 'crispr_type': 'PossibleCRISPR', 'start': '5066', 'end': '5145', 'dr_consensus_seq': 'TATTGCTTTTATTATTATGTCTATAT', 'dr_consensus_length': '26', 'elements': [{'element_type': 'CRISPRdr', 'start': '5066', 'end': '5091'}, {'element_type': 'CRISPRspacer', 'start': '5092', 'end': '5119', 'spacer_seq': 'ATAAATAGGAGACACTATGCAGTGTAAA', 'spacer_name': 'spacer_5092_28'}, {'element_type': 'CRISPRdr', 'start': '5120', 'end': '5145'}]}, {'crispr_id': 'AAAA02016417.1_14766_14847', 'genomic_accession': 'AAAA02016417.1', 'crispr_type': 'PossibleCRISPR', 'start': '14766', 'end': '14847', 'dr_consensus_seq': 'TACTCCTTCGTATTTAGCCCCGGTT', 'dr_consensus_length': '25', 'elements': [{'element_type': 'CRISPRdr', 'start': '14766', 'end': '14790'}, {'element_type': 'CRISPRspacer', 'start': '14791', 'end': '14822', 'spacer_seq': 'GGACGGCCATTACAGGTTGTCATAAGGAACTC', 'spacer_name': 'spacer_14791_32'}, {'element_type': 'CRISPRdr', 'start': '14823', 'end': '14847'}]}, {'crispr_id': 'AAAA02018790.1_16953_17048', 'genomic_accession': 'AAAA02018790.1', 'crispr_type': 'PossibleCRISPR', 'start': '16953', 'end': '17048', 'dr_consensus_seq': 'CAAGGCCACTAGATCTAATGTTG', 'dr_consensus_length': '23', 'elements': [{'element_type': 'CRISPRdr', 'start': '16953', 'end': '16975'}, {'element_type': 'CRISPRspacer', 'start': '16976', 'end': '17025', 'spacer_seq': 'GCGAGCTCGTGGCCATGGGATCCTGCATCGCGATGGCTGAAAGAGAGGAC', 'spacer_name': 'spacer_16976_50'}, {'element_type': 'CRISPRdr', 'start': '17026', 'end': '17048'}]}, {'crispr_id': 'AAAA02014824.1_9290_9397', 'genomic_accession': 'AAAA02014824.1', 'crispr_type': 'PossibleCRISPR', 'start': '9290', 'end': '9397', 'dr_consensus_seq': 'CGTCCTGCTGCTTCCGCCGCATCC', 'dr_consensus_length': '24', 'elements': [{'element_type': 'CRISPRdr', 'start': '9290', 'end': '9313'}, {'element_type': 'CRISPRspacer', 'start': '9314', 'end': '9373', 'spacer_seq': 'ACATCCTTCGCCGCCGCCGCCGCCGCCTGCGCCGATCCCCCCATCCATCCCGCGAGCGAG', 'spacer_name': 'spacer_9314_60'}, {'element_type': 'CRISPRdr', 'start': '9374', 'end': '9397'}]}, {'crispr_id': 'AAAA02007038.1_5522_5633', 'genomic_accession': 'AAAA02007038.1', 'crispr_type': 'PossibleCRISPR', 'start': '5522', 'end': '5633', 'dr_consensus_seq': 'CTGTAACAATTTATTATGTTATATCA', 'dr_consensus_length': '26', 'elements': [{'element_type': 'CRISPRdr', 'start': '5522', 'end': '5547'}, {'element_type': 'CRISPRspacer', 'start': '5548', 'end': '5607', 'spacer_seq': 'GTTTATTGCAAAAAGATCCTAATGTTCAGGTGAGAAAGATATAACCAATTCTGATGTGCT', 'spacer_name': 'spacer_5548_60'}, {'element_type': 'CRISPRdr', 'start': '5608', 'end': '5633'}]}, {'crispr_id': 'AAAA02008053.1_59432_59521', 'genomic_accession': 'AAAA02008053.1', 'crispr_type': 'PossibleCRISPR', 'start': '59432', 'end': '59521', 'dr_consensus_seq': 'CTGGAAACTGCAATGTGGTGCTG', 'dr_consensus_length': '23', 'elements': [{'element_type': 'CRISPRdr', 'start': '59432', 'end': '59454'}, {'element_type': 'CRISPRspacer', 'start': '59455', 'end': '59498', 'spacer_seq': 'TCAACTGTGGCTTCTGTCTCCATAATTACAGTGCTTGTTTCTCC', 'spacer_name': 'spacer_59455_44'}, {'element_type': 'CRISPRdr', 'start': '59499', 'end': '59521'}]}], 'num_crisprs': 14}
inputfile = "/Users/kianfaizi/temp/testAcsv.csv"


def update_test():
    sections = "10 sseqid slen sstart send sseq length sframe qseqid ppos pident qseq bitscore evalue saccver"  # from blast_formatter
    sections_list = sections.split()
    sections_list.pop(0) # remove "10" (csv filetype specifier) from sections_list

    # make blast csv into dataframe
    df_blast = pd.read_csv(inputfile, dtype=object, header=None, names=sections_list, index_col=None)

    # make parsed crisprfinder results into dataframe
    d_crispr = parsed_dict['crisprs_list']
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
                     'slen',
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
                       'crispr_type',
                       'crispr_id',
                       'DR_found',
                       'num_DRs',
                       'dr_consensus_seq',
                       'dr_consensus_length',
                       'slen',
                       'start',
                       'end',
                       'elements']]

    # rename columns more intuitively
    summary = summary.rename({'sseq': 'protein_sequence',
                              'slen': 'contig_length',
                              'qseqid': 'ortholog_match',
                              'dr_consensus_seq': 'DR_consensus_seq',
                              'dr_consensus_length': 'DR_consensus_length',
                              'start': 'array_start',
                              'end': 'array_end',
                              'elements': 'array'},
                             axis=1)

    # parse contigs for subsequent extraction, without overloading memory
    record_dict = SeqIO.index('/Users/kianfaizi/dev/output/A/in/A_dedupe.fa', 'fasta')

    # extract up to +/- 20kb flanking sequence around crispr locus
    extracted_list = []

    for row in summary.itertuples():
        extracted_raw = None
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
                extracted_raw = record_dict.get_raw(acc).decode()
                extracted_name = crispr_id + "_extracted"
                extracted_path = f"/Users/kianfaizi/dev/output/A/A_crisprs/{extracted_name}"  # update with Pathlib later

                # write extracted sequence to a new fasta file...
                with open(extracted_path, 'w+') as outfile:
                    SeqIO.write(extracted_seq, outfile, 'fasta')

                print(f"Wrote extracted CRISPR locus to {extracted_path}")

        except Exception:
            pass

        extracted_list.append(extracted_raw)

    # ...and add it to the summary dataframe
    summary.loc[:, 'extracted_sequence'] = extracted_list

    # write summary dataframe to csv
    with open('/Users/kianfaizi/temp/new_summary.csv', 'w') as f:
        summary.to_csv(f, mode='w', index=True, header=True)
