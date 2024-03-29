# make master csv (concatenate all summaries)
import pandas as pd
from pathlib import Path


def vertical_stack():
    outdir = Path("/mnt/md0/kfaizi/search/output/")
    NCBI_datestamp = "_jul25"
    final = Path("/mnt/md0/kfaizi/search/master.csv")

    names = ['A',
             'AAA',
             'AAB',
             'AAC',
             'AAD',
             'AAE',
             'AAF',
             'B',
             'C',
             'CAA',
             'CAB',
             'F',
             'J',
             'L',
             'M',
             'N',
             'O',
             'Q',
             'R',
             'V',
             'P',
             'S',
             'U',
             ]

    summary = None
    for name in names:
        name = name + NCBI_datestamp
        sumdir = Path(outdir, name)

        if summary is not None:
            summary = list(set(sumdir.glob("*summary_csv.csv")) - set(sumdir.glob(".*")))[0]  # exclude weird hidden binary excel files (from sshfs?)
            df_new = pd.read_csv(summary, dtype=object, index_col=False)
            df = pd.concat([df, df_new], axis=0)

        else:
            summary = list(set(sumdir.glob("*summary_csv.csv")) - set(sumdir.glob(".*")))[0]
            df = pd.read_csv(summary, dtype=object, index_col=False)

    with open(final, 'w') as f:
        df.to_csv(f, mode='w', index=False, header=True)


def winnow(path, out_csv, out_fasta):
    df = pd.read_csv(path, dtype=object, index_col=False)
    df = df[df.DR_found == 'yes']  # discard hits w/o crisprs
    df = df[df.protein_length.astype(int) > 500]  # discard proteins < 500 aa
    df = df[df.distance_protein_and_array.astype(float) < 20000]  # discard proteins outside 20kb flanking sequence
    df = df[df.pident.astype(float) < 90]  # discard most similar hits
###### NEED TO FIX FILTERS > NOT STRINGENT ENOUGH > STILL GETTING 1000s OF HITS
    # group rows by source id, sorted by % identity and crispr 'confidence', keeping only best
    df = df.sort_values(['extracted_header', 'pident', 'crispr_type'], ascending=[True, False, True]).groupby('extracted_header', as_index=False).head(1)

    # write filtered dataframe to new csv
    with open(out_csv, 'w') as f:
        df.to_csv(f, mode='w', index=False, header=True)

    # write filtered proteins to a fasta file
    with open(out_fasta, 'w') as f:
        for row in df.itertuples():
            f.write(">" + row.extracted_header + "\n")
            f.write(row.protein_sequence + "\n")


# vertical_stack()
# winnow(Path("/mnt/md0/kfaizi/search/master.csv"), Path('/mnt/md0/kfaizi/search/filter_master.csv'),  Path('/mnt/md0/kfaizi/search/filter_proteins.fasta'))


def fix(path, out_fasta):
    df = pd.read_csv(path, dtype=object, index_col=False)

    with open(out_fasta, 'w') as f:
        for row in df.itertuples():
            f.write(">" + row.extracted_header + "\n")
            f.write(row.protein_sequence + "\n")

fix(Path("/Users/kianfaizi/Desktop/best_test.csv"), Path("/Users/kianfaizi/Desktop/best_prots.fasta"))
