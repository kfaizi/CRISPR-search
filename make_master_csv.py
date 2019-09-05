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
             ]
# left: P, S, U
    summary = None

    for name in names:
        name = name + NCBI_datestamp
        sumdir = Path(outdir, name)

        if summary is not None:
            summary = list(set(sumdir.glob("*summary_csv.csv")) - set(sumdir.glob(".*")))[0]  # exclude weird hidden binary excel files (artefacts from sshfs?)
            df_new = pd.read_csv(summary, dtype=object, index_col=False)
            df = pd.concat([df, df_new], axis=0)

        else:
            summary = list(set(sumdir.glob("*summary_csv.csv")) - set(sumdir.glob(".*")))[0]
            df = pd.read_csv(summary, dtype=object, index_col=False)

    with open(final, 'w') as f:
        df.to_csv(f, mode='w', index=False, header=True)


def winnow(path, out):
    df = pd.read_csv(path, dtype=object, index_col=False)
    df = df[df.DR_found == 'yes']  # discard hits w/o crisprs
    df = df[df.protein_length.astype(int) > 500]  # discard proteins < 500 aa
    df = df[df.distance_protein_and_array.astype(float) < 20000] # discard proteins outside 20kb flanking sequence
    # df = df[df.pident.astype(float) != 100]  # discard perfect hits

    # group rows by source id, sorted by % identity and crispr 'confidence'
    df = df.sort_values(['extracted_header', 'crispr_id', 'pident', 'crispr_type'], ascending=[True, True, False, True])

    # keep only the best match (by identity) for each hit
    # df2 = df.groupby('extracted_header', as_index=False).head(1).reset_index(drop=True)
    df2 = df.groupby('extracted_header', as_index=False).head(1)

    with open(out, 'w') as f:
        df.to_csv(f, mode='w', index=False, header=True)

def other_winnow(path, out):
    df = pd.read_csv(path, dtype=object, index_col=False)
    df = df.astype({'pident': float})

    # group rows by source id, sorted by % identity and crispr 'confidence'; keep only the best match (by identity) for each hit
    df = df.sort_values(['extracted_header', 'crispr_id', 'pident'], ascending=[True, True, False]).groupby('extracted_header', as_index=False).head(1)

    #.reset_index(drop=True)

    with open(out, 'w') as f:
        df.to_csv(f, mode='w', index=False, header=True)


vertical_stack()
# other_winnow(Path('/Users/kianfaizi/Desktop/sep5_all_orthologs.csv'), Path('/Users/kianfaizi/Desktop/sep5_best_orthologs.csv'))
