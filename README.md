# Mining metagenomes for new CRISPR effectors
This is a Python pipeline and CLI for identifying novel orthologs of CRISPR-Cas effectors from publicly available metagenome data.

Given a set of known protein sequences, it:
1. queries local metagenomic BLAST databases for similar ORFs using TBLASTN 🧬
2. searches significant contigs for putative CRISPR arrays using CRISPRFinder<sup id="a1">[1](#f1)</sup> 🔍
3. ranks, cleans, and summarizes the results for subsequent synthesis and experimental characterization. 🧪

In addition, it includes helper functions for:
- sorting, formatting, and preprocessing raw sequence data into searchable BLAST databases
- efficiently deduplicating sequences
- multithreading
- logging output
- basic job scheduling, including customizable SMS alerts whenever a run finishes or fails.

Developed during my Summer 2019 research in the [Hsu Lab](http://patrickhsulab.org) at the Salk Institute.

## References
<b id="f1">1.</b> Grissa I, Vergnaud G, Pourcel C. CRISPRFinder: a web tool to identify clustered regularly interspaced short palindromic repeats. _Nucleic Acids Res._ (2007):W52-7. https://doi.org/10.1093/nar/gkm360 [↩](#a1)
