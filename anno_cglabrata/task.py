import argparse
import os

import anno_cglabrata



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ref",
        help="reference library of coding sequences of single-exon genes in fasta format",
        required=True
    )
    parser.add_argument(
        "--align",
        help="blastn alignment of single-exon genes to query genome in format 6",
        required=True
    )
    parser.add_argument(
        "--output_dir",
        help="location to store the output files",
        required=True
    )
    parser.add_argument(
        "--strain",
        help="Strain Identifier for gene ID in query genome",
        required=True
    )

    parser.add_argument(
        "--prefix",
        help="prefix of the output file",
        required=True
    )

    args = parser.parse_args()
    arguments = args.__dict__

    # Run the training job
    anno_cglabrata.annotate(arguments)
