import sys

import pysam


def reduce_edge_quals(in_path, out_path, fmt, edge_length, qual_reduction):
    if fmt == "sam":
        in_mode = "r"
        out_mode = "w"
    elif fmt == "bam":
        in_mode = "rb"
        out_mode = "wb"
    elif fmt == ".cram":
        in_mode = "rc"
        out_mode = "wc"

    in_sam = pysam.AlignmentFile(in_path, in_mode)
    out_sam = pysam.AlignmentFile(out_path, out_mode, template=in_sam)

    for read in in_sam.fetch(until_eof=True):
        # print(type(read), read.qqual, read.qual, read.query_qualities)
        qualities = read.query_qualities
        for idx in range(edge_length):
            new_qual = qualities[idx] - qual_reduction
            if new_qual < 0:
                new_qual = 0
            qualities[idx] = new_qual

            idx = -1 - idx
            new_qual = qualities[idx] - qual_reduction
            if new_qual < 0:
                new_qual = 0
            qualities[idx] = new_qual

        # Note that to set quality scores the sequence has to be set beforehand as this will determine the expected length of the quality score array.
        # Assigning to this attribute will invalidate any quality scores.
        read.query_sequence = read.query_sequence
        read.query_qualities = qualities
        out_sam.write(read)


if __name__ == "__main__":
    if len(sys.argv) == 4:
        in_path = sys.argv[1]
        out_path = sys.argv[2]
        fmt = sys.argv[3]
        edge_length = 3
        qual_reduction = 20
        reduce_edge_quals(in_path, out_path, fmt, edge_length, qual_reduction)
    else:
        print("Usage: trim_quals.py in_path out_path format")
