#!/usr/bin/env python
import sys
import re
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
#parser.add_argument("--label1", help="fraglabel")
#parser.add_argument("--label2", help="annolabel")
parser.add_argument("crop_file", help="cprops file")
parser.add_argument("fasta_file", help="fasta file")

args = parser.parse_args()

crop_file = args.crop_file
fasta_file = args.fasta_file

record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

fasta_id_list = []
# record the position of contigs
pos_dict = {}

for line in open(crop_file, "r"):
    items = re.split(":::|\s",line.strip())
    contig_id = items[0]
    contig_len = int(items[-1])

    if contig_id not in record_dict:
        continue

    if len(items) == 3:

        fasta_id_list.append(contig_id)

        seq = str(record_dict[contig_id].seq)
        print(">{}\n{}".format(contig_id, seq))
    else:

        if contig_id not in fasta_id_list:
            pos_dict[contig_id] = [contig_len]
            fasta_id_list.append(contig_id)
            seq_start = 0
            seq_end = seq_start + contig_len
            seq = str(record_dict[contig_id].seq)[seq_start:seq_end]
            if len(items) == 4:
                print(">{}:::{}\n{}".format(contig_id, items[1], seq))
            else:
                print(">{}:::{}:::{}\n{}".format(contig_id, items[1], items[2], seq))
        else:
            seq_start = 0 + sum(pos_dict[contig_id])
            seq_end = seq_start + contig_len
            pos_dict[contig_id].append(contig_len)
            seq = str(record_dict[contig_id].seq)[seq_start:seq_end]
            if len(items) == 4:
                print(">{}:::{}\n{}".format(contig_id, items[1], seq))
            else:
                print(">{}:::{}:::{}\n{}".format(contig_id, items[1], items[2], seq))
