#!/bin/python
import pandas as pd
import csv


def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1

    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)

    return kmers

# Insert genes sequence
Gene50mer=build_kmers('sequence',50)

with open('Gene50mer','w') as f:
    write=csv.writer(f)
    write.writerow(Gene50mer)
