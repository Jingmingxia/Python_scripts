#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
#from collections import OrderedDict
from Bio import SeqIO

def codon(files):
  #rna_file = sys.argv[1]
  len_list = []
  for record in SeqIO.parse(files,"fasta"):
    if len(record.seq) % 3 != 0:
      continue
    if "CDS" in record.description:
      len_list.append(len(record.seq)/3)
  return len_list

cds = codon(sys.argv[1])
num,max_cds,min_cds,avg_cds = len(cds),max(cds),min(cds),sum(cds)/len(cds)
with open (sys.argv[2],"w") as fh:
  fh.write("Total_cds\tLongest_cds\tShortest_cds\tAverage_cds\n")
  fh.write("%s\t%s\t%s\t%s\n" % (num,max_cds,min_cds,avg_cds))

window = 100
def plot_cds(cds,window):

    x = [(i+0.5)*window for i in range(int(max(cds)/window)+1)]
    y = [0 for _ in x]
    for i in cds:
        y[int(i/window)] += 1
    import matplotlib
    matplotlib.use("Agg")    
    from matplotlib import pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4.5))
    plt.subplots_adjust(top=0.95, left=0.10, right=0.95, bottom=0.10)
    ax.bar(x, y, width=window, linewidth=0.5, edgecolor="white", color="red")
    #plt.ylim([window*-0.1, plt.ylim()[1]])
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.xlabel("Protein length (aa)", fontsize=10, weight="bold")
    plt.ylabel("Number", weight="bold", fontsize=10,)
    plt.savefig("protein_length.pdf")
    plt.savefig("protein_length.png", dpi=900)

plot_cds(cds,window)
