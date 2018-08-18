#python script to pick randomly the fasta sequences

from Bio import SeqIO
from random import sample
with open("/home/tarakaramji/Downloads/shuff.fa") as f:
    seqs = SeqIO.parse(f, "fasta")
    samps = ((seq.name, seq.seq) for seq in  sample(list(seqs),544))
    for samp in samps:
        print(">{}\n{}".format(*samp))
