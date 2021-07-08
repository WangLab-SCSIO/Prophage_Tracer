#!/opt/sysoft/Python-3.7.2/bin/python3
import sys
import numpy as np
from numpy.random import choice

from Bio import SeqIO
from Bio.Seq import Seq

def biased_DNA_generator(GCRatio,size):
    DNAstr="".join(np.random.choice(["A","T","C","G"], size, p=[0.5-GCRatio/2,0.5-GCRatio/2,GCRatio/2,GCRatio/2]))
    return DNAstr


print (sys.argv[1])
type (sys.argv[1])
print (sys.argv[2])
type (sys.argv[2])
GC = float(sys.argv[1])
S  = int(sys.argv[2])
DNA = biased_DNA_generator(GC,S)
print (DNA)

if __name__ == '__main__':
  biased_DNA_generator()
