#!/usr/bin/python
import os
import sys
import glob
import string
import operator
from string import atof
from itertools import imap
from Bio import SeqIO

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def rc(seq):
  complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
  rcseq = seq.translate(complements)[::-1]
  return rcseq

def hamming(str1, str2):
  assert len(str1) == len(str2)
  return sum(imap(operator.ne, str1, str2))

def callmut(refseq):
  pass

def readrecord(R1file,R2file,refseq,bchash,outfiles):
  R2handle = SeqIO.parse(R2file,"fastq")
  countread = 0
  for record_R1 in SeqIO.parse(R1file,"fastq"):
    countread += 1
    if countread%1000000 == 0: print 'finish processing %d reads' % countread
    record_R2 = R2handle.next()
    ID_R1   = record_R1.id
    seq_R1  = str(record_R1.seq)
    Qual_R1 = record_R1.letter_annotations["phred_quality"]
    ID_R2   = record_R2.id
    seq_R2  = rc(str(record_R2.seq))
    Qual_R2 = record_R2.letter_annotations["phred_quality"][::-1]
    assert(ID_R1==ID_R2)
    
    #Quality Check
    MID_R1_A   = seq_R1[0:3]
    MID_R1_B   = rc(seq_R1[-3::])
    if MID_R1_A != MID_R1_B: continue
    if MID_R1_A not in bchash.keys(): continue
    MID_R1     = MID_R1_A

    #Parse
    ROI_R1     = seq_R1[26:35]+seq_R1[71:74]
    ROI_R2     = seq_R2[26:35]+seq_R2[71:74]
    Qual_R1    = Qual_R1[26:35]+Qual_R1[71:74]
    Qual_R2    = Qual_R2[26:35]+Qual_R2[71:74]
    if min(Qual_R1) < 0 or min(Qual_R2) < 0: continue
    if ROI_R1==ROI_R2: 
      outfiles[MID_R1].write(ROI_R1+"\n")

def hashinref(infile):
  refhash = {}
  handle = open(infile, "rU")
  for record in SeqIO.parse(handle, "fasta"):
    refhash[str(record.id)] = str(record.seq)
  handle.close()
  return refhash

def hashinbarcode(infile):
  bchash = {}
  infile = open(infile,'r')
  for line in infile.xreadlines():
    line = line.rstrip().rsplit("\t")
    bchash[line[1]] = line[0]
  return bchash

def main():
  refhash = hashinref("Fasta/SeqInfo.fa") 
  bchash  = hashinbarcode('Fasta/Barcode.fa')
  refseq  = refhash['ProteinG_AO2']
  Fprime  = refhash['ForwardPrimer']
  filenames = sorted(glob.glob('fastq/*R1*.fastq'))
  #filenames = sorted(glob.glob('fastq/*R1*.fastq'))
  outfiles = {}
  for bc in bchash.keys(): outfiles[bc] = open('paired/'+bchash[bc],'w')
  countfile = 0
  for R1file in filenames:
    countfile += 1
    print 'working on %d PE file' % countfile
    R2file = R1file.replace('_R1_','_R2_')
    readrecord(R1file,R2file,refseq,bchash,outfiles)
  for bc in bchash.keys(): outfiles[bc].close()  

if __name__ == "__main__":
  main()

