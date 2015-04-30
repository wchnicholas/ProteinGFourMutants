#!/usr/bin/python
import os
import sys

def main():
  WT  = ['39V','40D','41G','54V']
  pos = ['39','40','41','54']
  aas = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
  singlemut = [p+aa for p in pos for aa in aas if p+aa not in WT]
  doublemut = []
  for i in pos: 
    for j in pos:
      for aa1 in aas:
        for aa2 in aas:
          m1 = i+aa1
          m2 = j+aa2
          dmut = m1+'-'+m2
          if m1 in WT or m2 in WT or int(i) >= int(j): continue
          doublemut.append(dmut)
  index = ['avg']+singlemut+doublemut
  infile  = open('result/regression_coef.txt','r')
  n       = 0
  for line in infile.xreadlines():
    print index[n]+"\t"+line.rstrip()
    n += 1
  infile.close()
    
 

if __name__ == '__main__':
  main()
