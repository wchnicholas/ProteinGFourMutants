#!/usr/bin/python

def coloring(mut):
  p1  = int(mut.rsplit('-')[0][1:-1])
  p2  = int(mut.rsplit('-')[1][1:-1])
  poi = [38,39,40,53]
  if p1 in poi and p2 in poi: return 'red'
  else: return 'black'

def main():
  infile  = open('result/Epistasis','r')
  outfile = open('result/EpistasisCol','w')
  for line in infile.xreadlines():
    if 'Dmut' in line: outfile.write(line.rstrip()+"\t"+'color'+"\n")
    else: 
      line = line.rstrip().rsplit("\t")
      mut  = line[0] 
      epi  = line[1]
      col  = coloring(mut)
      outfile.write("\t".join([mut,epi,col])+"\n")
  outfile.close()

if __name__ == '__main__':
  main()

