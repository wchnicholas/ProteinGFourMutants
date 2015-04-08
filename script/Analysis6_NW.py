#!/usr/bin/python
import os
import sys
import glob
import matplotlib.pyplot as plt
import networkx as nx
import operator
from itertools import imap

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

def TsvWithHeader2Hash(fitfile):
  H = {}
  infile = open(fitfile,'r')
  countline = 0
  header = []
  for line in infile.xreadlines():
    countline += 1
    line = line.rstrip().rsplit("\t")
    if countline == 1: header = line; continue
    mut = line[0]
    H[mut] = {}
    for i in range(1,len(line)): H[mut][header[i]] = line[i]
  infile.close()
  return H

def filterfithash(fithash, condition):
  for mut in fithash.keys():
    if fithash[mut][condition] == 'NA' or '_' in mut:
      del fithash[mut]
  return fithash

def generatenodes(mut,WT,fithash):
  pathhash = {}
  nodes    = []
  for i in range(4): pathhash[i] = [WT[i],mut[i]]
  for i in range(16):
    index = str(bin(i)).rsplit('b')[1]
    index = ''.join(map(str,[0]*(4-len(index))))+index
    node  = ''
    for p in range(len(index)):
      node+=pathhash[p][int(index[p])]
    if fithash.has_key(node): nodes.append(node)
    else: nodes.append('NA')
  return list(set(nodes))

def buildgraph(nodes):
  G = nx.Graph()
  G.add_nodes_from(nodes)
  for n1 in nodes:
    for n2 in nodes:
      if hamming(n1,n2) == 1: G.add_edge(n1,n2)
  return G

def labelnode(var,fit):
  scale = int(round((2-fit)*35+30))
  if fit >= 2: return 'grey30'
  else: return 'grey'+str(scale)

def drawgraph(G,outfile,fithash,WT,condition):
  outfile=open(outfile,'w')
  outfile.write('strict digraph{'+"\n"+"\t"+'rankdir=LR'+"\n"+"\t"+'node [shape=box]'+"\n")
  for var in G.nodes():
    fit = float(fithash[var][condition])
    col = labelnode(var,fit)
    outfile.write("\t"+var+' [fillcolor='+col+', color=black, style="filled,rounded"];'+"\n")
  for E in G.edges():
    var1 = E[0]
    var2 = E[1]
    if hamming(WT,var1) < hamming(WT,var2):   s = var1; t = var2
    elif hamming(WT,var2) < hamming(WT,var1): s = var2; t = var1
    else: print 'Error in graph construction'; sys.exit()
    outfile.write("\t"+s+'->'+t+' [style=bold, color=black];'+"\n")
  outfile.write('}'+"\n")
  outfile.close()

def main():
  var_start  = sys.argv[1]
  var_end    = sys.argv[2]
  WT         = 'VDGV'
  fitfile    = 'result/Mutfit'
  if var_start == WT: outfile = 'xdot/WT2'+var_end+'.dot'
  else: outfile = 'xdot/'+var_start+'2'+var_end+'.dot'
  conditions = ['I10fit','I20fit','I90fit']
  condition  = 'I20fit'
  fithash    = TsvWithHeader2Hash(fitfile)
  fithash    = filterfithash(fithash, condition)
  nodes = generatenodes(var_end,var_start,fithash)
  G     = buildgraph(nodes)
  drawgraph(G,outfile,fithash,var_start,condition)
  os.system('dot -Tpng %s -o %s' % (outfile,outfile.replace('.dot','.png'))) 
  os.system('rm %s' % outfile)

if __name__ == '__main__':
  main()

