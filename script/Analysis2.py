#!/usr/bin/python
import os
import sys
import operator
import colorsys
import networkx as nx
import numpy as np
from math import exp
from itertools import imap

def floor(fit):
  if fit == 'NA': return 'NA'
  elif float(fit) < 0.01: return 0.01
  else: return float(fit)

def all_shortest_paths(G, source, target, weight=None):
    if weight is not None:
        pred,dist = nx.dijkstra_predecessor_and_distance(G,source,weight=weight)
    else:
        pred = nx.predecessor(G,source)
    if target not in pred:
        raise nx.NetworkXNoPath()
    stack = [[target,0]]
    top = 0
    while top >= 0:
        node,i = stack[top]
        if node == source:
            yield [p for p,n in reversed(stack[:top+1])]
        if len(pred[node]) > i:
            top += 1
            if top == len(stack):
                stack.append([pred[node][i],0])
            else:
                stack[top] = [pred[node][i],0]
        else:
            stack[top-1][1] += 1
            top -= 1

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

def fillinmissing(fithash,missfitfile,condition):
  infile = open(missfitfile,'r')
  for line in infile.xreadlines():
    if 'genotype_missing' in line: continue
    line = line.rstrip().rsplit("\t")
    mut  = line[0]
    fit  = exp(float(line[1]))
    fithash[mut] = {}
    fithash[mut][condition] = fit
    print mut, fit
  infile.close()
  return fithash

def labelnode(var,fit):
  high = float(7.5)
  low  = float(0)
  mid  = float(high-low)/2+low
  if fit > high:  return colorsys.rgb_to_hsv(1, 0, 0)
  elif fit > mid: return colorsys.rgb_to_hsv(1,(high-fit)/(high-mid),0)
  elif fit > low: return colorsys.rgb_to_hsv(1,1,(mid-fit)/(mid-low))
  elif fit <= low: return colorsys.rgb_to_hsv(1,1,1)
  else: print 'Something is wrong with the coloring function'; print fit; sys.exit()

def drawgraph(G,outfile,fithash,WT,mut,condition):
  outfile=open(outfile,'w')
  outfile.write('strict digraph{'+"\n"+"\t"+'rankdir=LR'+"\n"+"\t"+'node [shape=box]'+"\n")
  for var in G.nodes():
    fit = float(fithash[var][condition])
    col = labelnode(var,fit)
    col = ','.join(map(str,list(col)))
    outfile.write("\t"+var+' [fillcolor="'+col+'", color=black, style="filled,rounded"];'+"\n")
  for E in G.edges():
    var1 = E[0]
    var2 = E[1]
    fit1 = float(fithash[var1][condition])
    fit2 = float(fithash[var2][condition])
    if hamming(var1,var2) == 1:
      if fit1 < fit2: s = var1; t = var2
      elif fit2 < fit1: s = var2; t = var1
      else: print 'Error in graph construction'; sys.exit()
    '''
    if hamming(WT,var1) < hamming(WT,var2):   s = var1; t = var2
    elif hamming(WT,var2) < hamming(WT,var1): s = var2; t = var1
    elif hamming(WT,var2) == hamming(WT,var1):
      if hamming(mut,var1) < hamming(mut,var2): s= var2; t = var1
      elif hamming(mut,var2) < hamming(mut,var1): s= var1; t = var2
      else: print 'Error in graph construction'; sys.exit()
    else: print 'Error in graph construction'; sys.exit()
    '''
    outfile.write("\t"+s+'->'+t+' [style=bold, color=black];'+"\n")
  outfile.write('}'+"\n")
  outfile.close()

def buildgraph(nodes,fcutoff,fithash,condition):
  G = nx.Graph()
  [G.add_node(node) for node in nodes if float(fithash[node][condition]) >= float(fcutoff)] 
  for n1 in G.nodes():
    for n2 in G.nodes(): 
      if hamming(n1,n2) == 1: G.add_edge(n1,n2)
  return G

def main():
  WT          = 'VDGV'
  var         = sys.argv[1]
  fitfile     = 'result/Mutfit'
  missfitfile = 'result/regression_missing'
  condition   = 'I20fit'
  fcutoff     = float(-1)
  fithash     = TsvWithHeader2Hash(fitfile)
  print "Total # of variants in the raw data: %d" % len(fithash.keys())
  fithash     = filterfithash(fithash,condition)
  print "Total # of variants pass filter of raw data: %d" % len(fithash.keys())
  fithash     = fillinmissing(fithash,missfitfile,condition)
  muts        = fithash.keys()
  print "Total # of variants after fill in with regression: %d" % len(muts)
  sys.exit() ############
  '''
  G           = buildgraph(muts,fcutoff,fithash,condition)
  print "%d variants pass cutoff and used for graph building" % len(G.nodes())
  if nx.has_path(G,WT,var): paths = all_shortest_paths(G,WT,var)
  else: print 'There is no path exist between %s and %s' % (WT, var); sys.exit()
  for path in paths: 
    pathfits = [float(fithash[step][condition]) for step in path]
    print path, min(pathfits)
  '''
  nodes       = ['PIWI', 'FIWI', 'FIWV', 'FIFV', 'WIFV', 'WWFV', 'WWFG', 'WWLG', 'FWLG']
  outfile     = 'xdot/WT2PIWI_custom.dot'
  custom_G    = buildgraph(nodes,fcutoff,fithash,condition)
  drawgraph(custom_G,outfile,fithash,WT,var,condition)
  os.system('dot -Tpng %s -o %s' % (outfile,outfile.replace('.dot','.png')))
  #os.system('rm %s' % outfile)
  
if __name__ == '__main__':
  main()
