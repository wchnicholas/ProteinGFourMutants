#!/usr/bin/python
import os
import sys
import random
import operator
import colorsys
import copy
import numpy as np
import networkx as nx
from math import exp
import itertools 
from networkx.utils import open_file, make_str

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

def connected_components(G):
    if G.is_directed():
        raise nx.NetworkXError("""Not allowed for directed graph G.
              Use UG=G.to_undirected() to create an undirected graph.""")
    seen={}
    components=[]
    for v in G:
        if v not in seen:
            c=nx.single_source_shortest_path_length(G,v)
            components.append(list(c.keys()))
            seen.update(c)
    components.sort(key=len,reverse=True)
    return components

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(map(operator.ne, str1, str2))

def TsvWithHeader2Hash(fitfile):
  H = {}
  infile = open(fitfile,'r')
  countline = 0
  header = []
  for line in infile.readlines():
    countline += 1
    line = line.rstrip().rsplit("\t")
    if countline == 1: header = line; continue
    mut = line[0]
    H[mut] = {}
    for i in range(1,len(line)): H[mut][header[i]] = line[i]
  infile.close()
  return H

def filterfithash(fithash, condition):
  newhash = {}
  for mut in fithash.keys():
    if fithash[mut][condition] != 'NA' and '_' not in mut:
      newhash[mut]=fithash[mut]
  return newhash

def fillinmissing(fithash,missfitfile,condition):
  infile = open(missfitfile,'r')
  for line in infile.readlines():
    if 'genotype_' in line: continue
    line = line.rstrip().rsplit("\t")
    mut  = line[0]
    fit  = exp(float(line[1]))
    fithash[mut] = {}
    fithash[mut][condition] = fit
  infile.close()
  return fithash

def labelnode(var,fit):
  high = float(2)
  low  = float(0)
  mid  = float(high-low)/2+low
  if fit > high:  return colorsys.rgb_to_hsv(1, 0, 0)
  elif fit > mid: return colorsys.rgb_to_hsv(1,(high-fit)/(high-mid),0)
  elif fit > low: return colorsys.rgb_to_hsv(1,1,(mid-fit)/(mid-low))
  elif fit <= low: return colorsys.rgb_to_hsv(1,1,1)
  else: 
    print ("Something is wrong with the coloring function")
    print (var, fit)
    sys.exit()

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
    if hamming(WT,var1) < hamming(WT,var2):   s = var1; t = var2
    elif hamming(WT,var2) < hamming(WT,var1): s = var2; t = var1
    elif hamming(WT,var2) == hamming(WT,var1):
      if hamming(mut,var1) < hamming(mut,var2): s= var2; t = var1
      elif hamming(mut,var2) < hamming(mut,var1): s= var1; t = var2
      else: print ('Error in graph construction'); sys.exit()
    else: print ('Error in graph construction'); sys.exit()
    outfile.write("\t"+s+'->'+t+' [style=bold, color=black];'+"\n")
  outfile.write('}'+"\n")
  outfile.close()

def buildgraph(nodes,fithash,condition):
  G = nx.Graph()
  [G.add_node(node) for node in nodes] 
  for n1 in G.nodes():
    for n2 in G.nodes(): 
      if hamming(n1,n2) == 1: G.add_edge(n1,n2)
  return G

def drawgraph(G,outfile,fithash,condition):
  outfile=open(outfile,'w')
  outfile.write('strict graph{'+"\n"+"\t"+'rankdir=LR'+"\n"+"\t"+'node [shape=box]'+"\n")
  for var in G.nodes():
    fit = float(fithash[var][condition])
    col = labelnode(var,fit)
    col = ','.join(map(str,list(col)))
    outfile.write("\t"+var+' [fillcolor="'+col+'", color=black, style="filled,rounded", fixedsize=shape, width= 0.1];'+"\n")
  for E in G.edges():
    var1 = E[0]
    var2 = E[1]
    outfile.write("\t"+var1+'--'+var2+' [style=bold, color=black];'+"\n")
  outfile.write('}'+"\n")
  outfile.close()

def main():
  WT          = 'VDGV'
  fitfile     = 'result/Mutfit'
  outfile     = 'xdot/SciMimic_WTandBenNoFillin.dot'
  missfitfile = 'result/regression_missing'
  #missfitfile = 'result/regression_all_WT'
  condition   = 'I20fit'
  fprate      = float(0.00)
  finalsize   = float(1659)
  fithash     = TsvWithHeader2Hash(fitfile)
  print ("Total # of variants in the raw data: %s") % len(fithash.keys())
  fithash  = filterfithash(fithash,condition)
  print ("Total # of variants pass filter of raw data: %d") % len(fithash.keys())
  #fithash  = fillinmissing(fithash,missfitfile,condition)
  #print ("Total # of variants after fill in with regression: %d") % len(fithash.keys())
  muts = [mut for mut in fithash.keys() if float(fithash[mut][condition]) >= float(fithash[WT][condition])]
  muts.remove('IGEV')
  muts.remove('IGQV')
  muts.remove('WNWY')
  #muts = fithash.keys()                                       #For mimicking science paper
  #fits = [float(fithash[mut][condition]) for mut in muts]     #For mimicking science paper
  #muts = [x for (y,x) in sorted(zip(fits,muts),reverse=True)] #For mimicking science paper
  #muts = muts[0:int(finalsize/(1-fprate))]                    #For mimicking science paper
  #random.shuffle(muts)                                        #For mimicking science paper
  #muts = muts[0:int(finalsize)]                               #For mimicking science paper
  print ("%d variants pass cutoff and used for graph building") % len(muts)
  G    = buildgraph(muts,fithash,condition)
  print ("Number of connected components: %d") % len(connected_components(G))
  drawgraph(G,outfile,fithash,condition)
  #os.system('dot -Tpng %s -o %s' % (outfile,outfile.replace('.dot','.png')))
  
if __name__ == '__main__':
  main()
