################################################################################
# CS 224W (Fall 2017) - Project
# Weighted Sign Networks
# Author: vhying@stanford.edu
# Last Updated: Oct 29, 2017
################################################################################

import snap
import random
import numpy as np
import matplotlib.pyplot as plt


def loadSigns(filename):
    """
    :param - filename: undirected graph with associated edge sign

    return type: dictionary (key = node pair (a,b), value = sign)
    return: Return sign associated with node pairs. Both pairs, (a,b) and (b,a)
    are stored as keys. Self-edges are NOT included.
    """
    signs = {}
    with open(filename, 'r') as ipfile:
      for line in ipfile:
        if line[0] != '#':
            if filename[-3:] == "csv":
              line_arr = line.split(',')
            else:
              line_arr = line.split()
            if line_arr[0] == line_arr[1]:
              continue
            node1 = int(line_arr[0])
            node2 = int(line_arr[1])
            sign = float(line_arr[2])
            signs[(node1, node2)] = sign
            signs[(node2, node1)] = sign

    return signs


def selfEdgeDel(G):
    # Remove self-edges
    for edge in G.Edges():
      srcNId = edge.GetSrcNId()
      dstNId = edge.GetDstNId()
      if srcNId == dstNId: 
        G.DelEdge(srcNId, dstNId)

    return G


def computeBalanceTriads(G, signs):
    """
    :param - G: graph
    :param - signs: Dictionary of signs (key = node pair (a,b), value = sign)

    return type: List, each position representing count of t0, t1, t2, and t3, respectively.
    return: Return the counts for t0, t1, t2, and t3 triad types. Count each triad
    only once and do not count self edges.
    """
    triad_count = [0, 0, 0, 0] # each position represents count of t0, t1, t2, t3, respectively
    counted = set()

    for NI in G.Nodes():
        Id = NI.GetId()
        deg = NI.GetDeg()
        for nth in range(deg):
            nId = NI.GetNbrNId(nth)
            ndeg = G.GetNI(nId).GetDeg()
            for nnth in range(ndeg):
                nnId = G.GetNI(nId).GetNbrNId(nnth)
                if Id == nnId: 
                  continue

                if G.IsEdge(Id, nnId):
                  tup = tuple(sorted([Id, nId, nnId]))

                  # Check the Set
                  if (tup in counted):
                    continue

                  # Insert Triad into set
                  counted.add(tup)
    
                  # Count Triads
                  numPos = 0
                  if (signs[(Id, nId)] > 0):
                    numPos = numPos + 1
                  if (signs[(nId, nnId)] > 0):
                    numPos = numPos + 1
                  if (signs[(Id, nnId)] > 0):
                    numPos = numPos + 1

                  triad_count[numPos] = triad_count[numPos] + 1  

    return triad_count


def fracPosNeg(G, signs):
    """
    :param - G: graph
    :param - signs: Dictionary of signs (key = node pair (a,b), value = sign)

    Computes and prints the fraction of positive edges and negative edges,
        and the probability of each type of triad.
    """
    fracPos = 0
    fracNeg = 0

    for k,v in signs.iteritems():
      if v > 0:
        fracPos = fracPos + 1
      elif v < 0:
        fracNeg = fracNeg + 1

    fracPos = fracPos / 2
    fracNeg = fracNeg / 2
    total = fracPos + fracNeg
    fracPos = fracPos / float(total)
    fracNeg = fracNeg / float(total)

    print 'Fraction of Positive Edges: %0.4f' % (fracPos)
    print 'Fraction of Negative Edges: %0.4f' % (fracNeg)


def unbalancedTriads(G, signs):
    """
    :param - G: Graph
    :param - signs: Dictionary of signs (key = node pair (a,b), value = sign)

    return type: Boolean
    return: Returns whether G is balanced (True) or not (False).
    """
    isBalanced = False
    counted = set()

    for NI in G.Nodes():
        Id = NI.GetId()
        deg = NI.GetDeg()
        for nth in range(deg):
            nId = NI.GetNbrNId(nth)
            ndeg = G.GetNI(nId).GetDeg()
            for nnth in range(ndeg):
                nnId = G.GetNI(nId).GetNbrNId(nnth)
                if Id == nnId: 
                  continue

                if G.IsEdge(Id, nnId):
                  tup = tuple(sorted([Id, nId, nnId]))

                  # Check the Set
                  if (tup in counted):
                    continue

                  # Insert Triad into set
                  counted.add(tup)

    triads = list(counted)

    numUnbalanced = 0
    for idx in range(len(triads)):
      tup = triads[idx]

      e1 = (tup[0], tup[1])
      e2 = (tup[1], tup[2])
      e3 = (tup[0], tup[2])

      numNeg = 0
      if signs[e1] < 0:
        numNeg = numNeg + 1
      if signs[e2] < 0:
        numNeg = numNeg + 1
      if signs[e3] < 0:
        numNeg = numNeg + 1

      if (numNeg == 1) or (numNeg == 3):
        numUnbalanced = numUnbalanced + 1
      
    return (numUnbalanced, len(triads))


def incrementStatusTriad(G, signs, Id, nId, nnId, triad_count):
    X = G.GetNI(Id)
    A = G.GetNI(nId)
    B = G.GetNI(nnId)

    # t1, t2, t5, t6
    if X.IsInNId(nId) and X.IsOutNId(nnId):
      if (signs[(Id, nId)] > 0) and (signs[(Id, nnId)] > 0):
        triad_count[0] = triad_count[0] + 1  
      elif (signs[(Id, nId)] > 0) and (signs[(Id, nnId)] < 0):
        triad_count[1] = triad_count[1] + 1  
      elif (signs[(Id, nId)] < 0) and (signs[(Id, nnId)] > 0):
        triad_count[4] = triad_count[4] + 1  
      else:
        triad_count[5] = triad_count[5] + 1  
    # t3, t4, t7, t8
    elif X.IsInNId(nId) and X.IsInNId(nnId): 
      if (signs[(Id, nId)] > 0) and (signs[(Id, nnId)] > 0):
        triad_count[2] = triad_count[2] + 1  
      elif (signs[(Id, nId)] > 0) and (signs[(Id, nnId)] < 0):
        triad_count[3] = triad_count[3] + 1  
      elif (signs[(Id, nId)] < 0) and (signs[(Id, nnId)] > 0):
        triad_count[6] = triad_count[6] + 1  
      else:
        triad_count[7] = triad_count[7] + 1  
    # t9, t10, t13, t14
    elif X.IsOutNId(nId) and X.IsOutNId(nnId): 
      if (signs[(Id, nId)] > 0) and (signs[(Id, nnId)] > 0):
        triad_count[8] = triad_count[8] + 1  
      elif (signs[(Id, nId)] > 0) and (signs[(Id, nnId)] < 0):
        triad_count[9] = triad_count[9] + 1  
      elif (signs[(Id, nId)] < 0) and (signs[(Id, nnId)] > 0):
        triad_count[12] = triad_count[12] + 1  
      else:
        triad_count[13] = triad_count[13] + 1  
    # t11, t12, t15, t16
    elif X.IsOutNId(nId) and X.IsInNId(nnId): 
      if (signs[(Id, nId)] > 0) and (signs[(Id, nnId)] > 0):
        triad_count[10] = triad_count[10] + 1  
      elif (signs[(Id, nId)] > 0) and (signs[(Id, nnId)] < 0):
        triad_count[11] = triad_count[11] + 1  
      elif (signs[(Id, nId)] < 0) and (signs[(Id, nnId)] > 0):
        triad_count[14] = triad_count[14] + 1  
      else:
        triad_count[15] = triad_count[15] + 1  


def computeStatusTriads(G, signs):
    """
    :param - G: graph
    :param - signs: Dictionary of signs (key = node pair (a,b), value = sign)

    return type: List, each position representing count of t0, t1, t2, and t3, respectively.
    return: Return the counts for t0, t1, t2, and t3 triad types. Count each triad
    only once and do not count self edges.
    """
    # Each position represents count of t0, t1, t2, t3, ... respectively
    triad_count = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] 
    counted = set()

    for NI in G.Nodes():
        Id = NI.GetId()
        deg = NI.GetDeg()
        for nth in range(deg):
            nId = NI.GetNbrNId(nth)
            ndeg = G.GetNI(nId).GetDeg()
            for nnth in range(ndeg):
                nnId = G.GetNI(nId).GetNbrNId(nnth)
                if Id == nnId: 
                  continue

                if G.IsEdge(Id, nnId):
                  tup = tuple(sorted([Id, nId, nnId]))

                  # Check the Set
                  if (tup in counted):
                    continue

                  # Insert Triad into set
                  counted.add(tup)

                  # Count Triads
                  incrementStatusTriad(G, signs, Id, nId, nnId, triad_count)
                  incrementStatusTriad(G, signs, nnId, Id, nId, triad_count)
                  incrementStatusTriad(G, signs, nId, nnId, Id, triad_count)

    return triad_count


# Data Load
print "DATA LOAD"
print "Loading Unweighted Graphs:"
print "Loading epinions graph..."
epinionsGr = snap.TNGraph.New()
epinionsGr = snap.LoadEdgeList(snap.PNGraph, "Datasets/soc-sign-epinions.txt", 0, 1)
epinionsGr = selfEdgeDel(epinionsGr)
epinionsSigns = loadSigns("Datasets/soc-sign-epinions.txt")

print "Loading slashdot graph..."
slashdotGr = snap.TNGraph.New()
slashdotGr = snap.LoadEdgeList(snap.PNGraph, "Datasets/soc-sign-Slashdot090221.txt", 0, 1)
slashdotGr = selfEdgeDel(slashdotGr)
slashdotSigns = loadSigns("Datasets/soc-sign-Slashdot090221.txt")

#print "Loading wikipedia (requests for adminship) graph..."
#wikiGr = snap.TNGraph.New()
#wikiGr = snap.LoadEdgeList(snap.PNGraph, "Datasets/wiki-RfA.txt", 0, 1)

#print "Loading wikipedia (election) graph..."
#wikiGr = snap.TNGraph.New()
#wikiGr = snap.LoadEdgeList(snap.PNGraph, "Datasets/wikiElec.ElecBs3.txt", 0, 1)
#print

print "Loading Weighted Graphs:"
print "Loading wikipedia (requests for adminship) graph..."
wikiRFAGr = snap.TNGraph.New()
wikiRFAGr = snap.LoadEdgeList(snap.PNGraph, "Datasets/RFAnet.csv", 0, 1, ',')
wikiRFAGr = selfEdgeDel(wikiRFAGr)
wikiRFASigns = loadSigns("Datasets/RFAnet.csv")

#print "Loading wikipedia (edits) graph..."
#wikiEditsGr = snap.TNGraph.New()
#wikiEditsGr = snap.LoadEdgeList(snap.PNGraph, "Datasets/WikiSignedNet.csv", 0, 1, ',')
#wikiEditsGr = selfEdgeDel(wikiEditsGr)
#wikiEditsSigns = loadSigns("Datasets/WikiSignedNet.csv")

print "Loading bitcoin alpha graph..."
btcAlphaGr = snap.TNGraph.New()
btcAlphaGr = snap.LoadEdgeList(snap.PNGraph, "Datasets/soc-sign-bitcoinalpha.csv", 0, 1, ',')
btcAlphaGr = selfEdgeDel(btcAlphaGr)
btcAlphaSigns = loadSigns("Datasets/soc-sign-bitcoinalpha.csv")

print "Loading bitcoin otc graph..."
btcOTCGr = snap.TNGraph.New()
btcOTCGr = snap.LoadEdgeList(snap.PNGraph, "Datasets/soc-sign-bitcoinotc.csv", 0, 1, ',')
btcOTCGr = selfEdgeDel(btcOTCGr)
btcOTCSigns = loadSigns("Datasets/soc-sign-bitcoinotc.csv")
print


# Triad Count and Positive Negative edges
print "TRIAD COUNT"
print "Epinions graph:"
triad_count = computeBalanceTriads(epinionsGr, epinionsSigns)
total_triads = float(sum(triad_count)) if sum(triad_count) != 0 else 1
for i in range(4):
    print "Count of Triad t%d: %d" % (i, triad_count[i])
for i in range(4):
    print "Fraction of Triad t%d: %0.4f" % (i, triad_count[i]/total_triads)
fracPosNeg(epinionsGr, epinionsSigns)

print "Slashdot graph:"
triad_count = computeBalanceTriads(slashdotGr, slashdotSigns)
total_triads = float(sum(triad_count)) if sum(triad_count) != 0 else 1
for i in range(4):
    print "Count of Triad t%d: %d" % (i, triad_count[i])
for i in range(4):
    print "Fraction of Triad t%d: %0.4f" % (i, triad_count[i]/total_triads)

fracPosNeg(slashdotGr, slashdotSigns)

print "Wiki RFA graph:"
triad_count = computeBalanceTriads(wikiRFAGr, wikiRFASigns)
total_triads = float(sum(triad_count)) if sum(triad_count) != 0 else 1
for i in range(4):
    print "Count of Triad t%d: %d" % (i, triad_count[i])
for i in range(4):
    print "Fraction of Triad t%d: %0.4f" % (i, triad_count[i]/total_triads)
fracPosNeg(wikiRFAGr, wikiRFASigns)

#print "Wiki Edits graph:"
#triad_count = computeBalanceTriads(wikiEditsGr, wikiEditsSigns)
#total_triads = float(sum(triad_count)) if sum(triad_count) != 0 else 1
#for i in range(4):
#    print "Count of Triad t%d: %d" % (i, triad_count[i])
#for i in range(4):
#    print "Fraction of Triad t%d: %0.4f" % (i, triad_count[i]/total_triads)
#fracPosNeg(wikiEditsGr, wikiEditsSigns)

print "BTC Alpha graph:"
triad_count = computeBalanceTriads(btcAlphaGr, btcAlphaSigns)
total_triads = float(sum(triad_count)) if sum(triad_count) != 0 else 1
for i in range(4):
    print "Count of Triad t%d: %d" % (i, triad_count[i])
for i in range(4):
    print "Fraction of Triad t%d: %0.4f" % (i, triad_count[i]/total_triads)
fracPosNeg(btcAlphaGr, btcAlphaSigns)

print "BTC OTC graph:"
triad_count = computeBalanceTriads(btcOTCGr, btcOTCSigns)
total_triads = float(sum(triad_count)) if sum(triad_count) != 0 else 1
for i in range(4):
    print "Count of Triad t%d: %d" % (i, triad_count[i])
for i in range(4):
    print "Fraction of Triad t%d: %0.4f" % (i, triad_count[i]/total_triads)
fracPosNeg(btcOTCGr, btcOTCSigns)
print


# BALANCE
print "BALANCE"
print "Epinions graph:"
print "Unbalanced Triads = %d, Total Triads = %d" % unbalancedTriads(epinionsGr, epinionsSigns)

print "Slashdot graph:"
print "Unbalanced Triads = %d, Total Triads = %d" % unbalancedTriads(slashdotGr, slashdotSigns)

print "Wiki RFA graph:"
print "Unbalanced Triads = %d, Total Triads = %d" % unbalancedTriads(wikiRFAGr, wikiRFASigns)

#print "Wiki Edits graph:"
#print "Unbalanced Triads = %d, Total Triads = %d" % unbalancedTriads(wikiEditsGr, wikiEditsSigns)

print "BTC Alpha graph:"
print "Unbalanced Triads = %d, Total Triads = %d" % unbalancedTriads(btcAlphaGr, btcAlphaSigns)

print "BTC OTC graph:"
print "Unbalanced Triads = %d, Total Triads = %d" % unbalancedTriads(btcOTCGr, btcOTCSigns)
print


# STATUS
print "STATUS"
print "BTC Alpha graph:"
triad_count = computeStatusTriads(btcAlphaGr, btcAlphaSigns)
total_triads = float(sum(triad_count)) if sum(triad_count) != 0 else 1
for i in range(16):
    print "Count of Triad t%d: %d" % (i, triad_count[i])
for i in range(16):
    print "Fraction of Triad t%d: %0.4f" % (i, triad_count[i]/total_triads)

print "BTC OTC graph:"
triad_count = computeStatusTriads(btcOTCGr, btcOTCSigns)
total_triads = float(sum(triad_count)) if sum(triad_count) != 0 else 1
for i in range(16):
    print "Count of Triad t%d: %d" % (i, triad_count[i])
for i in range(16):
    print "Fraction of Triad t%d: %0.4f" % (i, triad_count[i]/total_triads)
print
