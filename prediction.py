################################################################################
# CS 224W (Fall 2017) - Project
# Weighted Sign Networks
# Author: akshayk@stanford.edu, vhying@stanford.edu
# Last Updated: Dec 06, 2017
################################################################################

import snap
import random
import numpy as np
import matplotlib.pyplot as plt
import powerlaw
from sklearn.model_selection import cross_val_score
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import VotingClassifier


def loadDirectedSigns(filename):
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
            sign = int(line_arr[2])
            signs[(node1, node2)] = sign

    return signs

def selfEdgeDel(G):
    # Remove self-edges
    for edge in G.Edges():
      srcNId = edge.GetSrcNId()
      dstNId = edge.GetDstNId()
      if srcNId == dstNId: 
        G.DelEdge(srcNId, dstNId)

    return G

def ModSign(A):
    B = {}
    for i in A.keys():
        B[i] = abs(A[i])
    return B

def NormalizeSign(A):
    M = max(ModSign(A).values())
    B = {}
    for i in A.keys():
        B[i] = 1.0 * A[i] / M
    return B

def UnweightedSign(A):
    B = {}
    for i in A.keys():
        if A[i] > 0:
            B[i] = 1
        else:
            B[i] = -1
    return B


# Data Load
print "Loading Weighted Graphs:"
print "Loading bitcoin alpha graph..."
btcAlphaGr = snap.TNGraph.New()
btcAlphaGr = snap.LoadEdgeList(snap.PNGraph, "Datasets/soc-sign-bitcoinalpha.csv", 0, 1, ',')
btcAlphaGr = selfEdgeDel(btcAlphaGr)
btcAlphaSigns = loadDirectedSigns("Datasets/soc-sign-bitcoinalpha.csv")

print "Loading bitcoin otc graph..."
btcOTCGr = snap.TNGraph.New()
btcOTCGr = snap.LoadEdgeList(snap.PNGraph, "Datasets/soc-sign-bitcoinotc.csv", 0, 1, ',')
btcOTCGr = selfEdgeDel(btcOTCGr)
btcOTCSigns = loadDirectedSigns("Datasets/soc-sign-bitcoinotc.csv")
print


# Graph nodes and edges
numNodes = btcAlphaGr.GetNodes()
numEdges = btcAlphaGr.GetEdges()
maxId = btcAlphaGr.GetMxNId()

# Compose training set for pair (feature vector, w(e)).
# Encode each edge as a 2n dimensional feature vector.
train_features = []
train_target = [] 

i = 0
for edge in btcAlphaGr.Edges():
  fVec = [0] * (maxId * maxId)

  x, y = edge.GetId()
  idx = (x - 1) * numNodes + y
  fVec[idx] = 1
  w = btcAlphaSigns[(x, y)]
  
  train_features.append(fVec)
  train_target.append(w)

  i = i + 1
  if i > 500:
    break

# Train a classification model.
X = train_features
y = train_target 

# Hard Voting
clf1 = LogisticRegression(random_state=1)
clf2 = RandomForestClassifier(random_state=1)
clf3 = GaussianNB()

eclf = VotingClassifier(estimators=[('lr', clf1), ('rf', clf2), ('gnb', clf3)], voting='hard')

for clf, label in zip([clf1, clf2, clf3, eclf], ['Logistic Regression', 'Random Forest', 'naive Bayes', 'Ensemble']):
  scores = cross_val_score(clf, X, y, cv=5, scoring='accuracy')
  print("Accuracy: %0.2f (+/- %0.2f) [%s]" % (scores.mean(), scores.std(), label))

#Soft Voting
"""
clf1 = DecisionTreeClassifier(max_depth=4)
clf2 = KNeighborsClassifier(n_neighbors=7)
clf3 = SVC(kernel='rbf', probability=True)

eclf = VotingClassifier(estimators=[('lr', clf1), ('rf', clf2), ('gnb', clf3)], voting='soft', weights=[2,1,2])

clf1 = clf1.fit(X,y)
clf2 = clf2.fit(X,y)
clf3 = clf3.fit(X,y)
eclf = eclf.fit(X,y)
"""


# Predict signed weights of remaining directed edges.

# Number of training examples is m.

# Look at clusters (communities) in the network.

# Prune away non-relevant features for a given edge.

