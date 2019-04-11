
# coding: utf-8

# In[2]:


import snap
import random
import numpy as np
import matplotlib.pyplot as plt
import powerlaw


# In[3]:


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
            sign = int(line_arr[2])
            signs[(node1, node2)] = sign
            signs[(node2, node1)] = sign

    return signs


# In[53]:


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


# In[4]:


def selfEdgeDel(G):
    # Remove self-edges
    for edge in G.Edges():
      srcNId = edge.GetSrcNId()
      dstNId = edge.GetDstNId()
      if srcNId == dstNId: 
        G.DelEdge(srcNId, dstNId)

    return G


# In[5]:


def GetUnweightedDegreeDistribution(G):
    Deg = [N.GetDeg() for N in G.Nodes()]
    DegDis = [0 for i in range(max(Deg) + 1)]
    for i in Deg:
        DegDis[i] += 1

    DegDisX = []
    DegDisY = []
    for i in range(len(DegDis)):
        if DegDis[i] > 0:
            DegDisX.append(i)
            DegDisY.append(DegDis[i])

    results = powerlaw.Fit(DegDis)
    print(results.power_law.alpha)
    print(results.power_law.xmin)

    return DegDisX, DegDisY


# In[48]:


GUU = selfEdgeDel(snap.LoadEdgeList(snap.PUNGraph, "Datasets/soc-sign-bitcoinotc.csv", 0, 1, ','))
UUX, UUY = GetUnweightedDegreeDistribution(GUU)

plt.loglog(UUX, UUY, color = 'r', label = "Degree distribution for undirected unweighted")
plt.xlabel('Degree')
plt.ylabel('Count of degree')
plt.title('Degree Distribution')
plt.legend()
plt.show()
plt.show()


# In[49]:


GDU = selfEdgeDel(snap.LoadEdgeList(snap.PNGraph, "Datasets/soc-sign-bitcoinotc.csv", 0, 1, ','))
DUX, DUY = GetUnweightedDegreeDistribution(GDU)

plt.loglog(DUX, DUY, color = 'g', label = "Degree distribution for directed unweighted")
plt.xlabel('Degree')
plt.ylabel('Count of degree')
plt.title('Degree Distribution')
plt.legend()
plt.show()
plt.show()


# In[6]:


def GetUndirectedWeightedDegreeDistribution(G, signs):
    PDeg = []
    NDeg = []
    for N in G.Nodes():
        Degree = N.GetDeg()
        weightP = 0
        weightN = 0
        for i in range(Degree):
            weight = signs[(N.GetId(), N.GetNbrNId(i))]
            if weight >= 0:
                weightP += weight
            else:
                weightN -= weight
        PDeg.append(weightP)
        NDeg.append(weightN)
    
    minPW = min(PDeg)
    maxPW = max(PDeg)
    results = powerlaw.Fit(PDeg)
    print(results.power_law.alpha)
    print(results.power_law.xmin)
    DegDisPY, DegDisPX = np.histogram(PDeg, (maxPW - minPW + 1)/5)
    
    minNW = min(NDeg)
    maxNW = max(NDeg)
    results = powerlaw.Fit(NDeg)
    print(results.power_law.alpha)
    print(results.power_law.xmin)
    DegDisNY, DegDisNX = np.histogram(NDeg, (maxNW - minNW + 1)/5)

            
    return [DegDisPX[:-1], DegDisPY, DegDisNX[:-1], DegDisNY]


# In[76]:


GD = selfEdgeDel(snap.LoadEdgeList(snap.PNGraph, "Datasets/soc-sign-bitcoinotc.csv", 0, 1, ','))
signs = loadSigns("Datasets/soc-sign-bitcoinotc.csv")

R = GetUndirectedWeightedDegreeDistribution(GD, signs)
plt.loglog(R[0], R[1], color = 'g', label = "Positive weight distribution for undirected weighted")
plt.xlabel('Weight')
plt.ylabel('Count of weight')
plt.title('Weight Distribution')
plt.legend()
plt.show()

plt.loglog(R[2], R[3], color = 'r', label = "Negative weight distribution for undirected weighted")
plt.xlabel('Weight')
plt.ylabel('Count of weight')
plt.title('Weight Distribution')
plt.legend()
plt.show()


# In[7]:


def fitPowerLaw(A):
    results = powerlaw.Fit(A)
    print "alpha = " + str(results.power_law.alpha)
    print "x_min = " + str(results.power_law.xmin)


# In[8]:


def GetDirectedWeightedDegreeDistribution(G, signs):
    PIDeg = []
    PODeg = []
    NIDeg = []
    NODeg = []
    
    for N in G.Nodes():
        OutDegree = N.GetOutDeg()
        weightP = 0
        weightN = 0
        for i in range(OutDegree):
            weight = signs[(N.GetId(), N.GetOutNId(i))]
            if weight >= 0:
                weightP += weight
            else:
                weightN -= weight
        PODeg.append(weightP)
        NODeg.append(weightN)
        
        InDegree = N.GetInDeg()
        weightP = 0
        weightN = 0
        for i in range(InDegree):
            weight = signs[(N.GetInNId(i), N.GetId())]
            if weight >= 0:
                weightP += weight
            else:
                weightN -= weight
        PIDeg.append(weightP)
        NIDeg.append(weightN)
    
    print "Directed weighted signed positive links"
    YO, XO = np.histogram(PODeg, (max(PODeg) - min(PODeg) +1)/5)
    YI, XI = np.histogram(PIDeg, (max(PIDeg) - min(PIDeg) +1)/5)
    plt.loglog(XO[:-1], YO, color = 'r', label = "Outgoing links")
    plt.loglog(XI[:-1], YI, color = 'g', label = "Incoming links")
    plt.xlabel('Weight')
    plt.ylabel('Count of weight')
    plt.title('Directed weighted signed positive links')
    plt.legend()
    plt.show()
    print "Outgoing:"
    fitPowerLaw(PODeg)
    print "Incoming:"
    fitPowerLaw(PIDeg)
    
    print "Directed weighted signed negative links"
    YO, XO = np.histogram(NODeg, (max(NODeg) - min(NODeg) +1)/5)
    YI, XI = np.histogram(NIDeg, (max(NIDeg) - min(NIDeg) +1)/5)
    plt.loglog(XO[:-1], YO, color = 'r', label = "Outgoing links")
    plt.loglog(XI[:-1], YI, color = 'g', label = "Incoming links")
    plt.xlabel('Weight')
    plt.ylabel('Count of weight')
    plt.title('Directed weighted signed negative links')
    plt.legend()
    plt.show()
    print "Outgoing:"
    fitPowerLaw(NODeg)
    print "Incoming:"
    fitPowerLaw(NIDeg)


# In[97]:


G = selfEdgeDel(snap.LoadEdgeList(snap.PNGraph, "Datasets/soc-sign-bitcoinotc.csv", 0, 1, ','))
signs = loadSigns("Datasets/soc-sign-bitcoinotc.csv")
GetDirectedWeightedDegreeDistribution(G, signs)


# In[9]:


def GetDirectedUnweightedDegreeDistribution(G, signs):
    PIDeg = []
    PODeg = []
    NIDeg = []
    NODeg = []
    
    for N in G.Nodes():
        OutDegree = N.GetOutDeg()
        weightP = 0
        weightN = 0
        for i in range(OutDegree):
            weight = signs[(N.GetId(), N.GetOutNId(i))]
            if weight >= 0:
                weightP += 1
            else:
                weightN -= -1
        PODeg.append(weightP)
        NODeg.append(weightN)
        
        InDegree = N.GetInDeg()
        weightP = 0
        weightN = 0
        for i in range(InDegree):
            weight = signs[(N.GetInNId(i), N.GetId())]
            if weight >= 0:
                weightP += 1
            else:
                weightN -= -1
        PIDeg.append(weightP)
        NIDeg.append(weightN)
    
    print "Directed unweighted signed positive links"
    YO, XO = np.histogram(PODeg, (max(PODeg) - min(PODeg) +1)/5)
    YI, XI = np.histogram(PIDeg, (max(PIDeg) - min(PIDeg) +1)/5)
    plt.loglog(XO[:-1], YO, color = 'r', label = "Outgoing links")
    plt.loglog(XI[:-1], YI, color = 'g', label = "Incoming links")
    plt.xlabel('Weight')
    plt.ylabel('Count of weight')
    plt.title('Directed unweighted signed positive links')
    plt.legend()
    plt.show()
    print "Outgoing:"
    fitPowerLaw(PODeg)
    print "Incoming:"
    fitPowerLaw(PIDeg)
    
    print "Directed unweighted signed negative links"
    YO, XO = np.histogram(NODeg, (max(NODeg) - min(NODeg) +1)/5)
    YI, XI = np.histogram(NIDeg, (max(NIDeg) - min(NIDeg) +1)/5)
    plt.loglog(XO[:-1], YO, color = 'r', label = "Outgoing links")
    plt.loglog(XI[:-1], YI, color = 'g', label = "Incoming links")
    plt.xlabel('Weight')
    plt.ylabel('Count of weight')
    plt.title('Directed unweighted signed negative links')
    plt.legend()
    plt.show()
    print "Outgoing:"
    fitPowerLaw(NODeg)
    print "Incoming:"
    fitPowerLaw(NIDeg)


# In[95]:


GD = selfEdgeDel(snap.LoadEdgeList(snap.PNGraph, "Datasets/soc-sign-bitcoinotc.csv", 0, 1, ','))
signs = loadSigns("Datasets/soc-sign-bitcoinotc.csv")
GetDirectedUnweightedDegreeDistribution(G, signs)


# In[15]:


def GetAdjacencyMatrix(G, signs):
    maps = [N.GetId() for N in G.Nodes()]
    maps.sort()
    print len(maps)
        
    A = []
    AU = []
    AUU = []
    for i in maps:
        if i % 1000 == 0:
            print i,
        a = []
        au = []
        auu = []
        for j in maps:
            if G.IsEdge(i, j):
                a.append(signs[(i, j)])
                if signs[(i, j)] > 0:
                    au.append(1)
                    auu.append(1)
                else:
                    au.append(-1)
                    auu.append(1)
            else:
                a.append(0)
                au.append(0)
                auu.append(0)
        A.append(a)
        AU.append(au)
        AUU.append(auu)
    
#     return np.matrix(A), np.matrix(AU), np.matrix(AUU)
    return A, AU, AUU


# In[21]:


def SquareMatrix(A):
    n = len(A)
    b = [0 for i in range(n)]
    B = [b for i in range(n)]
    
    for i in range(n):
        print i
        for j in range(n):
            for k in range(n):
                B[i][j] += A[i][k] * A[k][j]
                
    return B


# In[44]:


def ModSign(A):
    B = {}
    for i in A.keys():
        B[i] = abs(A[i])
    return B


# In[49]:


def NormalizeSign(A):
    M = max(ModSign(A).values())
    B = {}
    for i in A.keys():
        B[i] = 1.0 * A[i] / M
    return B


# In[64]:


def UnweightedSign(A):
    B = {}
    for i in A.keys():
        if A[i] > 0:
            B[i] = 1
        else:
            B[i] = -1
    return B


# In[101]:


def CC_Undirected_Unsigned_Unweighted(filename):
    G = selfEdgeDel(snap.LoadEdgeList(snap.PUNGraph, filename, 0, 1, ','))
    signs = loadSigns(filename)
    
    C = []
    for X in G.Nodes():
        c = 0
        Deg = X.GetDeg()
        for y in range(Deg):
            Y = G.GetNI(X.GetNbrNId(y))
            YDeg = Y.GetDeg()
            for z in range(YDeg):
                ZId = Y.GetNbrNId(z)
                if ZId != X.GetId() and Y.GetId() < ZId:
                    if G.IsEdge(X.GetId(), ZId):
                        c += 1
        if Deg != 0 and Deg != 1:
            C.append(2.0 * c / Deg / (Deg - 1))
            
    print sum(C)/len(C)
    
# CC_Undirected_Unsigned_Unweighted()

# def CC_Undirected_Unsigned_Unweighted():
#     G = selfEdgeDel(snap.LoadEdgeList(snap.PUNGraph, "Datasets/soc-sign-bitcoinotc.csv", 0, 1, ','))
#     signs = loadSigns("Datasets/soc-sign-bitcoinotc.csv")
    
#     n = 0
#     d = 0
#     for X in G.Nodes():
#         Deg = X.GetDeg()
#         for y in range(Deg):
#             Y = G.GetNI(X.GetNbrNId(y))
#             for z in range(Deg):
#                 Z = G.GetNI(X.GetNbrNId(z))
#                 # Unique triangles
#                 if Y.GetId() < X.GetId() and Z.GetId() < Y.GetId():
#                     if G.IsEdge(Y.GetId(), Z.GetId()):
#                         n += 1
#                     d += 1
            
#     print 1.0 * n / d
    
# CC_Undirected_Unsigned_Unweighted()


# In[110]:


def CC_Directed_Unsigned_Unweighted(filename):
    G = selfEdgeDel(snap.LoadEdgeList(snap.PNGraph, filename, 0, 1, ','))
    signs = loadDirectedSigns(filename)
    
    C = []
    for X in G.Nodes():
        c = 0
        Deg = X.GetOutDeg()
        for y in range(Deg):
            Y = G.GetNI(X.GetOutNId(y))
            YDeg = Y.GetOutDeg()
            for z in range(YDeg):
                ZId = Y.GetOutNId(z)
                if ZId != X.GetId():
                    if G.IsEdge(X.GetId(), ZId):
                        c += 1
        if Deg != 0 and Deg != 1:
            C.append(1.0 * c / Deg / (Deg - 1))
            
    print sum(C)/len(C)
    
# CC_Directed_Unsigned_Unweighted()

# def CheckTriangleConsistency(G, X, Y, Z):
#     if G.IsEdge(X, Y) and G.IsEdge(Y, X):
#         return True
#     if G.IsEdge(Z, Y) and G.IsEdge(Y, Z):
#         return True
#     if G.IsEdge(Z, X) and G.IsEdge(X, Z):
#         return True
#     if G.IsEdge(X, Y) and G.IsEdge(Y, Z) and G.IsEdge(Z, X):
#         return False
#     if G.IsEdge(Y, X) and G.IsEdge(Z, Y) and G.IsEdge(X, Z):
#         return False
#     return True


# def CC_Directed_Unsigned_Unweighted():
#     G = selfEdgeDel(snap.LoadEdgeList(snap.PNGraph, "Datasets/soc-sign-bitcoinotc.csv", 0, 1, ','))
#     signs = loadDirectedSigns("Datasets/soc-sign-bitcoinotc.csv")
        
#     G = snap.TNGraph.New()
#     [G.AddNode(i) for i in [0,1,2]]
#     [G.AddEdge(x,y) for (x,y) in [(0,1),(1,2),(0,2)]]
#     n = 0
#     d = 0
#     for X in G.Nodes():
#         Deg = X.GetDeg()
#         for y in range(Deg):
#             Y = G.GetNI(X.GetNbrNId(y))
#             for z in range(Deg):
#                 Z = G.GetNI(X.GetNbrNId(z))
#                 # Unique triangles
#                 if Y.GetId() < X.GetId() and Z.GetId() < Y.GetId():
#                     if G.IsEdge(Y.GetId(), Z.GetId()) or G.IsEdge(Z.GetId(), Y.GetId()):
#                         if CheckTriangleConsistency(G, X.GetId(), Y.GetId(), Z.GetId()):
#                             n += 1
#                     d += 1
    
#     print 1.0 * n / d

# CC_Directed_Unsigned_Unweighted()


# In[103]:


def CC_Undirected_Unsigned_Weighted(filename):
    G = selfEdgeDel(snap.LoadEdgeList(snap.PUNGraph, filename, 0, 1, ','))
    signs = loadSigns(filename)
    signs = NormalizeSign(ModSign(signs))
    
    C = []
    for X in G.Nodes():
        n = 0
        d = 0
        Deg = X.GetDeg()
        for y in range(Deg):
            Y = G.GetNI(X.GetNbrNId(y))
            YDeg = Y.GetDeg()
            for z in range(YDeg):
                ZId = Y.GetNbrNId(z)
                if ZId != X.GetId():
                    if G.IsEdge(X.GetId(), ZId):
                        n += signs[(X.GetId(), Y.GetId())] * signs[(Y.GetId(), ZId)] * signs[(X.GetId(), ZId)]
            
            for z in range(Deg):
                ZId = X.GetNbrNId(z)
                if ZId != Y.GetId():
                    d += signs[(X.GetId(), Y.GetId())] * signs[(X.GetId(), ZId)]
        if Deg != 0 and Deg != 1:
            C.append(1.0 * n / d)
            
    print sum(C)/len(C)
    
# CC_Undirected_Unsigned_Weighted()


# In[104]:


def CC_Directed_Unsigned_Weighted(filename):
    G = selfEdgeDel(snap.LoadEdgeList(snap.PNGraph, filename, 0, 1, ','))
    signs = loadDirectedSigns(filename)
    signs = NormalizeSign(ModSign(signs))
    
    C = []
    for X in G.Nodes():
        n = 0
        d = 0
        Deg = X.GetOutDeg()
        for y in range(Deg):
            Y = G.GetNI(X.GetOutNId(y))
            YDeg = Y.GetOutDeg()
            for z in range(YDeg):
                ZId = Y.GetOutNId(z)
                if ZId != X.GetId():
                    if G.IsEdge(X.GetId(), ZId):
                        n += signs[(X.GetId(), Y.GetId())] * signs[(Y.GetId(), ZId)] * signs[(X.GetId(), ZId)]

            for z in range(Deg):
                ZId = X.GetOutNId(z)
                if ZId != Y.GetId():
                    d += signs[(X.GetId(), Y.GetId())] * signs[(X.GetId(), ZId)]
        if Deg != 0 and Deg != 1:
            C.append(1.0 * n / d)
            
    print sum(C)/len(C)
    
# CC_Directed_Unsigned_Weighted()


# In[105]:


def CC_Undirected_Signed_Unweighted(filename):
    G = selfEdgeDel(snap.LoadEdgeList(snap.PUNGraph, filename, 0, 1, ','))
    signs = loadSigns(filename)
    signs = UnweightedSign(signs)
    
    C = []
    for X in G.Nodes():
        n = 0
        d = 0
        Deg = X.GetDeg()
        for y in range(Deg):
            Y = G.GetNI(X.GetNbrNId(y))
            YDeg = Y.GetDeg()
            for z in range(YDeg):
                ZId = Y.GetNbrNId(z)
                if ZId != X.GetId() and ZId < Y.GetId():
                    if G.IsEdge(X.GetId(), ZId):
                        if signs[(X.GetId(), Y.GetId())] * signs[(Y.GetId(), ZId)] == signs[(X.GetId(), ZId)]:
                            n += 1
            
            d = Deg * (Deg - 1) / 2
        if Deg != 0 and Deg != 1:
            C.append(1.0 * n / d)
            
    print sum(C)/len(C)
    
# CC_Undirected_Signed_Unweighted()


# In[106]:


def CC_Directed_Signed_Unweighted(filename):
    G = selfEdgeDel(snap.LoadEdgeList(snap.PNGraph, filename, 0, 1, ','))
    signs = loadDirectedSigns(filename)
    normalizedsigns = NormalizeSign(ModSign(signs))
    signs = UnweightedSign(signs)
    
    C = []
    for X in G.Nodes():
        n = 0
        d = 0
        Deg = X.GetOutDeg()
        for y in range(Deg):
            Y = G.GetNI(X.GetOutNId(y))
            YDeg = Y.GetOutDeg()
            for z in range(YDeg):
                ZId = Y.GetOutNId(z)
                if ZId != X.GetId():
                    if G.IsEdge(X.GetId(), ZId):
                        if signs[(X.GetId(), Y.GetId())] * signs[(Y.GetId(), ZId)] == signs[(X.GetId(), ZId)]:
                            n += 1
            
            d = Deg * (Deg - 1)
        if Deg != 0 and Deg != 1:
            C.append(1.0 * n / d)
            
    print sum(C)/len(C)
    
# CC_Directed_Signed_Unweighted()


# In[107]:


def CC_Undirected_Signed_Weighted(filename):
    G = selfEdgeDel(snap.LoadEdgeList(snap.PUNGraph, filename, 0, 1, ','))
    signs = loadSigns(filename)
    normalizedSigns = NormalizeSign(ModSign(signs))
    signs = UnweightedSign(signs)
    
    C = []
    for X in G.Nodes():
        n = 0
        d = 0
        Deg = X.GetDeg()
        for y in range(Deg):
            Y = G.GetNI(X.GetNbrNId(y))
            YDeg = Y.GetDeg()
            for z in range(YDeg):
                ZId = Y.GetNbrNId(z)
                if ZId != X.GetId() and ZId < Y.GetId():
                    if G.IsEdge(X.GetId(), ZId):
                        if signs[(X.GetId(), Y.GetId())] * signs[(Y.GetId(), ZId)] == signs[(X.GetId(), ZId)]:
                            n += normalizedSigns[(X.GetId(), Y.GetId())] * normalizedSigns[(Y.GetId(), ZId)] * normalizedSigns[(X.GetId(), ZId)]
            
            for z in range(Deg):
                ZId = X.GetNbrNId(z)
                if ZId != Y.GetId() and ZId < Y.GetId():
                    d += normalizedSigns[(X.GetId(), Y.GetId())] * normalizedSigns[(X.GetId(), ZId)]
        if Deg != 0 and Deg != 1:
            C.append(1.0 * n / d)
            
    print sum(C)/len(C)
    
# CC_Undirected_Signed_Weighted()


# In[108]:


def CC_Directed_Signed_Weighted(filename):
    G = selfEdgeDel(snap.LoadEdgeList(snap.PNGraph, filename, 0, 1, ','))
    signs = loadDirectedSigns(filename)
    normalizedSigns = NormalizeSign(ModSign(signs))
    signs = UnweightedSign(signs)
    
    C = []
    for X in G.Nodes():
        n = 0
        d = 0
        Deg = X.GetOutDeg()
        for y in range(Deg):
            Y = G.GetNI(X.GetOutNId(y))
            YDeg = Y.GetOutDeg()
            for z in range(YDeg):
                ZId = Y.GetOutNId(z)
                if ZId != X.GetId():
                    if G.IsEdge(X.GetId(), ZId):
                        if signs[(X.GetId(), Y.GetId())] * signs[(Y.GetId(), ZId)] == signs[(X.GetId(), ZId)]:
                            n += normalizedSigns[(X.GetId(), Y.GetId())] * normalizedSigns[(Y.GetId(), ZId)] * normalizedSigns[(X.GetId(), ZId)]
            
            for z in range(Deg):
                ZId = X.GetNbrNId(z)
                if ZId != Y.GetId():
                    d += normalizedSigns[(X.GetId(), Y.GetId())] * normalizedSigns[(X.GetId(), ZId)]
        if Deg != 0 and Deg != 1:
            C.append(1.0 * n / d)
            
    print sum(C)/len(C)
    
# CC_Directed_Signed_Weighted()


# In[111]:


for filename in ["Datasets/soc-sign-bitcoinotc.csv", "Datasets/soc-sign-bitcoinalpha.csv"]:
    print filename
    print "CC_Undirected_Unsigned_Unweighted",
    CC_Undirected_Unsigned_Unweighted(filename)
    print "CC_Directed_Unsigned_Unweighted",
    CC_Directed_Unsigned_Unweighted(filename)
    print "CC_Undirected_Unsigned_Weighted",
    CC_Undirected_Unsigned_Weighted(filename)
    print "CC_Directed_Unsigned_Weighted",
    CC_Directed_Unsigned_Weighted(filename)
    print "CC_Undirected_Signed_Unweighted",
    CC_Undirected_Signed_Unweighted(filename)
    print "CC_Directed_Signed_Unweighted",
    CC_Directed_Signed_Unweighted(filename)
    print "CC_Undirected_Signed_Weighted",
    CC_Undirected_Signed_Weighted(filename)
    print "CC_Directed_Signed_Weighted",
    CC_Directed_Signed_Weighted(filename)


# In[22]:


GD = selfEdgeDel(snap.LoadEdgeList(snap.PNGraph, "Datasets/soc-sign-bitcoinotc.csv", 0, 1, ','))
signs = loadSigns("Datasets/soc-sign-bitcoinotc.csv")
A, AU, AUU = GetAdjacencyMatrix(GD, signs)


# In[23]:


print GD.GetNodes(), GD.GetEdges()


# In[27]:


G = snap.TNGraph.New()
[G.AddNode(i) for i in [0,1,2]]
[G.AddEdge(i,j) for (i,j) in [(0,1),(1,2),(2,0)]]

for N in G.Nodes():
    OutDeg = N.GetOutDeg()
    for i in range(OutDeg):
        print str(N.GetId()) + " -> " + str(N.GetOutNId(i))
    InDeg = N.GetInDeg()
    for i in range(InDeg):
        print str(N.GetId()) + " <- " + str(N.GetInNId(i))


# In[51]:


d = {}
d[(1)] = -4
d[(2)] = -8
print ModSign(d), NormalizeSign(d), ModSign(NormalizeSign(d))


# In[ ]:


# G = selfEdgeDel(snap.LoadEdgeList(snap.PNGraph, "Datasets/soc-sign-bitcoinotc.csv", 0, 1, ','))
# signs = loadDirectedSigns("Datasets/soc-sign-bitcoinotc.csv")

# R = []
# for N in G.Nodes():
#     for i in range(N.GetOutDeg()):
#         M = G.GetNI(N.GetOutNId(i))
        

