import snap
import random
import numpy as np
import matplotlib.pyplot as plt
import powerlaw
import math
import networkx as nx
import Queue


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


def GetNetworkxG(G, signs):
    Gnx = nx.DiGraph()
    Gnxu = nx.Graph()
    for N in G.Nodes():
        Gnx.add_node(N.GetId())
        Gnxu.add_node(N.GetId())
    for E in G.Edges():
        src = E.GetSrcNId()
        dst = E.GetDstNId()
        Gnx.add_weighted_edges_from([(src, dst, signs[(src, dst)])])
        Gnxu.add_weighted_edges_from([(src, dst, signs[(src, dst)])])

    return Gnx, Gnxu


def makeSignsPositive(signs):
    minSign = min(signs.values())
    modMinSign = abs(minSign)
    newSigns = {}
    for i in signs.keys():
        newSigns[i] = signs[i] + modMinSign + 1
    return newSigns


def GetCNCoeff(G, signs, x, y):
    # signs > 0
    # x -> y

    xo = GetNeighO(G, signs, x)
    xi = GetNeighI(G, signs, x)
    yo = GetNeighO(G, signs, y)
    yi = GetNeighI(G, signs, y)

    # wxy
    woo = sum([xo[i] if i in yo.keys() else 0 for i in xo.keys()])
    wio = sum([xi[i] if i in yo.keys() else 0 for i in xi.keys()])
    woi = sum([xo[i] if i in yi.keys() else 0 for i in xo.keys()])
    wii = sum([xi[i] if i in yi.keys() else 0 for i in xi.keys()])

    return woo, wio, woi, wii


def GetJaccardCoeff(G, signs, x, y):
    # signs > 0
    # x -> y

    xo = GetNeighO(G, signs, x)
    xi = GetNeighI(G, signs, x)
    yo = GetNeighO(G, signs, y)
    yi = GetNeighI(G, signs, y)

    # nxy
    noo = sum([xo[i] if i in yo.keys() else 0 for i in xo.keys()])
    nio = sum([xi[i] if i in yo.keys() else 0 for i in xi.keys()])
    noi = sum([xo[i] if i in yi.keys() else 0 for i in xo.keys()])
    nii = sum([xi[i] if i in yi.keys() else 0 for i in xi.keys()])

    # dxy
    doo = sum(xo.values()) + sum(yo.values())
    dio = sum(xi.values()) + sum(yo.values())
    doi = sum(xo.values()) + sum(yi.values())
    dii = sum(xi.values()) + sum(yi.values())

    foo = 1.0 * noo / doo if doo != 0 else 0
    fio = 1.0 * nio / dio if dio != 0 else 0
    foi = 1.0 * noi / doi if doi != 0 else 0
    fii = 1.0 * nii / dii if dii != 0 else 0

    return foo, fio, foi, fii


def PACoeff(G, signs, x, y):
    # signs > 0
    # x -> y

    xo = GetNeighO(G, signs, x)
    xi = GetNeighI(G, signs, x)
    yo = GetNeighO(G, signs, y)
    yi = GetNeighI(G, signs, y)

    woo = sum(xo.values()) * sum(yo.values())
    wio = sum(xi.values()) * sum(yo.values())
    woi = sum(xo.values()) * sum(yi.values())
    wii = sum(xi.values()) * sum(yi.values())

    return woo, wio, woi, wii


def AACoeff(G, signs, x, y, VolO, VolI, Vol):
    xo = GetNeighO(G, signs, x)
    xi = GetNeighI(G, signs, x)
    yo = GetNeighO(G, signs, y)
    yi = GetNeighI(G, signs, y)

    woo = sum([1.0 * (xo[z] + yo[z]) / math.log(1 + VolI[z]) if z in yo.keys() else 0 for z in xo.keys()])
    wio = sum([1.0 * (xi[z] + yo[z]) / math.log(1 + Vol[z]) if z in yo.keys() else 0 for z in xi.keys()])
    woi = sum([1.0 * (xo[z] + yi[z]) / math.log(1 + Vol[z]) if z in yi.keys() else 0 for z in xo.keys()])
    wii = sum([1.0 * (xi[z] + yi[z]) / math.log(1 + VolO[z]) if z in yi.keys() else 0 for z in xi.keys()])

    return woo, wio, woi, wii


def RACoeff(G, signs, x, y, VolO, VolI, Vol):
    xo = GetNeighO(G, signs, x)
    xi = GetNeighI(G, signs, x)
    yo = GetNeighO(G, signs, y)
    yi = GetNeighI(G, signs, y)

    woo = sum([1.0 * (xo[z] + yo[z]) / VolI[z] if z in yo.keys() else 0 for z in xo.keys()])
    wio = sum([1.0 * (xi[z] + yo[z]) / Vol[z] if z in yo.keys() else 0 for z in xi.keys()])
    woi = sum([1.0 * (xo[z] + yi[z]) / Vol[z] if z in yi.keys() else 0 for z in xo.keys()])
    wii = sum([1.0 * (xi[z] + yi[z]) / VolO[z] if z in yi.keys() else 0 for z in xi.keys()])

    return woo, wio, woi, wii


def CCCoeff(G, signs, x, y, C):
    return C[x] + C[y]


def GetNeighO(G, signs, NId):
    N = G.GetNI(NId)
    Neighs = {}
    OutDeg = N.GetOutDeg()
    for i in range(OutDeg):
        Neigh = N.GetOutNId(i)
        Weight = signs[(N.GetId(), Neigh)]
        Neighs[Neigh] = Weight
    return Neighs


def GetNeighI(G, signs, NId):
    N = G.GetNI(NId)
    Neighs = {}
    InDeg = N.GetInDeg()
    for i in range(InDeg):
        Neigh = N.GetInNId(i)
        Weight = signs[(Neigh, N.GetId())]
        Neighs[Neigh] = Weight
    return Neighs


def GetVolO(G, signs, NId):
    return sum(GetNeighO(G, signs, NId).values())


def GetVolI(G, signs, NId):
    return sum(GetNeighI(G, signs, NId).values())


def GetVol(G, signs, NId):
    return GetVolO(G, signs, NId) + GetVolI(G, signs, NId)


def GetVectors(G, signs):
    # G = snap.LoadEdgeList(snap.PNGraph, "Datasets/soc-sign-bitcoinotc.csv", 0, 1, ',')
    # signs = loadDirectedSigns("Datasets/soc-sign-bitcoinotc.csv")
    pSigns = makeSignsPositive(signs)
    Gnx, Gnxu = GetNetworkxG(G, signs)

    c = 0
    Vol = {}
    VolI = {}
    VolO = {}
    CC = {}
    for N in G.Nodes():
        CC[N.GetId()] = nx.clustering(Gnxu, N.GetId())
        VolO[N.GetId()] = GetVolO(G, pSigns, N.GetId())
        VolI[N.GetId()] = GetVolI(G, pSigns, N.GetId())
        Vol[N.GetId()] = GetVol(G, pSigns, N.GetId())

    M = []
    Y = []
    for e in G.Edges():
        if c % 100 == 0:
            print c,
        f = []
        src = e.GetSrcNId()
        dst = e.GetDstNId()

        c1, c2, c3, c4 = GetCNCoeff(G, pSigns, src, dst)
        f.append(c1)
        f.append(c2)
        f.append(c3)
        f.append(c4)

        c1, c2, c3, c4 = GetJaccardCoeff(G, pSigns, src, dst)
        f.append(c1)
        f.append(c2)
        f.append(c3)
        f.append(c4)

        c1, c2, c3, c4 = PACoeff(G, pSigns, src, dst)
        f.append(c1)
        f.append(c2)
        f.append(c3)
        f.append(c4)

        c1, c2, c3, c4 = AACoeff(G, pSigns, src, dst, VolO, VolI, Vol)
        f.append(c1)
        f.append(c2)
        f.append(c3)
        f.append(c4)

        c1, c2, c3, c4 = RACoeff(G, pSigns, src, dst, VolO, VolI, Vol)
        f.append(c1)
        f.append(c2)
        f.append(c3)
        f.append(c4)

        c1 = CCCoeff(G, pSigns, src, dst, CC)
        f.append(c1)

        f.append(1)

        M.append(f)
        Y.append(pSigns[src, dst])

        c += 1

    return np.array(M), np.array(Y)

def RestrictGraph(G, x, y):
    depths = {}
    Q = Queue.Queue()
    Q.put(x)
    Q.put(y)
    depths[x] = 0
    depths[y] = 0

    while not Q.empty():
        NId = Q.get()
        N = G.GetNI(NId)
        Deg = N.GetDeg()
        for i in range(Deg):
            Neigh = N.GetNbrNId(i)
            if Neigh not in depths.keys():
                Q.put(Neigh)
                depths[Neigh] = depths[NId] + 1

    print len(depths.keys()), G.GetNodes()
    MaxDepth = 2
    Gr = snap.TNGraph.New()
    for N in G.Nodes():
        if (N.GetId() in depths.keys()) and (depths[N.GetId()] <= MaxDepth):
            Gr.AddNode(N.GetId())
    for E in G.Edges():
        if Gr.IsNode(E.GetSrcNId()) and Gr.IsNode(E.GetDstNId()):
            Gr.AddEdge(E.GetSrcNId(), E.GetDstNId())

    return Gr

def GetRestrictedSigns(Gr, signs):
    rSigns = {}
    for (x, y) in signs.keys():
        if Gr.IsEdge(x, y):
            rSigns[(x, y)] = signs[(x, y)]
    return rSigns

G = snap.LoadEdgeList(snap.PNGraph, "Datasets/soc-sign-bitcoinotc.csv", 0, 1, ',')
signs = loadDirectedSigns("Datasets/soc-sign-bitcoinotc.csv")
(x,y) = random.choice([(E.GetSrcNId(), E.GetDstNId()) for E in G.Edges()])
Gr = RestrictGraph(G, x, y)
rSigns = GetRestrictedSigns(Gr, signs)
M, Y = GetVectors(G, rSigns)
print M.shape, Y.shape
