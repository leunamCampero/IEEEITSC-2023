import geopandas as gpd
import matplotlib.pyplot as plt
import momepy
import networkx as nx
from contextily import add_basemap
from libpysal import weights
from shapely.geometry import LineString
import pandas as pd
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import matplotlib as mpl
import os
from itertools import combinations
from scipy.sparse import csr_matrix
from scipy import linalg
from collections import Counter


clear = lambda: os.system('cls')
clear()


def inverse_weight(G):
    '''
    ----------
     Parameters
     ----------
     dw: dictionary
         links as the keys and weights as the values

     --------
     Returns:
     --------
     G: networkx weighted graph.
     dw: dictionary
         Keys are the edges and the values is safety road classification of the corresponding edge.

    '''
    weight = nx.get_edge_attributes(G, 'weight')
    Gi = G.copy()
    m = max(list(weight.values()))
    for u, v, d in Gi.edges(data=True):
        d['weight'] = m+1-d['weight']
    return Gi

def max_weight(G):
    '''
    ----------
     Parameters
     ----------
     dw: dictionary
         links as the keys and weights as the values

     --------
     Returns:
     --------
     G: networkx weighted graph.
     dw: dictionary
         Keys are the edges and the values is safety road classification of the corresponding edge.

    '''
    weight = nx.get_edge_attributes(G, 'weight')
    Gw = G.copy()
    m = max(list(weight.values()))
    for u, v, d in Gw.edges(data=True):
        d['weight'] = m
    return Gw


#-------------------------------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------------------------------------
def current_graph(G):
    '''
    ----------
     Parameters
     ----------
     dw: dictionary
         links as the keys and weights as the values

     --------
     Returns:
     --------
     G: networkx weighted graph.
     dw: dictionary
         Keys are the edges and the values is safety road classification of the corresponding edge.

    '''
    weight = nx.get_edge_attributes(G, 'weight')
    Gw = nx.Graph()
    for i in list(weight.keys()):
        if (weight[i] != 0):
            Gw.add_edge(i[0], i[1], weight = weight[i])
    return Gw
#-------------------------------------------------------------------------------------------------------------------------------------------------------------

def duplicate_d(dw):
    iterator_k = list(dw.keys())
    s = set()
    for i in list(iterator_k):
        s.add(i[0])
        s.add(i[1])
    ddw = {}
    for i in list(s):
        for j in list(s):
            if ((i,j) in iterator_k):
                ddw[(i,j)] = dw[(i,j)]
            if ((j, i) in iterator_k):
                ddw[(i, j)] = dw[(j, i)]
        ddw[i,i] = 0
    return ddw
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
def degree_G(G):
    '''
    ----------
     Parameters
     ----------
     G: networkx weighted graph.

     dw: dictionary
         links as the keys and weights as the values
     --------
     Returns:
     --------
     dd: dictionary
         keys as the nodes, and degrees of the respectively nodes as values
    '''
    H = G.degree(weight='weight')
    dd = {}
    for i in list(H):
        dd[i[0]] = i[1]
    return dd
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
def centrality_degree_plot(G, IG):
    # edge weight labels
    actual = degree_G(G)
    safe = degree_G(IG)
    p = {}
    s_actual = sum(list(actual.values()))
    s_safe = sum(list(safe.values()))
    for i in list(actual.keys()):
        p[i] = (actual[i]/safe[i])*100
    T = (s_actual/s_safe)*100
    return T
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
def sub_list(sov):
    su = []
    for i in range(20):
        sua = []
        if (i != 19):
            for j in range((i+1)*(int(len(sov)/20))+1):
                sua.append(sov[j])
        else:
            sua = sov
        su.append(sua)
    return su
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
def sub_graph(G, su):
    Gu = []
    for isu in su:
        Ga = G.copy()
        for u, v, d in Ga.edges(data=True):
            if ((u,v) in isu):
                d['weight'] = 4
        Gu.append(Ga)
    return Gu
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#df = gpd.GeoDataFrame(['a', 'b', 'c', 'd', 'e'], geometry=[l1, l2, l3, l4, l5])
bikes = gpd.read_file(r'C:\Users\ali_saleme\Documents\Phd_Inria\PyCharm\venv\velo_mobilites_m.geojson')
bikes = momepy.extend_lines(bikes, 0.00001)
bikes.geometry = momepy.close_gaps(bikes, 0.00001)
df1 = pd.DataFrame(bikes)
list(bikes['geometry'][0].coords)
#bikes.plot(figsize=(10, 10)).set_axis_off()
#bikes = momepy.extend_lines(bikes, 0.001)
#bikes_e = momepy.extend_lines(bikes, 0.0000)
#bikes_extended.plot(figsize=(10, 10)).set_axis_off()
#bikes_e.geometry = momepy.close_gaps(bikes_e, 0.000)
#bikes_e = momepy.remove_false_nodes(bikes)
#bikes = momepy.extend_lines(bikes, 1)
def coords(geom):
    return list(geom.coords)
coords = bikes.apply(lambda row: coords(row.geometry), axis=1)
coordn = coords.to_numpy()
pandadata = pd.DataFrame(bikes)
numpydata = pandadata.to_numpy()

final_points = []
G = nx.Graph()
c = 0
dw = {}
for i in coordn:
    final_points.append(i[0])
    final_points.append(i[len(i)-1])
    for j in range(len(i)-1):
        if (numpydata[c][4] == 'chronovelo'):
            G.add_edge(i[j],i[j+1], weight = 4)
            dw[(i[j],i[j+1])] = 4
            dw[(i[j + 1], i[j])] = 4
        if (numpydata[c][4] == 'veloseparatif'):
            G.add_edge(i[j],i[j+1], weight = 3)
            dw[(i[j], i[j + 1])] = 3
            dw[(i[j + 1], i[j])] = 3
        if (numpydata[c][4] == 'veloconseille'):
            G.add_edge(i[j],i[j+1], weight = 2)
            dw[(i[j], i[j + 1])] = 2
            dw[(i[j + 1], i[j])] = 2
        if (numpydata[c][4] == 'velodifficile'):
            G.add_edge(i[j],i[j+1], weight = 1)
            dw[(i[j], i[j + 1])] = 1
            dw[(i[j + 1], i[j])] = 1
    c = c + 1
# dw = duplicate_d(dw)
for i in list(G.nodes):
    nh = [n for n in G.neighbors(i)]
    if ((len(nh) == 2) and (i not in final_points)):
        G.add_edge(nh[0], nh[1], weight = dw[(i,nh[0])])
        G.remove_node(i)
for e in list(G.edges):
    if (e[0] == e[1]):
        G.remove_edge(*e)
pos = {n: [n[0], n[1]] for n in list(G.nodes)}


CG = G.copy()
IG = inverse_weight(G)

bca = nx.edge_betweenness_centrality(IG, k=None, normalized=True, weight='weight', seed=None)
vk = list(bca.keys())
s3 = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] == 3]
ind_f = [vk.index(i) for i in s3]
bci = [vk[i] for i in ind_f]
ov = {}
for i in bci:
    ov[i] = bca[i]
rsov = dict(sorted(ov.items(), key=lambda item: item[1]))
sov = list(reversed(list(rsov.keys())))
#sov = list(rsov.keys())
su = sub_list(sov)
Gu = sub_graph(G, su)
strength_vector = []
for Gi in Gu:
    CG = Gi.copy()
    MG = max_weight(Gi)
    strength_vector.append(centrality_degree_plot(CG, MG))
print(strength_vector)

f = open("Second_scenary_strengthtxt.txt", "w")
for d in strength_vector:
    f.write(f"{d}\n")
f.close()

