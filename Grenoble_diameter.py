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



clear = lambda: os.system('cls')
clear()



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
    G = nx.Graph()
    m = max(list(weight.values()))
    for i in list(weight.keys()):
        if (weight[i] != 0):
            G.add_edge(i[0], i[1], weight=m)
    return G


#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
def subgraph_max_weight(G):
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
    Gs = nx.Graph()
    m = max(list(weight.values()))
    for i in list(weight.keys()):
        if (weight[i] == m):
            Gs.add_edge(i[0], i[1], weight=m)
    return Gs



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

#-------------------------------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------------------------------------
def duplicate_d(dw):
    iterator_k = list(dw.keys())
    ddw = {}
    for i in iterator_k:
        ddw[(i[0], i[1])] = dw[i]
        ddw[(i[1], i[0])] = dw[i]
    return ddw

#-------------------------------------------------------------------------------------------------------------------------------------------------------------
def edge_clustering(G, dd, e):
    '''
    ----------
     Parameters
     ----------
     G: networkx weighted graph.

     dw: dictionary
         links as the keys and weights as the values
     e: list of two elements
        edge of the graph G
     --------
     Returns:
     --------
     nce/dce: float
         edge clustering for weighted graphs
    '''
    Nu = set()
    Nv = set()
    for i in G.neighbors(e[0]):
        Nu.add(i)
    for i in G.neighbors(e[1]):
        Nv.add(i)
    iN = Nu.intersection(Nv)
    uN = Nu.union(Nv)
    dd = degree_G(G, dw)
    dw = duplicate_d(dd)
    nce = 0
    for i in iN:
        nce = nce + dw[(i,e[0])] + dw[(i,e[1])]
    dce = 0
    for i in uN:
        dce = dce + dw[(i, e[0])] + dw[(i, e[1])]
    dce = dce - 2*ddw[e]
    return nce/dce



#-------------------------------------------------------------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------------------------------------------------------
def edge_clustering_coefficient(G, dw):
    de = {}
    for e in list(G.edges):
        sui = set()
        svi = set()
        suu = set()
        svu = set()
        dw = duplicate_d(dw)
        for i in list(G[e[0]]):
            suu.add(i)
            if (dw[(e[0],i)] >= 3):
                sui.add(i)
        for i in list(G[e[1]]):
            svu.add(i)
            if (dw[(e[1], i)] >= 2):
                svi.add(i)
        #if (weight[e] <=2):
        #    ph_i = 1
        #else:
            ph_i = 0
        nume = len(list(sui.intersection(svi)))
        denom = len(list(suu.union(svu)))
        if (denom - 2> 0):
            de[e] = nume/(denom - 2)
        else:
            de[e] = 0
    return de
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
def diameter_plot(CG, SsG, SG, pos):
    # edge weight labels
    actual = diameter_connected_c(CG, SsG)
    safe = diameter_connected_c(SG, SG)
    #largest_cc = max(nx.connected_components(SsG), key=len)
    #color_map = []
    #for node in G:
    #    if node in largest_cc:
    #        color_map.append('red')
    #    else:
    #        color_map.append('black')
    plt.figure(figsize=(8, 10))
    #s1 = [(u, v) for (u, v, d) in CG.edges(data=True) if d["weight"] == 1]
    #s2 = [(u, v) for (u, v, d) in CG.edges(data=True) if d["weight"] == 2]
    #s3 = [(u, v) for (u, v, d) in CG.edges(data=True) if d["weight"] == 3]
    #s4 = [(u, v) for (u, v, d) in CG.edges(data=True) if d["weight"] == 4]

    #nx.draw_networkx_edges(CG, pos, edgelist=s1, width=1, edge_color="r")
    #nx.draw_networkx_edges(CG, pos, edgelist=s2, width=1, edge_color="tab:orange")
    #nx.draw_networkx_edges(CG, pos, edgelist=s3, width=1, edge_color="y")
    #nx.draw_networkx_edges(CG, pos, edgelist=s4, width=2, edge_color="c")
    hsv_modified = cm.get_cmap('nipy_spectral_r', 256)  # create new hsv colormaps in range of 0.3 (green) to 0.7 (blue)
    newcmp = ListedColormap(hsv_modified(np.linspace(0.05, 0.50, 256)))  # show figure
    cmap = newcmp
    # node labels
    # nx.draw_networkx_edge_labels(G, pos, font_size=1, font_family="sans-serif", alpha=0.7)
    e_attributes = nx.get_edge_attributes(CG, 'weight')
    #nx.draw_networkx_edge_labels(CG, pos, edge_labels=edge_labels, font_size=7, font_family="sans-serif", alpha=1)
    nx.draw_networkx_nodes(
        CG,
        pos,
        nodelist=CG.nodes(),
        node_size=0,
        #node_color=color_map
        # cmap=plt.cm.Blues,
        # cmap=plt.cm.Reds_r,
        # cmap=cmap,
    )
    edges = nx.draw_networkx_edges(CG, pos, edge_color=list(e_attributes.values()), width=1,
                                   edge_cmap=cmap)
    plt.colorbar(edges)
    ax = plt.gca()
    ax.set_axis_off()
    T = (actual/safe)*100
    #plt.xlim(-5.7, 3.8)
    #plt.ylim(-1.8, 3.9)
    plt.axis("off")
    plt.title('Safest diameter = ' + str(T)+'%')
    plt.savefig('diameter_metric.eps', format='eps')
    plt.show()
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
def diameter_connected_c(G, SG):
    largest_cc = max(nx.connected_components(SG), key=len)
    H = G.subgraph(list(largest_cc))
    return nx.diameter(H) + len(largest_cc)


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




CG = current_graph(G)
SsG = subgraph_max_weight(G)
SG = max_weight(G)

diameter_plot(CG, SsG, SG, pos)
