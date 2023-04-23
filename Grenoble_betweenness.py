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
    Gi = nx.Graph()
    m = max(list(weight.values()))
    for i in list(weight.keys()):
        if (weight[i] != 0):
            Gi.add_edge(i[0], i[1], weight = m+1-weight[i])
    return Gi


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
def plot_edge_betweenness_centrality(CG, IG, pos):
    bca = nx.edge_betweenness_centrality(IG, k=None, normalized=True, weight='weight', seed=None)
    bc = {}
    s_bca = sum(list(bca.values()))
    for i in list(bca.keys()):
        bc[i] = (bca[i]/s_bca)*100
    plt.figure(figsize=(12, 8))
    nx.draw_networkx_edges(CG, pos, width=1, edge_color="k")
    # node labels
    # nx.draw_networkx_edge_labels(G, pos, font_size=1, font_family="sans-serif", alpha=0.7)
    # edge_labels = nx.get_edge_attributes(Gw, 'weight')
    # nx.draw_networkx_edge_labels(Gw, pos, edge_labels=edge_labels, font_size=7, font_family="sans-serif", alpha=1)
    hsv_modified = cm.get_cmap('nipy_spectral_r', 256)  # create new hsv colormaps in range of 0.3 (green) to 0.7 (blue)
    newcmp = ListedColormap(hsv_modified(np.linspace(0.05, 0.50, 256)))  # show figure
    cmap = newcmp
    nx.draw_networkx_nodes(
        CG,
        pos,
        nodelist=G.nodes(),
        node_size=0,
        node_color="k"
        # cmap=plt.cm.Blues,
        # cmap=plt.cm.Reds_r,
        #cmap=cmap,
    )
    edges = nx.draw_networkx_edges(CG, pos, edge_color=list(bc.values()), width=1,
                                   edge_cmap=cmap)
    plt.colorbar(edges)
    #plt.xlim(-5.7, 3.8)
    #plt.ylim(-1.8, 3.9)
    plt.axis("off")
    # ax = plt.gca()
    # ax.set_title(sum(list(cc.values())))
    plt.savefig('edge_betweenness.eps', format='eps')
    plt.show()
    #plt.show()
    #pc = mpl.collections.PatchCollection(edges, cmap=cmap)
    #pc.set_array(bc)

    #ax = plt.gca()
    #ax.set_axis_off()
    #plt.colorbar(pc, ax=ax)
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
            G.add_edge(i[j],i[j+1], weight = 3)
            dw[(i[j], i[j + 1])] = 3
            dw[(i[j + 1], i[j])] = 3
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
IG = inverse_weight(G)

plot_edge_betweenness_centrality(CG, IG, pos)
