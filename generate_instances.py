import networkx as nx
import matplotlib.pyplot as plt
import math
import csv

def print_matrix(filename, pos, metric):
    with open('./instances/'+filename, 'wb') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        for _i1, c in pos.iteritems():
            line = []
            for _i2, c2 in pos.iteritems():
                line.append(metric(c, c2)*100)
            writer.writerow(line)


def euclidean_distance((a, b), (c, d)):
    return math.sqrt((a - c)**2 + (b - d)**2)

def manhattan_distance((a, b), (c, d)):
    return abs(a-c) + abs(b-d)


dims = [9, 16, 25, 36, 49, 64]

#
# Generate random layouts
#

for dim in dims:
    for i in range(0, 32):
        G = nx.empty_graph(dim)
        pos = nx.random_layout(G)
        print_matrix("rnd_"+str(dim)+"_ed_"+str(i), pos, euclidean_distance)
        print_matrix("rnd_"+str(dim)+"_md_"+str(i), pos, manhattan_distance)


#
# Generate uniform layouts
#

for dim in dims:
    G = nx.grid_2d_graph(int(math.sqrt(dim)), int(math.sqrt(dim)))
    pos = nx.graphviz_layout(G, prog='neato')
    print_matrix("uni_"+str(dim)+"_ed", pos, euclidean_distance)
    print_matrix("uni_"+str(dim)+"_md", pos, manhattan_distance)


#
# Generate semi-clustered layouts
#

for dim in dims:
    for i in range(0,32):
        G = nx.barabasi_albert_graph(dim, int(math.sqrt(dim)))
        pos = nx.graphviz_layout(G, prog='dot')
        print_matrix("cls_"+str(dim)+"_ed_"+str(i), pos, euclidean_distance)
        print_matrix("cls_"+str(dim)+"_md_"+str(i), pos, manhattan_distance)

