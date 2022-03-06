import random
import csv
import argparse
import networkx as nx
import matplotlib.pyplot as plt
import math

def generate_random_fully_connected_graph(num_nodes = 5, objective_bounds = [[50,200],[50,200]], outdir="./input/", suffix=""):
    edge_list = []

    G=nx.Graph()
    for i in range(num_nodes):
        G.add_node(i)

    for i in range(num_nodes):
        for j in range(num_nodes):
            if j > i:
                cost_ = []
                for e in objective_bounds:
                    cost_.append(random.randint(e[0],e[1]))
                ch1 = i
                ch2 = j
                edge_list.append([ch1, ch2, *cost_])
                edge_list.append([ch2, ch1, *cost_])
                G.add_edge(ch1, ch2, weight=cost_, dist=sum(cost_))
 
    pos = nx.spring_layout(G, weight="dist", k=5/math.sqrt(G.order()))
    if num_nodes <=5:
        nx.draw(G, pos, node_color='green', with_labels=True)
    else:
        nx.draw(G, pos, node_color='red', with_labels=True)
    labels = nx.get_edge_attributes(G,'weight')
    nx.draw_networkx_edge_labels(G,pos,edge_labels=labels,font_size=8)
    plt.savefig("./input/tsp_final"+suffix+".png")
    plt.clf()

    return edge_list

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=int, default=5, help="Number of nodes in generated graph")
    args = parser.parse_args()

    outdir = "./input/"

    for i in range(10):
        filename = "tsp_final_edges_5_"+str(i)+".csv"
        edges = generate_random_fully_connected_graph(num_nodes=5, outdir=outdir, suffix="_5_"+str(i))

        with open(outdir+filename, 'w', newline='') as f:
            write = csv.writer(f, delimiter=',')
            write.writerows(edges)

    for i in range(10):
        filename = "tsp_final_edges_10_"+str(i)+".csv"
        edges = generate_random_fully_connected_graph(num_nodes=10, outdir=outdir, suffix="_10_"+str(i))

        with open(outdir+filename, 'w', newline='') as f:
            write = csv.writer(f, delimiter=',')
            write.writerows(edges)

if __name__ == '__main__':
    main()