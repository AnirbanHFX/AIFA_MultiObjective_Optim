import random
import csv
import argparse
import networkx as nx
import matplotlib.pyplot as plt

def generate_random_fully_connected_graph(num_nodes = 5, objective_bounds = [[5,100],[5,100]]):
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
                G.add_edge(ch1, ch2, weight=cost_)
 
    pos = nx.circular_layout(G)
    nx.draw(G, pos)
    labels = nx.get_edge_attributes(G,'weight')
    nx.draw_networkx_edge_labels(G,pos,edge_labels=labels)
    plt.savefig("./input/tsp_final.png")

    return edge_list

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=int, default=5, help="Number of nodes in generated graph")
    args = parser.parse_args()

    outdir = "./input/"
    filename = "tsp_final_edges.csv"
    edges = generate_random_fully_connected_graph(num_nodes=args.n)

    with open(outdir+filename, 'w', newline='') as f:
        write = csv.writer(f, delimiter=',')
        write.writerows(edges)

if __name__ == '__main__':
    main()