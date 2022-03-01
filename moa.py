import csv

class graph:

    def __init__(self, vertex_path=None, edge_path=None):
        self._adjacency_list = {}
        self._vertex_list = {}
        if (vertex_path is not None) and (edge_path is not None):
            self.read_graph(vertex_path, edge_path)

    @property
    def adjacency_list(self):
        return self.adjacency_list

    @adjacency_list.setter
    def adjacency_list(self, adjacency_list):
        self._adjacency_list = adjacency_list

    @property
    def vertex_list(self):
        return self._vertex_list

    @vertex_list.setter
    def vertex_list(self, vertex_list):
        self._vertex_list = vertex_list

    def add_vertex(self, vertex, heuristic_vector):
        """
        Add a vertex with specified heuristic vector
        """
        assert isinstance(vertex, int)
        #assert len(heuristic_vector) == 2
        assert (vertex not in self.vertex_list)
        self.vertex_list[vertex] = {
            'H' : heuristic_vector,
            'G' : [],
            'F' : [],
            'path' : []     # Note elements in G, F, path must be ordered as they correspond to one another
        }

    def add_directed_edge(self, vertex_A, vertex_B, cost_vector):
        """
        Add a directed edge from vertex_A to vertex_B specified cost_vector
        """
        assert isinstance(vertex_A, int)
        assert isinstance(vertex_B, int)
        assert (vertex_A in self.vertex_list)
        assert (vertex_B in self.vertex_list)
        #assert len(cost_vector) == 2
        self.adjacency_list[vertex_A] = {
            'next' : vertex_B,
            'cost' : cost_vector
        }

    def read_graph(self, vertex_path, edge_path):
        """
        Read vertices and edges from csv files

        Vertex file format -

            <vertex #>, <h1>, <h2>
            ...

        Edge file format -

            <vertex 1>, <vertex 2>, <cost 1>, <cost 2>
            ...
        
        """

        with open(vertex_path, 'r') as f1:
            reader = csv.reader(f1)
            for row in reader:
                heuristic_vector = tuple([int(i) for i in row[1:]])
                self.add_vertex(
                    vertex=int(row[0]),
                    heuristic_vector=heuristic_vector
                )

        with open(edge_path, 'r') as f2:
            reader = csv.reader(f2)
            for row in reader:
                cost_vector = tuple([int(i) for i in row[2:]])
                self.add_directed_edge(
                    vertex_A=int(row[0]),
                    vertex_B=int(row[1]),
                    cost_vector=cost_vector
                )

def main():
    
    vertex_path = "./input/toy_graph_vertex.csv"
    edge_path = "./input/toy_graph_edges.csv"

    G = graph(vertex_path=vertex_path, edge_path=edge_path)

    for vertex in G.vertex_list:
        print(G.vertex_list[vertex])

if __name__ == '__main__':
    main()