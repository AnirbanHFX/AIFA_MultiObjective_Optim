import csv

class Graph:

    def __init__(self, vertex_path=None, edge_path=None):
        self._adjacency_list = {}
        self._vertex_list = {}
        if (vertex_path is not None) and (edge_path is not None):
            self.read_graph(vertex_path, edge_path)

    @property
    def adjacency_list(self):
        return self._adjacency_list

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

class MOAsolver:

    def __init__(self, graph, start_node, goal_list):

        # Store graph
        assert isinstance(graph, Graph)
        self.graph = graph

        # Create Open list initialized with the start node
        assert start_node in graph.vertex_list
        self.OPEN = [start_node]
        self.graph.vertex_list[start_node]['G'].append((0, 0))
        # For each G associated with the node, assign an F[i] = G[i] + H
        self.graph.vertex_list[start_node]['F'].append(tuple([sum(x) for x in zip(
            self.graph.vertex_list[start_node]['G'][0],
            self.graph.vertex_list[start_node]['H']
        )]))

        # Initialize empty lists for CLOSED, SOLUTION, SOLUTION_COSTS, LABEL
        self.SOLUTION_GOALS = []
        self.SOLUTION_COSTS = []
        self.CLOSED = []
        self.SOLUTION = []
        self.LABEL = []

    def cost_dominates(self, F1, F2):
        """
        If all elements of F1 >= F2 and at least 1 element of F1 > F2
            return True
        else
            return False
        """
        temp = False
        for elem1, elem2 in zip(F1, F2):
            if elem2 > elem1:
                return False    # all elem1 >= elem2 condition failed
            elif elem1 > elem2:
                temp = True     # at least 1 elem1 > elem2 satisfied
        return temp

    def _find_non_dominated(self, vertices):

        ND = []

        def not_dominated(vertex1, vertex2):
            """
            IF there is at least 1 cost vector in vertex1 that is not dominated by
            any cost vector in vertex2
            """
            if isinstance(vertex2, int):
                second_iterator = self.graph.vertex_list[vertex2]['F']
            else:
                second_iterator = list[vertex2]

            for F1 in self.graph.vertex_list[vertex1]['F']:
                flag = True
                for F2 in second_iterator:
                    if self.cost_dominates(F2, F1):
                        # F1 is fully dominated by F2
                        # So this F1 does not satisfy the criteria (there exists F1 which is not dominmated by any F2)
                        flag = False
                        break
                if flag is True:
                    # An F1 has been found which no F2 dominates
                    return True
                # else continue searching with next F1
            
            return False # No F1 has been found so return False

        # Check whether other a vertex is not dominated by any other potential solution represented by another vertex
        temp_ND = []
        for vertex in vertices:
            nd_flag = True
            for other_vertex in vertices:
                if vertex != other_vertex:
                    if not not_dominated(vertex, other_vertex):
                        # If vertex1 is not not-dominated by the other vertex, vertex1 does not belong in ND
                        nd_flag = False
                        break
            if nd_flag is True:
                temp_ND.append(vertex)

        # Check whether a vertex that qualified the previous round is also not dominated by any existing solution
        for vertex in temp_ND:
            flag = True
            for cost in self.SOLUTION_COSTS:    #TODO: Check data structure of SOLUTION_COSTS and cost
                                                #      Current assumption - cost is a single vector (tuple)
                if not not_dominated(vertex, [cost]):
                    # vertex does not have a single F that is not dominated by cost
                    flag = False
                    break
            if flag is True:
                ND.append(vertex)

        return ND

    def MOA(self):

        while(True):

            ###### FIND ND ######
            ND = self._find_non_dominated(self.OPEN)

    def test_ND(self):

        self.OPEN.append(1)
        self.graph.vertex_list[1]['G'].append((5, 6))
        # For each G associated with the node, assign an F[i] = G[i] + H
        self.graph.vertex_list[1]['F'].append(tuple([sum(x) for x in zip(
            self.graph.vertex_list[1]['G'][0],
            self.graph.vertex_list[1]['H']
        )]))

        self.OPEN.append(2)
        self.graph.vertex_list[2]['G'].append((4, 5))
        # For each G associated with the node, assign an F[i] = G[i] + H
        self.graph.vertex_list[2]['F'].append(tuple([sum(x) for x in zip(
            self.graph.vertex_list[2]['G'][0],
            self.graph.vertex_list[2]['H']
        )]))

        self.OPEN.append(3)
        self.graph.vertex_list[3]['G'].append((8, 10))
        self.graph.vertex_list[3]['G'].append((7, 9))
        # For each G associated with the node, assign an F[i] = G[i] + H
        for i in range(len(self.graph.vertex_list[3]['G'])):
            self.graph.vertex_list[3]['F'].append(tuple([sum(x) for x in zip(
                self.graph.vertex_list[3]['G'][i],
                self.graph.vertex_list[3]['H']
            )]))

        print(self.OPEN)
        print(self._find_non_dominated(self.OPEN))

def main():
    
    vertex_path = "./input/toy_graph_vertex.csv"
    edge_path = "./input/toy_graph_edges.csv"

    G = Graph(vertex_path=vertex_path, edge_path=edge_path)
    moa = MOAsolver(G, 0, [3])

    moa.test_ND()

if __name__ == '__main__':
    main()