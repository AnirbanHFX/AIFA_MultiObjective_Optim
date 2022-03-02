import csv
import copy

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
        if vertex_A not in self.adjacency_list:
            self.adjacency_list[vertex_A] = []
        self.adjacency_list[vertex_A].append({
            'next' : vertex_B,
            'cost' : cost_vector
        })

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

        # Keep list of goal nodes
        self.goal_nodes = set(goal_list)

        # Create Open list initialized with the start node
        assert start_node in graph.vertex_list
        self.OPEN = set([start_node])
        self.graph.vertex_list[start_node]['G'].append((0, 0))
        # For each G associated with the node, assign an F[i] = G[i] + H
        self.graph.vertex_list[start_node]['F'].append(tuple([sum(x) for x in zip(
            self.graph.vertex_list[start_node]['G'][0],
            self.graph.vertex_list[start_node]['H']
        )]))
        self.graph.vertex_list[start_node]['path'].append([start_node])

        # Initialize empty lists for CLOSED, SOLUTION, LABEL
        self.SOLUTION = {}
        self.CLOSED = set([])
        self.LABEL = {}

    def cost_dominates(self, F1, F2):
        """
        If all elements of F1 <= F2 and at least 1 element of F1 < F2
            return True
        else
            return False
        """
        temp = False
        for elem1, elem2 in zip(F1, F2):
            if elem2 < elem1:
                return False    # all elem1 >= elem2 condition failed
            elif elem1 < elem2:
                temp = True     # at least 1 elem1 > elem2 satisfied
        return temp

    def _find_non_dominated(self, vertices):

        ND = []

        def not_dominated(vertex1, vertex2, verbose=0):
            """
            IF there is at least 1 cost vector in vertex1 that is not dominated by
            any cost vector in vertex2
            """
            if isinstance(vertex2, int):
                second_iterator = self.graph.vertex_list[vertex2]['F']
            else:
                second_iterator = list(vertex2)

            for F1 in self.graph.vertex_list[vertex1]['F']:
                flag = True
                for F2 in second_iterator:
                    if verbose == 1:
                        print("F2 : ", F2, "F1 : ", F1)
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
                    #### DEBUG #####
                    # if (vertex == 3) and (nd_flag is True):
                    #     print("DEBUG STEP")
                    #     print(other_vertex)
                    #     print(not_dominated(vertex, other_vertex, verbose=1))
                    #     z = input()
                    ################
            if nd_flag is True:
                #### DEBUG ####
                # if (vertex == 3):
                #     print("3: ",self.graph.vertex_list[vertex])
                #     for other in vertices:
                #         print(other, ": ", self.graph.vertex_list[other])
                #     z = input()
                ###############
                temp_ND.append(vertex)

        # Check whether a vertex that qualified the previous round is also not dominated by any existing solution
        for vertex in temp_ND:
            flag = True
            for node in self.SOLUTION:    #TODO: Check data structure of SOLUTION_COSTS and cost
                                                #      Current assumption - cost is a single vector (tuple)
                cost = self.SOLUTION[node]['G']    # Should return a list
                assert isinstance(cost, list)
                if not not_dominated(vertex, cost):
                    # vertex does not have a single F that is not dominated by cost
                    flag = False
                    break
            if flag is True:
                ND.append(vertex)

        return ND

    def _terminate(self):
        print(self.SOLUTION)

    def _choose_from_ND(self, ND):
        # TODO: Use domain specific heuristic to choose from ND
        return ND[0]

    def _remove_dominated_solutions(self):

        for vertex1 in self.SOLUTION:
            iterator1 = copy.deepcopy(self.SOLUTION[vertex1]['G'])
            for G1 in iterator1:
                for vertex2 in self.SOLUTION:
                    iterator2 = copy.deepcopy(self.SOLUTION[vertex2]['G'])
                    for G2 in iterator2:
                        if self.cost_dominates(G2, G1):
                            # G1 is fully dominated by G2
                            # So G1 must be removed
                            idx = self.SOLUTION[vertex1]['G'].index(G1)
                            self.SOLUTION[vertex1]['G'].pop(idx)
                            self.SOLUTION[vertex1]['path'].pop(idx)

    def _expand(self, node):

        def accrued_non_dominated_paths(n):
            ND = {'G' : [], 'path' : []}
            for G1, path1 in zip(self.graph.vertex_list[n]['G'], self.graph.vertex_list[n]['path']):
                flag = True
                for G2 in self.graph.vertex_list[n]['G']:
                    if G1 != G2:
                        if self.cost_dominates(G2, G1):
                            # G1 is fully dominated by G2
                            # So G1 must be removed
                            flag = False
                            break
                if flag is True:
                    ND['G'].append(G1)
                    ND['path'].append(path1)
            return ND

        if len(self.graph.adjacency_list[node]) == 0:
            return

        for next_node_dict in self.graph.adjacency_list[node]:
            next_node = next_node_dict['next']
            cost = next_node_dict['cost']

            if next_node not in self.OPEN and next_node not in self.CLOSED:
                # Add costs
                for G in self.graph.vertex_list[node]['G']:
                    nextG = tuple([sum(x) for x in zip(
                        G,
                        cost
                    )])
                    assert (len(nextG) == 2)
                    self.graph.vertex_list[next_node]['G'].append(nextG)
                # Establish backpointers
                for path in self.graph.vertex_list[node]['path']:
                    self.graph.vertex_list[next_node]['path'].append([*path, next_node])
                # Set LABEL(next_node, node)
                self.LABEL[(next_node, node)] = accrued_non_dominated_paths(next_node)  #TODO: Check correctness of LABEL(n,n') and LABEL(n',n)
                assert isinstance(self.LABEL[(next_node, node)]['G'], list)
                assert isinstance(self.LABEL[(next_node, node)]['path'], list)
                self.graph.vertex_list[next_node]['G'] = copy.deepcopy(self.LABEL[(next_node, node)]['G'])
                self.graph.vertex_list[next_node]['path'] = copy.deepcopy(self.LABEL[(next_node, node)]['path'])

                self.graph.vertex_list[next_node]['F'] = []
                for G in self.graph.vertex_list[next_node]['G']:
                    self.graph.vertex_list[next_node]['F'].append(tuple([sum(x) for x in zip(
                        G,
                        self.graph.vertex_list[next_node]['H']
                    )]))

                self.OPEN.add(next_node)

            else:
                for G, path in zip(self.graph.vertex_list[node]['G'], self.graph.vertex_list[node]['path']):
                    # Check new cost vectors 'nextG'
                    nextG = tuple([sum(x) for x in zip(
                        G,
                        cost
                    )])
                    assert (len(nextG) == 2)

                    if (next_node, node) not in self.LABEL:
                        self.LABEL[(next_node, node)] = {'G' : [], 'path' : []}
                    if nextG not in self.LABEL[(next_node, node)]['G']:
                        flag = True
                        for G2 in self.LABEL[(next_node, node)]['G']:
                            if self.cost_dominates(G2, nextG):
                                # G2 dominates nextG
                                # nextG is not a candidate to be added to LABEL
                                flag = False
                                break
                        if flag is True:
                            self.LABEL[(next_node, node)]['G'].append(nextG)
                            self.LABEL[(next_node, node)]['path'].append([*path, next_node])

                            assert nextG not in self.graph.vertex_list[next_node]['G']
                            if nextG not in self.graph.vertex_list[next_node]['G']:
                                self.graph.vertex_list[next_node]['G'].append(nextG)
                                self.graph.vertex_list[next_node]['path'].append([*path, next_node])

                            iterator = copy.deepcopy(self.LABEL[(next_node, node)]['G'])
                            for G2 in iterator:
                                if self.cost_dominates(nextG, G2):
                                    # new G strictly dominates G2
                                    # So remove G2
                                    idx = self.LABEL[(next_node, node)]['G'].index(G2)
                                    self.LABEL[(next_node, node)]['G'].pop(idx)
                                    self.LABEL[(next_node, node)]['path'].pop(idx)
                            
                            if next_node in self.CLOSED:
                                self.CLOSED.remove(next_node)
                                self.OPEN.add(next_node)
                
    def MOA(self, verbose=0):

        iter = 0
        while(True):

            if verbose==1:
                print("Iter: ", iter)
                print("OPEN: ", self.OPEN)
                print("CLOSED: ", self.CLOSED)
                print("SOLUTION: ", self.SOLUTION)

            iter += 1

            ###### FIND ND ######
            ND = self._find_non_dominated(self.OPEN)

            if verbose==1:
                print("ND: ", ND)

            ###### Terminate if ND is empty ######
            if len(ND) == 0:
                if verbose==1:
                    print("******************************")
                self._terminate()    #TODO: Complete this function
                break

            ###### Select ######
            n = self._choose_from_ND(ND)
            self.OPEN.remove(n)         # Remove from OPEN
            self.CLOSED.add(n)          # Add to CLOSED
            #TODO: Bookkeeping

            ####### Identify Solutions #######
            if n in self.goal_nodes:
                # If n is a goal node
                # Add its cost vectors to SOLUTIONS COSTS
                if n not in self.SOLUTION:
                    self.SOLUTION[n] = {'G' : [], 'path' : []}
                    assert isinstance(self.graph.vertex_list[n]['G'], list)
                    assert isinstance(self.graph.vertex_list[n]['path'], list)
                    assert len(self.graph.vertex_list[n]['G']) == len(self.graph.vertex_list[n]['path'])
                    self.SOLUTION[n]['G'].extend(self.graph.vertex_list[n]['G'])
                    self.SOLUTION[n]['path'].extend(self.graph.vertex_list[n]['path'])
                else:
                    assert isinstance(self.graph.vertex_list[n]['G'], list)
                    assert isinstance(self.graph.vertex_list[n]['path'], list)
                    assert len(self.graph.vertex_list[n]['G']) == len(self.graph.vertex_list[n]['path'])
                    self.SOLUTION[n]['G'].extend(self.graph.vertex_list[n]['G'])
                    self.SOLUTION[n]['path'].extend(self.graph.vertex_list[n]['path'])
                self._remove_dominated_solutions()
                if verbose == 1:
                    print("EXPANDED: None")
            else:
                # Expand
                if verbose == 1:
                    print("EXPANDED: ", n)
                self._expand(n)

            if verbose==1:
                print("******************************")



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
    
    vertex_path = "./input/test_graph_vertex.csv"
    edge_path = "./input/test_graph_edges.csv"

    G = Graph(vertex_path=vertex_path, edge_path=edge_path)
    moa = MOAsolver(G, 0, [8, 9, 10])

    F1 = (3, 2)
    F2 = (1, 1)
    print(moa.cost_dominates(F1, F2))

    moa.MOA(verbose=1)

if __name__ == '__main__':
    main()