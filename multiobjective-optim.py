import csv
import copy
from xml.etree.ElementPath import find

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

    def reset_costs(self):
        """
        Reset costs and paths accrued over the course of some prior search algorithm
        """
        for vertex in self.vertex_list:
            assert isinstance(self.vertex_list[vertex], dict)
            assert 'G' in self.vertex_list[vertex]
            assert 'F' in self.vertex_list[vertex]
            assert 'H' in self.vertex_list[vertex]
            assert 'path' in self.vertex_list[vertex]
            self.vertex_list[vertex]['G'] = []
            self.vertex_list[vertex]['F'] = []
            self.vertex_list[vertex]['path'] = []

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
            if nd_flag is True:
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
        print("SOLUTIONS FOUND -")
        for vertex in self.SOLUTION:
            for G, path in zip(self.SOLUTION[vertex]['G'], self.SOLUTION[vertex]['path']):
                print("Goal: ", vertex, "\tCost: ", G, "\tPath: ", path)

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

class DFBBsolver:

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

        # Initialize empty lists for CLOSED, SOLUTION, BEST_COST
        self.BEST_COSTS = [{'node' : -1, 'G' : (1e10, 1e10), 'path' : []}]
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

    def cost_equivalent(self, F1, F2):
        """
        Check if two costs are equivalent
        """
        return F1 == F2

    def non_dominated_best_costs(self):
        """
        Find set of non-dominated costs in BEST_COSTS
        """
        ND = []
        for item1 in self.BEST_COSTS:
            nd_flag = True
            for item2 in self.BEST_COSTS:
                if item1 != item2:
                    if self.cost_dominates(item2['G'], item1['G']):
                        # item2 dominates item1
                        # item1 is not a candidate for ND
                        nd_flag = False
                        break
            if nd_flag is True:
                # No other item in BEST_COSTS dominates item1
                # item1 is non-dominated
                ND.append(item1)
        return ND

    def backtrack(self, node):
        """
        Backtrack if there does not exist any F in node that is not-dominated by any Best-Cost
        """
        for F in self.graph.vertex_list[node]['F']:
            flag = True
            for cost in self.BEST_COSTS:
                bestG = cost['G']
                if self.cost_dominates(bestG, F):
                    # BestG dominates this F
                    # So search for another F
                    flag = False
            if flag is True:
                # Found an F that is not dominated by any BestG
                # Should Expand
                return False # Expand
        # Could not find an F that is not dominated by any BestG
        # Backtrack
        return True

    def terminate(self):
        self.BEST_COSTS = self.non_dominated_best_costs()
        print("SOLUTIONS FOUND -")
        for item in self.BEST_COSTS:
            print("Goal: ", item['node'], "\tCost: ", item['G'], "\tPath: ", item['path'])

    def non_dominated_best_costs(self):
        """
        Find set of non-dominated costs in BEST_COSTS
        """
        ND = []
        for item1 in self.BEST_COSTS:
            nd_flag = True
            for item2 in self.BEST_COSTS:
                if item1 != item2:
                    if self.cost_dominates(item2['G'], item1['G']):
                        # item2 dominates item1
                        # item1 is not a candidate for ND
                        nd_flag = False
                        break
            if nd_flag is True:
                # No other item in BEST_COSTS dominates item1
                # item1 is non-dominated
                ND.append(item1)
        pop_idx = set()
        clean_ND = []
        for i in range(len(ND)):
            for j in range(i+1, len(ND)):
                if ND[i]['G'] == ND[j]['G'] and ND[i]['path'] == ND[j]['path']:
                    pop_idx.add(i)
        for i in range(len(ND)):
            if i not in pop_idx:
                clean_ND.append(ND[i])
        return clean_ND

    def DFBB(self, verbose=0):

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

        iter = 0
        while True:

            if verbose==1:
                print("Iter: ", iter)
                print("OPEN: ", self.OPEN)
                print("CLOSED: ", self.CLOSED)
                print("BEST_COSTS: ", self.BEST_COSTS)

            iter += 1

            # Break if Open list is empty
            if len(self.OPEN) == 0:
                if verbose==1:
                    print("******************************")
                self.terminate()
                break

            # Pick a node from OPEN and put it in CLOSED
            node = self.OPEN.pop()
            self.CLOSED.add(node)

            if node in self.goal_nodes:
                for G, path in zip(self.graph.vertex_list[node]['G'], self.graph.vertex_list[node]['path']):
                    self.BEST_COSTS.append({
                        'node' : node, 
                        'G' : G, 
                        'path' : path})

            # Decide whether to backtrack from node
            if self.backtrack(node):
                if verbose==1:
                    print("BACKTRACK: ", node)
                    print("BACKTRACKED COST: ", self.graph.vertex_list[node]['F'])
                    print("******************************")
                continue
            else:
                if node not in self.graph.adjacency_list:
                    # Leaf node detected
                    # Backtrack
                    if verbose == 1:
                        print("BACKTRACK: ", node)
                        print("BACKTRACKED COST: ", self.graph.vertex_list[node]['F'])
                        print("******************************")
                    continue
                else:
                    if verbose==1:
                        print("EXPAND: ", node)

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

class IDMOAsolver:

    def __init__(self, graph, start_node, goal_list):

        # Store graph
        assert isinstance(graph, Graph)
        self.graph = graph
        self.start_node = start_node

        # Keep list of goal nodes
        self.goal_nodes = set(goal_list)

        # Create Open list initialized with the start node
        assert start_node in graph.vertex_list

        self.SOLUTIONS = []
        self.LABEL = {}

        self.reset()

    def reset(self):
        # Initialize empty sets for closed, open
        self.CLOSED = set([])
        self.OPEN = set([self.start_node])
        self.graph.reset_costs()
        self.graph.vertex_list[self.start_node]['G'].append((0, 0))
        # For each G associated with the node, assign an F[i] = G[i] + H
        self.graph.vertex_list[self.start_node]['F'].append(tuple([sum(x) for x in zip(
            self.graph.vertex_list[self.start_node]['G'][0],
            self.graph.vertex_list[self.start_node]['H']
        )]))
        self.graph.vertex_list[self.start_node]['path'].append([self.start_node])

    def non_dominated_best_costs(self):
        """
        Find set of non-dominated costs in BEST_COSTS
        """
        ND = []
        for item1 in self.SOLUTIONS:
            nd_flag = True
            for item2 in self.SOLUTIONS:
                if item1 != item2:
                    if self.cost_dominates(item2['G'], item1['G'], objective_idx=1):
                        # item2 dominates item1
                        # item1 is not a candidate for ND
                        nd_flag = False
                        break
            if nd_flag is True:
                # No other item in BEST_COSTS dominates item1
                # item1 is non-dominated
                ND.append(item1)
        pop_idx = set()
        clean_ND = []
        for i in range(len(ND)):
            for j in range(i+1, len(ND)):
                if ND[i]['G'] == ND[j]['G'] and ND[i]['path'] == ND[j]['path']:
                    pop_idx.add(i)
        for i in range(len(ND)):
            if i not in pop_idx:
                clean_ND.append(ND[i])
        return clean_ND

    def backtrack(self, node, threshold, objective_idx):
        """
        Decide whether to backtrack from a node based on threshold
        """
        assert objective_idx <= 1 and objective_idx >= 0
        assert isinstance(threshold, int)
        if objective_idx == 0:
            #assert len(self.graph.vertex_list[node]['G']) == 1
            cost = self.graph.vertex_list[node]['G'][0][objective_idx]
            if cost > threshold:
                return True
            else:
                return False
        else:
            if len(self.graph.vertex_list[node]['G']) == 1:
                cost = self.graph.vertex_list[node]['G'][0][objective_idx]
                if cost > threshold:
                    return True
                else:
                    return False
            else:
                cntr = 0
                remove_G = []
                for G in self.graph.vertex_list[node]['G']:
                    cost = G[objective_idx]
                    if cost > threshold:
                        remove_G.append(G)
                        cntr += 1
                if cntr == len(self.graph.vertex_list[node]['G']):
                    # All of the paths are worthless
                    # self.graph.vertex_list[node]['G'] = []
                    # self.graph.vertex_list[node]['path'] = []
                    return True      # So backtrack
                else:
                    # Some good paths exist, so don't backtrack but remove bad paths
                    for G in remove_G:
                        idx = self.graph.vertex_list[node]['G'].index(G)
                        self.graph.vertex_list[node]['G'].pop(idx)
                        self.graph.vertex_list[node]['path'].pop(idx)
                    return False    # Do not backtrack

    def cost_dominates(self, F1, F2, objective_idx):
        """
        Return True if F1 dominates F2 based on objective idx
        """
        if objective_idx == 0:
            return F1[objective_idx] < F2[objective_idx]
        else:
            temp = False
            for elem1, elem2 in zip(F1, F2):
                if elem2 < elem1:
                    return False    # all elem1 >= elem2 condition failed
                elif elem1 < elem2:
                    temp = True     # at least 1 elem1 > elem2 satisfied
            return temp

    def search(self, objective_idx, threshold, find_all=False, verbose=1):

        def accrued_non_dominated_paths(n, objective_idx):
            ND = {'G' : [], 'path' : []}
            for G1, path1 in zip(self.graph.vertex_list[n]['G'], self.graph.vertex_list[n]['path']):
                flag = True
                for G2 in self.graph.vertex_list[n]['G']:
                    if G1 != G2:
                        if self.cost_dominates(G2, G1, objective_idx):
                            # G1 is fully dominated by G2
                            # So G1 must be removed
                            flag = False
                            break
                if flag is True:
                    ND['G'].append(G1)
                    ND['path'].append(path1)
            return ND

        BACKTRACKED = []

        iter = 1
        while True:

            if verbose==1:
                print("Iter: ", iter)
                print("OPEN: ", self.OPEN)
                print("CLOSED: ", self.CLOSED)
                iter += 1

            if len(self.OPEN) == 0:
                if verbose==1:
                    print("******************************")
                # Open list is empty, break
                break

            node = self.OPEN.pop()
            self.CLOSED.add(node)

            if self.backtrack(node, threshold, objective_idx):
                # Node exceeds threshold
                # Must backtrack
                if verbose==1:
                    print("BACKTRACK: ", node)
                    print("BACKTRACKED COST: ", self.graph.vertex_list[node]['F'])
                    print("******************************")
                BACKTRACKED.append({
                    'node' : node, 
                    'F' : self.graph.vertex_list[node]['F'],
                    'path' : self.graph.vertex_list[node]['path']
                })
                continue

            # Check if next node is goal node
            if node in self.goal_nodes and not find_all:
                # Goal node found
                # Break
                BACKTRACKED.append({
                    'node' : node, 
                    'F' : self.graph.vertex_list[node]['F'],
                    'path' : self.graph.vertex_list[node]['path']
                })
                break
            elif node in self.goal_nodes:
                # Goal node found
                # But do not break
                BACKTRACKED.append({
                    'node' : node, 
                    'F' : self.graph.vertex_list[node]['F'],
                    'path' : self.graph.vertex_list[node]['path']
                })

            # Continue if node is a leaf
            if node not in self.graph.adjacency_list:
                if verbose==1:
                    print("BACKTRACK: ", node)
                    print("BACKTRACKED COST: ", self.graph.vertex_list[node]['F'])
                    print("******************************")
                BACKTRACKED.append({
                    'node' : node, 
                    'F' : self.graph.vertex_list[node]['F'],
                    'path' : self.graph.vertex_list[node]['path']
                })
                continue

            if verbose==1:
                print("EXPAND: ", node)
                print("******************************")

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
                    self.LABEL[(next_node, node)] = accrued_non_dominated_paths(next_node, objective_idx)  #TODO: Check correctness of LABEL(n,n') and LABEL(n',n)
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
                                if self.cost_dominates(G2, nextG, objective_idx):
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
                                    if self.cost_dominates(nextG, G2, objective_idx):
                                        # new G strictly dominates G2
                                        # So remove G2
                                        idx = self.LABEL[(next_node, node)]['G'].index(G2)
                                        self.LABEL[(next_node, node)]['G'].pop(idx)
                                        self.LABEL[(next_node, node)]['path'].pop(idx)
                                
                                if next_node in self.CLOSED:
                                    self.CLOSED.remove(next_node)
                                    self.OPEN.add(next_node)

        return BACKTRACKED

    def IDA_1(self, verbose=0):

        threshold = self.graph.vertex_list[self.start_node]['F'][0][0]

        while True:

            self.reset()
            assert isinstance(threshold, int)
            BACKTRACKED = self.search(0, threshold, False, verbose = verbose)
            assert len(BACKTRACKED) > 0

            min_threshold = 1e10
            min_idx = -1
            idx = 0
            for item in BACKTRACKED:
                if item['node'] in self.goal_nodes:
                    return item
                assert len(item['F']) == 1
                if item['F'][0][0] < min_threshold:
                    min_threshold = item['F'][0][0]
                    min_idx = idx
                idx += 1

            threshold = min_threshold

    def IDA_2(self, max_threshold, verbose=0):

        SOLUTIONS = []

        self.reset()
        assert isinstance(max_threshold, int)
        BACKTRACKED = self.search(1, max_threshold, True, verbose=verbose)
        assert len(BACKTRACKED) > 0

        for item in BACKTRACKED:
            if item['node'] in self.goal_nodes:
                SOLUTIONS.append(item)

        return SOLUTIONS

    def IDMOA(self, verbose=0):

        # Step 1: Perform IDA* on the first objective
        # to obtain the optimum solution in terms of the same
        item = self.IDA_1(verbose=verbose)

        # Get the threshold for the second objective
        new_threshold = item['F'][0][1]

        self.reset()
        unfiltered_solutions = self.IDA_2(new_threshold, verbose=verbose)

        for item in unfiltered_solutions:
            for F, path in zip(item['F'], item['path']):
                self.SOLUTIONS.append({
                    'node' : item['node'],
                    'G' : F,
                    'path' : path
                })

        self.SOLUTIONS = self.non_dominated_best_costs()

        print("SOLUTIONS FOUND -")
        for item in self.SOLUTIONS:
            print("Goal: ", item['node'], "\tCost: ", item['G'], "\tPath: ", item['path'])


def main():
    
    vertex_path = "./input/normal_vertex_test.csv"
    edge_path = "./input/normal_graph_test.csv"

    G1 = Graph(vertex_path=vertex_path, edge_path=edge_path)
    moa = MOAsolver(G1, 0, [8, 9])

    moa.MOA(verbose=0)

    G2 = Graph(vertex_path=vertex_path, edge_path=edge_path)
    dfbb = DFBBsolver(G2, 0, [8, 9])

    dfbb.DFBB(verbose=0)

    G3 = Graph(vertex_path=vertex_path, edge_path=edge_path)
    idmoa = IDMOAsolver(G3, 0, [8, 9])

    idmoa.IDMOA(verbose=0)

if __name__ == '__main__':
    main()