import csv
import copy

class Graph:

    def __init__(self, start_node, num_vertex, edge_path=None):
        self._adjacency_list = {}
        self._state_list = {}
        self.start_node = start_node
        self._num_vertices = num_vertex
        if (edge_path is not None):
            self.read_graph(edge_path)

        self._heuristic_mem = {}

    @property
    def adjacency_list(self):
        return self._adjacency_list

    @adjacency_list.setter
    def adjacency_list(self, adjacency_list):
        self._adjacency_list = adjacency_list

    @property
    def state_list(self):
        return self._state_list

    @state_list.setter
    def state_list(self, state_list):
        self._state_list = state_list

    def add_state(self, state, heuristic_vector):
        """
        Add a state with specified heuristic vector
        """
        assert isinstance(state, tuple)
        #assert len(heuristic_vector) == 2
        assert (state not in self.state_list)
        self.state_list[state] = {
            'H' : heuristic_vector,
            'G' : [],
            'F' : []
        }

    def add_directed_edge(self, vertex_A, vertex_B, cost_vector):
        """
        Add a directed edge from vertex_A to vertex_B specified cost_vector
        """
        assert isinstance(vertex_A, int)
        assert isinstance(vertex_B, int)
        assert vertex_A < self._num_vertices and vertex_B < self._num_vertices
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
        self.state_list = {}
        self._heuristic_mem = {}

    def get_heuristic(self, state):
        """
        Get heuristic for a state
        Input -
            state : partial path traversed by travelling salesperson
        Methodology -
            If a new partial path is encountered which has a set of nodes that have not been seen before
                Compute shortest paths from all sources to all destinations through {(set of nodes not in partial path) V (last node in partial path)}
                Add these shortest paths to memory
                Return shortest path from (last node in partial path) to (goal node)
            Else
                Return shortest path which has been stored in memory
        """

        inf = 1e10
        
        state_set = frozenset(state)

        if state_set in self._heuristic_mem:
            if state[-1] in self._heuristic_mem[state_set]:
                return tuple(self._heuristic_mem[state_set][state[-1]][state[-1]][self.start_node])
        
        # Get the list of remaining vertices
        if state[-1] == self.start_node:
            remaining_vertices = [self.start_node]
        else:
            remaining_vertices = [self.start_node, state[-1]]
        for i in range(self._num_vertices):
            if i not in state_set:
                remaining_vertices.append(i)

        dist = {}

        # Initialize table of shortest paths
        for elem1 in remaining_vertices:
            assert isinstance(elem1, int)
            dist[elem1] = {}
            for elem2 in remaining_vertices:
                assert isinstance(elem2, int)
                if elem1 == elem2:
                    dist[elem1][elem2] = [0, 0]
                else:
                    flag = False
                    for next_dict in self.adjacency_list[elem1]:
                        if next_dict['next'] == elem2:
                            dist[elem1][elem2] = list(next_dict['cost'])
                            assert len(dist[elem1][elem2]) == 2
                            flag = True
                            break
                    if not flag:    # Direct path does not exist
                        dist[elem1][elem2] = [inf, inf]

        for elem3 in remaining_vertices:
            for elem1 in remaining_vertices:
                for elem2 in remaining_vertices:
                    if dist[elem1][elem3][0] + dist[elem3][elem2][0] < dist[elem1][elem2][0]:
                         dist[elem1][elem2][0] = dist[elem1][elem3][0] + dist[elem3][elem2][0]
                    if dist[elem1][elem3][1] + dist[elem3][elem2][1] < dist[elem1][elem2][1]:
                         dist[elem1][elem2][1] = dist[elem1][elem3][1] + dist[elem3][elem2][1]

        if state_set in self._heuristic_mem:
            # state[-1] not in self._heuristic_mem[state_set]
            self._heuristic_mem[state_set][state[-1]] = dist
        else:
            self._heuristic_mem[state_set] = {state[-1] : dist}

        return tuple(dist[state[-1]][self.start_node])


    def generate_state(self, state):
        """
        Appends state to state dict
        """
        assert isinstance(state, tuple)
        if state not in self.state_list:
            heuristic = self.get_heuristic(state)
            assert isinstance(heuristic, tuple)
            for elem in heuristic:
                assert elem < 1e8
            self.add_state(state, self.get_heuristic(state))

    def get_next_states(self, state):
        """
        Find next states
        Returns a list of dicts containing the following keys -
            'next' : a tuple describing the next state (partial path)
            'cost' : a tuple containing the bi-objective cost
        """
        assert isinstance(state, tuple)
        assert len(state) > 0
        if len(state) < self._num_vertices:
            state = list(state)
            state_set = set(state)
            current_vertex = state[-1]
            next_state_list = []
            for adjacent in self._adjacency_list[current_vertex]:
                if adjacent['next'] not in state_set:
                    temp = []
                    temp = copy.deepcopy(state)
                    temp.append(adjacent['next'])
                    temp = tuple(temp)
                    next_state_list.append({'next' : temp, 'cost' : adjacent['cost']})
            state = tuple(state)
            return next_state_list
        else:
            current_vertex = state[-1]
            for adjacent in self._adjacency_list[current_vertex]:
                if adjacent['next'] == self.start_node:
                    state = list(state)
                    temp = []
                    temp = copy.deepcopy(state)
                    temp.append(self.start_node)
                    temp = tuple(temp)
                    state = tuple(state)
                    return [{'next' : temp, 'cost' : adjacent['cost']}]
            return []

    def read_graph(self, edge_path):
        """
        Read edges from csv files

        Edge file format -

            <vertex 1>, <vertex 2>, <cost 1>, <cost 2>
            ...
        
        """

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

    def __init__(self, graph):

        # Store graph
        assert isinstance(graph, Graph)
        self.graph = graph

        self.start_node = self.graph.start_node
        start_node = self.start_node

        # Create Open list initialized with the start node
        assert start_node in graph.adjacency_list
        start_state = tuple([start_node])
        self.OPEN = set([start_state])

        self.graph.generate_state(start_state)
        self.graph.state_list[start_state]['G'].append((0, 0))
        self.graph.state_list[start_state]['H']
        # For each G associated with the node, assign an F[i] = G[i] + H
        self.graph.state_list[start_state]['F'].append(tuple([sum(x) for x in zip(
            self.graph.state_list[start_state]['G'][0],
            self.graph.state_list[start_state]['H']
        )]))

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

    def _find_non_dominated(self, states):

        ND = []

        def not_dominated(state1, state2, verbose=0):
            """
            IF there is at least 1 cost vector in state1 that is not dominated by
            any cost vector in state2
            """
            if isinstance(state2, tuple):
                second_iterator = self.graph.state_list[state2]['F']
            else:
                second_iterator = list(state2)

            for F1 in self.graph.state_list[state1]['F']:
                flag = True
                for F2 in second_iterator:
                    if verbose == 1:
                        print("F2 : ", F2, "F1 : ", F1)
                    if self.cost_dominates(F2, F1):
                        # F1 is fully dominated by F2
                        # So this F1 does not satisfy the criteria (there exists F1 which is not dominated by any F2)
                        flag = False
                        break
                if flag is True:
                    # An F1 has been found which no F2 dominates
                    return True
                # else continue searching with next F1
            
            return False # No F1 has been found so return False

        # Check whether a state is not dominated by any other
        # potential solution represented by another state
        temp_ND = []
        for state in states:
            nd_flag = True
            for other_state in states:
                if state != other_state:
                    if not not_dominated(state, other_state):
                        # If state1 is not not-dominated by the other state,
                        # state1 does not belong in ND
                        nd_flag = False
                        break
            if nd_flag is True:
                temp_ND.append(state)

        # Check whether a state that qualified the previous round 
        # is also not dominated by any existing solution
        for state in temp_ND:
            flag = True
            for goal in self.SOLUTION: 
                cost = self.SOLUTION[goal]['G']
                assert isinstance(cost, list)
                if not not_dominated(state, cost):
                    # state does not have a single F that is not dominated by cost
                    flag = False
                    break
            if flag is True:
                ND.append(state)

        return ND

    def _terminate(self):
        print("SOLUTIONS FOUND -")
        for state in self.SOLUTION:
            print("\tCost: ", self.SOLUTION[state]['G'][0], "\tPath: ", list(state))

    def _choose_from_ND(self, ND):
        # Return a random node from ND
        return ND[0]

    def _remove_dominated_solutions(self):

        pop_keys = []
        for state1 in self.SOLUTION:
            iterator1 = copy.deepcopy(self.SOLUTION[state1]['G'])
            for G1 in iterator1:
                for state2 in self.SOLUTION:
                    iterator2 = copy.deepcopy(self.SOLUTION[state2]['G'])
                    for G2 in iterator2:
                        if self.cost_dominates(G2, G1):
                            # G1 is fully dominated by G2
                            # So G1 must be removed
                            idx = self.SOLUTION[state1]['G'].index(G1)
                            self.SOLUTION[state1]['G'].pop(idx)
            if len(self.SOLUTION[state1]['G']) == 0:
                pop_keys.append(state1)
        for key in pop_keys:
            self.SOLUTION.pop(key)

    def _expand(self, node):

        def accrued_non_dominated_paths(n):
            ND = {'G' : []}
            for G1 in self.graph.state_list[n]['G']:
                flag = True
                for G2 in self.graph.state_list[n]['G']:
                    if G1 != G2:
                        if self.cost_dominates(G2, G1):
                            # G1 is fully dominated by G2
                            # So G1 must be removed
                            flag = False
                            break
                if flag is True:
                    ND['G'].append(G1)
            return ND

        next_state_list = self.graph.get_next_states(node)
        if len(next_state_list) == 0:
            return

        for next_node_dict in next_state_list:
            next_node = next_node_dict['next']
            cost = next_node_dict['cost']

            if next_node not in self.OPEN and next_node not in self.CLOSED:

                # Generate new TSP state at runtime
                self.graph.generate_state(next_node)

                # Add costs
                for G in self.graph.state_list[node]['G']:
                    nextG = tuple([sum(x) for x in zip(
                        G,
                        cost
                    )])
                    assert (len(nextG) == 2)
                    self.graph.state_list[next_node]['G'].append(nextG)

                self.LABEL[(next_node, node)] = accrued_non_dominated_paths(next_node)
                assert isinstance(self.LABEL[(next_node, node)]['G'], list)
                self.graph.state_list[next_node]['G'] = copy.deepcopy(self.LABEL[(next_node, node)]['G'])

                self.graph.state_list[next_node]['F'] = []
                for G in self.graph.state_list[next_node]['G']:
                    self.graph.state_list[next_node]['F'].append(tuple([sum(x) for x in zip(
                        G,
                        self.graph.state_list[next_node]['H']
                    )]))

                self.OPEN.add(next_node)

            else:
                for G in self.graph.state_list[node]['G']:
                    # Check new cost vectors 'nextG'
                    nextG = tuple([sum(x) for x in zip(
                        G,
                        cost
                    )])
                    assert (len(nextG) == 2)

                    if (next_node, node) not in self.LABEL:
                        self.LABEL[(next_node, node)] = {'G' : []}
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

                            assert nextG not in self.graph.state_list[next_node]['G']
                            if nextG not in self.graph.state_list[next_node]['G']:
                                self.graph.state_list[next_node]['G'].append(nextG)

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
                        
    def is_goal(self, state):

        assert(state[0] == self.start_node)
        if (len(state) == self.graph._num_vertices+1) and state[-1] == self.start_node:
            return True
        else:
            return False

                
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
                self._terminate()
                break

            ###### Select ######
            n = self._choose_from_ND(ND)
            self.OPEN.remove(n)         # Remove from OPEN
            self.CLOSED.add(n)          # Add to CLOSED

            ####### Identify Solutions #######
            if self.is_goal(n):
                # If n is a goal state
                # Add its cost vectors to SOLUTION
                if n not in self.SOLUTION:
                    self.SOLUTION[n] = {'G' : []}
                    assert isinstance(self.graph.state_list[n]['G'], list)
                    self.SOLUTION[n]['G'].extend(self.graph.state_list[n]['G'])
                else:
                    assert isinstance(self.graph.state_list[n]['G'], list)
                    self.SOLUTION[n]['G'].extend(self.graph.state_list[n]['G'])
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

    def __init__(self, graph):

        # Store graph
        assert isinstance(graph, Graph)
        self.graph = graph

        self.start_node = self.graph.start_node
        start_node = self.start_node

        # Create Open list initialized with the start state
        assert start_node in graph.adjacency_list
        start_state = tuple([start_node])
        self.OPEN = set([start_state])

        self.graph.generate_state(start_state)
        self.graph.state_list[start_state]['G'].append((0, 0))
        self.graph.state_list[start_state]['H']
        # For each G associated with the state, assign an F[i] = G[i] + H
        self.graph.state_list[start_state]['F'].append(tuple([sum(x) for x in zip(
            self.graph.state_list[start_state]['G'][0],
            self.graph.state_list[start_state]['H']
        )]))

        # Initialize empty lists for CLOSED, SOLUTION, BEST_COST
        self.BEST_COSTS = [{'node' : (-1,), 'G' : (1e10, 1e10)}]
        self.CLOSED = set([])

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
        Backtrack if there does not exist any F in node that is not-dominated by all Best-Costs
        """
        for F in self.graph.state_list[node]['F']:
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
        print("SOLUTIONS FOUND -")
        for item in self.BEST_COSTS:
            print("\tCost: ", item['G'], "\tPath: ", list(item['node']))

    def is_goal(self, state):

        assert(state[0] == self.start_node)
        if (len(state) == self.graph._num_vertices+1) and state[-1] == self.start_node:
            return True
        else:
            return False

    def DFBB(self, verbose=0):

        def accrued_non_dominated_paths(n):
            ND = {'G' : []}
            for G1 in self.graph.state_list[n]['G']:
                flag = True
                for G2 in self.graph.state_list[n]['G']:
                    if G1 != G2:
                        if self.cost_dominates(G2, G1):
                            # G1 is fully dominated by G2
                            # So G1 must be removed
                            flag = False
                            break
                if flag is True:
                    ND['G'].append(G1)
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

            # Pick a state from OPEN and put it in CLOSED
            node = self.OPEN.pop()
            self.CLOSED.add(node)

            # Decide whether to backtrack from state
            if self.backtrack(node):
                if verbose==1:
                    print("BACKTRACK: ", node)
                    print("BACKTRACKED COST: ", self.graph.state_list[node]['F'])
                    print("******************************")
                continue
            else:
                if node[-1] not in self.graph.adjacency_list:
                    # Leaf city detected
                    # Not possible in TSP hence must assert false
                    raise Exception("Graph cannot be traversed by salesperson")
                    # Backtrack
                    if verbose == 1:
                        print("BACKTRACK: ", node)
                        print("BACKTRACKED COST: ", self.graph.state_list[node]['F'])
                        print("******************************")
                    continue
                else:
                    if verbose==1:
                        print("EXPAND: ", node)

            # Get next states
            next_state_list = self.graph.get_next_states(node)

            # Expand this state if we do not backtrack
            for next_node_dict in next_state_list:

                next_node = next_node_dict['next']
                cost = next_node_dict['cost']

                if next_node not in self.graph.state_list:
                    # Generate new TSP state at runtime
                    self.graph.generate_state(next_node)

                # Compute G for next state
                for G in self.graph.state_list[node]['G']:
                    nextG = tuple([sum(x) for x in zip(
                        G,
                        cost
                    )])
                    assert (len(nextG) == 2)
                    self.graph.state_list[next_node]['G'].append(nextG)

                # Find non-dominated costs
                nd_costs = accrued_non_dominated_paths(next_node)
                assert isinstance(nd_costs['G'], list)
                self.graph.state_list[next_node]['G'] = copy.deepcopy(nd_costs['G'])

                # Update F for next_node
                self.graph.state_list[next_node]['F'] = []
                for G in self.graph.state_list[next_node]['G']:
                    self.graph.state_list[next_node]['F'].append(tuple([sum(x) for x in zip(
                        G,
                        self.graph.state_list[next_node]['H']
                    )]))

                self.OPEN.add(next_node)
                if next_node in self.CLOSED:
                    self.CLOSED.remove(next_node)

                # Update best costs if next_node is goal state
                if self.is_goal(next_node):
                    for G in self.graph.state_list[next_node]['G']:
                        self.BEST_COSTS.append(
                            {'node' : next_node, 
                            'G' : G
                            }
                        )
                    # Only keep those best costs that are not dominated by any other
                    self.BEST_COSTS = self.non_dominated_best_costs()

            if verbose==1:
                print("******************************")

class IDMOAsolver:

    def __init__(self, graph):

        # Store graph
        assert isinstance(graph, Graph)
        self.graph = graph
        self.start_node = self.graph.start_node
        self.start_state = tuple([self.start_node])

        self.SOLUTIONS = []

        self.reset()

    def reset(self):
        # Initialize empty sets for closed, open
        self.CLOSED = set([])
        start_state = tuple([self.start_node])
        self.OPEN = set([start_state])
        self.graph.reset_costs()
    
        self.graph.generate_state(start_state)
        self.graph.state_list[start_state]['G'].append((0, 0))
        self.graph.state_list[start_state]['H']
        # For each G associated with the node, assign an F[i] = G[i] + H
        self.graph.state_list[start_state]['F'].append(tuple([sum(x) for x in zip(
            self.graph.state_list[start_state]['G'][0],
            self.graph.state_list[start_state]['H']
        )]))

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
                if ND[i]['G'] == ND[j]['G'] and ND[i]['node'] == ND[j]['node']:
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
            assert len(self.graph.state_list[node]['F']) == 1
            cost = self.graph.state_list[node]['F'][0][objective_idx]
            if cost > threshold:
                return True
            else:
                return False
        else:
            if len(self.graph.state_list[node]['F']) == 1:
                cost = self.graph.state_list[node]['F'][0][objective_idx]
                if cost > threshold:
                    return True
                else:
                    return False
            else:
                cntr = 0
                remove_F = []
                for F in self.graph.state_list[node]['F']:
                    cost = F[objective_idx]
                    if cost > threshold:
                        remove_F.append(F)
                        cntr += 1
                if cntr == len(self.graph.state_list[node]['F']):
                    # All of the paths are worthless
                    self.graph.state_list[node]['F'] = []
                    return True      # So backtrack
                else:
                    # Some good paths exist, so don't backtrack but remove bad paths
                    for F in remove_F:
                        idx = self.graph.state_list[node]['F'].index(F)
                        self.graph.state_list[node]['F'].pop(idx)
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

    def is_goal(self, state):

        assert(state[0] == self.start_node)
        if (len(state) == self.graph._num_vertices+1) and state[-1] == self.start_node:
            return True
        else:
            return False

    def search(self, objective_idx, threshold, find_all=False, verbose=1):

        def accrued_non_dominated_paths(n, objective_idx):
            ND = {'G' : []}
            for G1 in self.graph.state_list[n]['G']:
                flag = True
                for G2 in self.graph.state_list[n]['G']:
                    if G1 != G2:
                        if self.cost_dominates(G2, G1, objective_idx):
                            # G1 is fully dominated by G2
                            # So G1 must be removed
                            flag = False
                            break
                if flag is True:
                    ND['G'].append(G1)
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

            flag = True

            # Check if node is goal node
            if self.is_goal(node) and not find_all:
                # Goal node found
                # Break
                BACKTRACKED.append({
                    'node' : node, 
                    'F' : self.graph.state_list[node]['F'],
                })
                break
            elif self.is_goal(node):
                # Goal node found
                # But do not break
                flag = False
                BACKTRACKED.append({
                    'node' : node, 
                    'F' : self.graph.state_list[node]['F'],
                })

            if self.backtrack(node, threshold, objective_idx) and flag:
                # Node exceeds threshold
                # Must backtrack
                if verbose==1:
                    print("BACKTRACK: ", node)
                    print("BACKTRACKED COST: ", self.graph.state_list[node]['F'])
                    print("******************************")
                BACKTRACKED.append({
                    'node' : node, 
                    'F' : self.graph.state_list[node]['F'],
                })
                continue

            # Continue if node is a leaf
            if node[-1] not in self.graph.adjacency_list:
                # Not possible for TSP so assert False
                raise Exception("Graph cannot be traversed by salesperson")
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

            # Get next nodes
            next_state_list = self.graph.get_next_states(node)

            for next_node_dict in next_state_list:

                next_node = next_node_dict['next']
                cost = next_node_dict['cost']

                if next_node not in self.graph.state_list:
                    # Generate new TSP state at runtime
                    self.graph.generate_state(next_node)

                # Compute G for next node
                for G in self.graph.state_list[node]['G']:
                    nextG = tuple([sum(x) for x in zip(
                        G,
                        cost
                    )])
                    assert (len(nextG) == 2)
                    self.graph.state_list[next_node]['G'].append(nextG)

                # Find non-dominated costs
                nd_costs = accrued_non_dominated_paths(next_node, objective_idx)
                assert isinstance(nd_costs['G'], list)
                self.graph.state_list[next_node]['G'] = copy.deepcopy(nd_costs['G'])

                # Update F for next_node
                self.graph.state_list[next_node]['F'] = []
                for G in self.graph.state_list[next_node]['G']:
                    self.graph.state_list[next_node]['F'].append(tuple([sum(x) for x in zip(
                        G,
                        self.graph.state_list[next_node]['H']
                    )]))

                self.OPEN.add(next_node)
                if next_node in self.CLOSED:
                    self.CLOSED.remove(next_node)

        return BACKTRACKED

    def IDA_1(self, verbose=0):

        threshold = self.graph.state_list[self.start_state]['F'][0][0]

        while True:

            self.reset()
            assert isinstance(threshold, int)
            BACKTRACKED = self.search(0, threshold, False, verbose = verbose)       ###########################################################
            assert len(BACKTRACKED) > 0

            min_threshold = 1e10
            min_idx = -1
            idx = 0
            for item in BACKTRACKED:
                if self.is_goal(item['node']):
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
            if self.is_goal(item['node']):
                SOLUTIONS.append(item)

        return SOLUTIONS

    def IDMOA(self, verbose=0):

        # Step 1: Perform IDA* on the first objective
        # to obtain the optimum solution in terms of the same
        item = self.IDA_1(verbose=verbose)

        # Get the threshold for the second objective
        new_threshold = int(item['F'][0][1])

        self.reset()
        unfiltered_solutions = self.IDA_2(new_threshold, verbose=verbose)

        for item in unfiltered_solutions:
            for F in item['F']:
                self.SOLUTIONS.append({
                    'node' : item['node'],
                    'G' : F,
                })

        self.SOLUTIONS = self.non_dominated_best_costs()

        print("SOLUTIONS FOUND -")
        for item in self.SOLUTIONS:
            print("\tCost: ", item['G'], "\tPath: ", list(item['node']))

def main():
    
    edge_path = "./input/tsp_final_edges.csv"
    num_vertices = 8

    G1 = Graph(0, num_vertices, edge_path=edge_path)
    moa = MOAsolver(G1)

    print("Using MOA*")
    moa.MOA(verbose=0)

    G2 = Graph(0, num_vertices, edge_path=edge_path)
    dfbb = DFBBsolver(G2)

    print("\nUsing DFBB")
    dfbb.DFBB(verbose=0)

    G3 = Graph(0, num_vertices, edge_path=edge_path)
    idmoa = IDMOAsolver(G3)

    print("\nUsing IDMOA*")
    idmoa.IDMOA(verbose=0)

if __name__ == '__main__':
    main()