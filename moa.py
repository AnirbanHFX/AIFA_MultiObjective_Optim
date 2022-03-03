#TODO copy.deepcopy()
import copy
def not_nd(cv_list_1, cv_list_2):
    #if all f_2 in cv_list_2 gets dominated by some f_1 in cv_list_1 returns true

    for f_1 in cv_list_1:
        all_dominated = True
        for f_2 in cv_list_2:
            if not f_1.dominates(f_2):
                all_dominated = False
                break
        if all_dominated:
            return True
    return False

def remove_dominated(solution):
    # solution list of {"cost_vector":None,"path": [1,2,3,some goal node] }
    remove_list = []
    for i in range(len(solution)):
        for j in range(len(solution)):
            if i==j:
                continue
            if solution[i]["cost_vector"].dominates(solution[j]["cost_vector"]):
                if not (j in remove_list):
                    remove_list.append(solution[j])
    for x in remove_list:
        print("-----")
        print(x)
        print(solution)
        if x in solution:
            solution.remove(x)
    return solution
class cost_vector:
    def __init__(self, x = 1, y = 1):
        self.x = x
        self.y = y
    def dominates(self,other): #a.dominates(b) means a is better than b 
        if self.x<=other.x and self.y <= other.y:
            if self.x == other.x and self.y == other.y:
                return False
            else:# (<,=) (=,<), (<,<)
                return True
        else:
            return False
    def __eq__(self, other):
        if not isinstance(other, cost_vector):
            return NotImplemented
        return self.x == other.x and self.y == other.y
    def __repr__(self):
        rep = 'cost_vector(' + str(self.x) + ',' + str(self.y) + ')'
        return rep
class edge:
    def __init__(self, u, v, cv = cost_vector()):
        self.cv = cv
        self.u = u
        self.v = v
class graph:
    def __init__(self, adj_list = None, num_nodes = 5, H = None, G = None, F = None, start = 0, goals = [1]):
        self.start = start
        self.goals = goals
        self.num_nodes = num_nodes
        if adj_list == None:
            self.adj_list = [ [edge(i,(i+1)%self.num_nodes, cost_vector(1,1))] for i in range(num_nodes)]
            # self.adj_list = [[edge(0,1,cost_vector(5,6)), edge(0,2,cost_vector(4,5))], [edge(1,3, cost_vector(3,4))],[edge(2,3, cost_vector(3,4))] ]
        else:
            self.adj_list = adj_list
        if H == None:
            self.H = [cost_vector(0,0) for i in range(self.num_nodes)]
        else:
            self.H = H
        if G == None:
            self.G = [[] for i in range(self.num_nodes)] # for each node[[path, cost_vector]
            self.G[self.start].append({"path": [self.start], "cost_vector": cost_vector(0,0)})
        else:
            self.G = G
        if F == None:
            self.F = [[] for i in range(self.num_nodes)] #for each node [(path, cost_vector)]
            self.F[self.start].append({"path": [self.start], "cost_vector": self.H[self.start]})
        else:
            self.F = F
        self.OPEN = [self.start]
        self.CLOSED = []
        self.SOLUTION = []
        
    
    def find_nd(self):
        nd = []
        for node in self.OPEN:
            flag = 0 #0-> belongs to nd, 1->doesn't belong to nd
            for node_2 in self.OPEN:
                if node == node_2:
                    continue
                if not_nd([f["cost_vector"] for f in self.F[node_2]],[f["cost_vector"] for f in self.F[node]]): #if all f in F[node] gets dominated by some f2 in F[node_2]
                    flag = 1
            if not_nd([f["cost_vector"] for f in self.SOLUTION], [f["cost_vector"] for f in self.F[node]]):
                flag = 1
            if flag == 0:
                nd.append(node)
        return nd
    def chose_from_nd(self, nd):
        min_node = nd[0]
        #TODO chose min
        return min_node
    
    def terminate(self):
        print("here")
        for f in self.SOLUTION:
            print(f["path"])
            print(f["cost_vector"].x, f["cost_vector"].y)
        pass
    
    def expand(self,node):
    # [[edge(0,1,cost_vector(5,6)), edge(0,2,cost_vector(4,5))],
        for outgoing_edge in self.adj_list[node]:
            u = outgoing_edge.u
            v = outgoing_edge.v
            assert(u == node)
            cv = outgoing_edge.cv
            if (not v in self.OPEN) and (not v in self.CLOSED):
                # parent ka G usme cv add karenge G(v) usme heuristics add karenge
                #G[node] = [{"path": [n1, n2], "cost_vector": cost_vector(1,2)}, {"path": [n3, n2], "cost_vector": cost_vector(2,2)}]
                for g in self.G[node]:
                    new_g = copy.deepcopy(g)
                    new_g["path"].append(v)
                    new_g["cost_vector"]=cost_vector(g["cost_vector"].x+cv.x, g["cost_vector"].y+cv.y)
                    self.G[v].append(new_g)
                remove_dominated(self.G[v])
                for g in self.G[v]:
                    new_f = copy.deepcopy(g)
                    new_f["cost_vector"].x += self.H[v].x
                    new_f["cost_vector"].y += self.H[v].y
                    self.F[v].append(new_f)
                assert(v not in self.OPEN)
                self.OPEN.append(v)
            else:
                for g in self.G[node]:
                    new_g = copy.deepcopy(g)
                    new_g["path"].append(v)
                    new_g["cost_vector"]=cost_vector(g["cost_vector"].x+cv.x, g["cost_vector"].y+cv.y)
                    self.G[v].append(new_g)
                remove_dominated(self.G[v])
                for g in self.G[v]:
                    new_f = copy.deepcopy(g)
                    new_f["cost_vector"].x += self.H[v].x
                    new_f["cost_vector"].y += self.H[v].y
                    self.F[v].append(new_f)
                if v in self.CLOSED:
                    self.CLOSED.remove(v)
                    self.OPEN.append(v)

if __name__=="__main__":
    num_nodes = 11
    adj_list = []
    adj_list.append([edge(0,1,cost_vector(2,3)), edge(0,2,cost_vector(4,2))])
    adj_list.append([edge(1,3,cost_vector(3,2)), edge(1,4,cost_vector(2,3))])
    adj_list.append([edge(2,4,cost_vector(3,2)), edge(2,5,cost_vector(3,2))])
    adj_list.append([edge(3,6,cost_vector(5,7))])
    adj_list.append([edge(4,6,cost_vector(1,3)), edge(4,7,cost_vector(1,1))])
    adj_list.append([edge(5,7,cost_vector(5,2))])
    adj_list.append([edge(6,8,cost_vector(2,5)), edge(6,9,cost_vector(9,8))])
    adj_list.append([edge(7,9,cost_vector(6,8)), edge(7,10,cost_vector(4,3))])
    adj_list.append([])
    adj_list.append([])
    adj_list.append([])
    H = []
    H.append(cost_vector(2,3))
    H.append(cost_vector(2,3))
    H.append(cost_vector(3,2))
    H.append(cost_vector(5,7))
    H.append(cost_vector(1,1))
    H.append(cost_vector(5,2))
    H.append(cost_vector(2,5))
    H.append(cost_vector(4,3))
    H.append(cost_vector(0,0))
    H.append(cost_vector(0,0))
    H.append(cost_vector(0,0))
    start = 0
    goals = [8,9,10]
    # G = []
    # G.append({"path": [], "cost_vector": cost_vector(0,0)})
    # G.append({"path": [], "cost_vector": cost_vector(5,6)})
    # G.append({"path": [], "cost_vector": cost_vector(0,0)})
    # G.append({"path": [], "cost_vector": cost_vector(0,0)})
    g = graph(num_nodes = num_nodes, adj_list = adj_list, H = H, start = start, goals = goals) 
    while(True):
        ND = g.find_nd()
        if ND == []:
            g.terminate()
            break
        else:
            n = g.chose_from_nd(ND)
            g.OPEN.remove(n)
            assert((n in g.OPEN) == False)
            g.CLOSED.append(n)
            #Do bookkeeping???
            if n in g.goals:
                for f in g.F[n]:
                    if f not in g.SOLUTION:
                        print(f)
                        print(g.SOLUTION)
                        g.SOLUTION.append(f)
                g.SOLUTION = remove_dominated(g.SOLUTION)
            else:
                g.expand(n)
#print solution etc