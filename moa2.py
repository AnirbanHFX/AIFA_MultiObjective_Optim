#TODO copy.deepcopy()
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

class edge:
    def __init__(self, u, v, cost_vector = cost_vector()):
        self.cost_vector = cost_vector
        self.u = u
        self.v = v
class graph:
    def __init__(self, adj_list = None, num_nodes = 5, H = None, G = None, F = None, start = 0, goals = [1]):
        self.start = start
        self.goals = goals
        self.num_nodes = num_nodes
        if H == None:
            self.H = [cost_vector(0,0) for i in range(graph.num_nodes)]
        else:
            self.H = H
        if G == None:
            self.G = [[] for i in range(graph.num_nodes)] # for each node[[path, cost_vector]
        else:
            self.G = G
        if F == None:
            self.F = [[] for i in range(graph.num_nodes)] #for each node [(path, cost_vector)]
        else:
            self.F = F
        self.OPEN = [g.start]
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

    
    
if __name__=="__main__":
    g = graph() 
    while(True):
        ND = g.find_nd()
        if ND == []:
            g.terminate()
        else:
            n = g.chose_from_nd(ND)
            g.OPEN.remove(n)
            assert((n in g.OPEN) == False)
            g.CLOSED.append(n)
            #Do bookkeeping???
            if n in g.goals:
                for f in g.F[n]:
                    g.SOLUTION.append(f)
                g.SOLUTION = remove_dominated(g.SOLUTION)
            else:
                g.expand(n)
#print solution etc