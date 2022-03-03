from glob import glob1


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
    def get(self, idx):
        if idx == 0:
            return self.x
        else:
            return self.y

class edge:
    def __init__(self, u, v, cv = cost_vector()):
        self.cv = cv
        self.u = u
        self.v = v


class graph:
    def __init__(self, adj_list = None, num_nodes = 5, start = 0, goals = [1], H = None, F = None, G = None):
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
            self.G = [cost_vector(1e20,1e20) for i in range(self.num_nodes)]
            self.G[start] = cost_vector(0,0)
        else:
            self.G = G
        if F == None:
            self.F = [cost_vector(self.G[i].x+self.H[i].x, self.G[i].y+self.H[i].y) for i in range(self.num_nodes)]
        else:
            self.F = F

    def clear_all(self):
        self.G = [cost_vector(1e20,1e20) for i in range(self.num_nodes)]
        self.G[start] = cost_vector(0,0)
        self.F = [cost_vector(self.G[i].x+self.H[i].x, self.G[i].y+self.H[i].y) for i in range(self.num_nodes)]

    def dfbb(self, thresh, node, objective_no, backtracked_nodes, goals_discovered):
        print(node, self.F[node].get(objective_no), thresh)
        # abc = input("abc")
        #True => backtracked
        if self.F[node].get(objective_no)>thresh:
            backtracked_nodes.append(node)
            return False
        if node in self.goals:
            goals_discovered.add(node)
            return True
        goalFound = False
        for edg in self.adj_list[node]:
            u = edg.u
            assert(u == node)
            v = edg.v
            cv = edg.cv
            new_g = cost_vector(self.G[node].x+cv.x,self.G[node].y+cv.y)
            new_f = cost_vector(new_g.x+cv.x, new_g.y+cv.y)
            if new_f.get(objective_no)<self.F[v].get(objective_no):
                self.F[v] = new_f
                self.G[v] = new_g
            goalFound = self.dfbb(thresh, v, objective_no, backtracked_nodes, goals_discovered)
            if goalFound and objective_no == 0:
                return True
            
            return goalFound
    def idmoa(self):
        thresh = self.H[self.start].x
        backtracked_nodes = []
        goal_found = False
        goals_discovered = set()
        while not goal_found:
            goal_found = self.dfbb(thresh, self.start, 0, backtracked_nodes, goals_discovered)
            if goal_found:
                break
            assert(len(goals_discovered) == 0)
            thresh = min([self.F[x].get(0) for x in backtracked_nodes])
            backtracked_nodes = []
            self.clear_all()
        
        max_thresh = min([self.F[x].get(1) for x in self.goals])
        print(max_thresh)
        abc = input("abc")
        backtracked_nodes = []
        self.clear_all()
        self.dfbb(max_thresh, self.start, 1, backtracked_nodes, goals_discovered)
        goal_costs = [self.F[x] for x in goals_discovered]
        return goal_costs

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
    print(g.idmoa())