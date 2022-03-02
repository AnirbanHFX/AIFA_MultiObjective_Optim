#
# Name: Vishal Raj
# Roll No.: 17EE35027
#

INF = 1e19
IN_OPEN = 0
IN_CLOSED = 1
UNDISCOVERED = 2

# A node of the graph
class Node:

    def __init__(self, id = -1):
        self.id = id
        self.adj_list = [] # [[child1,path_cost1],....]
        self.g_n = [] # [[parent1,cost1],[parent2,cost2],....]
        self.h_n = [0,0]
        self.f_n = []
        self.is_goal_node = False
        self.status = UNDISCOVERED

    def print_(self):
        print(self.id)
        print(self.adj_list)
        print(self.g_n)
        print(self.h_n)
        print(self.f_n)
        print(self.is_goal_node)
        print(self.status)

# Checks if v1 strictly dominates v2
def is_strictly_dominating(v1, v2):
    if v1[0] < v2[0] and v1[1] <= v2[1]:
        return True
    if v1[0] <= v2[0] and v1[1] < v2[1]:
        return True
    return False

# Returms sum of v1 and v2
def cost_sum(v1, v2):
    return [v1[0]+v2[0], v1[1]+v2[1]]

# Tells if current node is to be considered 
# based on existing cost of solution nodes in best_cost_list
def consider_this_node(node, best_cost_list):
    for e1 in best_cost_list:
        c1 = e1[1]
        tot = len(node.f_n)
        cnt = 0
        for e2 in node.f_n:
            c2 = e2[1]
            cnt += is_strictly_dominating(c1,c2)
        if cnt == tot:
            return False
    return True

# Returns F(n) given G(n) and H(n)
def get_f_n(g_n, h_n):
    return [[e[0],cost_sum(e[1],h_n)] for e in g_n]

# Returns non dominating subset of set cost_list
def get_non_dominating_subset(cost_list):
    ans = []
    if len(cost_list) == 1:
        return cost_list
    sz = len(cost_list)

    # Find the non-dominating subset
    for i in range(0,sz):
        consider_this = True
        c1 = cost_list[i][1]
        for j in range(0,sz):
            if i == j:
                continue
            c2 = cost_list[j][1]
            if is_strictly_dominating(c2,c1):
                consider_this = False
                break
        if consider_this:
            ans.append(cost_list[i])
    ans_fin = []
    sz = len(ans)

    # Remove duplicates
    for i in range(0,sz-1):
        cnt = 0
        a1 = ans[i]
        for j in range(i+1,sz):
            a2 = ans[j]
            if a1[0] == a2[0] and a1[1][0] == a2[1][0] and a1[1][1] == a2[1][1]:
                cnt += 1
        if cnt == 0:
            ans_fin.append(ans[i])
    ans_fin.append(ans[sz-1])
    return ans_fin

# Floyd Warshall Algorithm for calculating H(n)
def floyd_warshall(node_list,goal_nodes):
    num_nodes = len(node_list)
    h_n_list = [[0,0] for i in range(num_nodes)]
    dist = [[[INF,INF] for i in range(num_nodes)] for j in range(num_nodes)]
    for i in range(num_nodes):
        dist[i][i] = [0,0]
    for i in range(num_nodes):
        node = node_list[i]
        for e in node.adj_list:
            ch_id = e[0]
            path_cost = e[1]
            dist[i][ch_id] = path_cost

    # Floyd Warshal being run for the graph
    for k in range(num_nodes):
        for i in range(num_nodes):
            for j in range(num_nodes):
                for l in range(2):
                    dist[i][j][l] = min(dist[i][j][l],dist[i][k][l]+dist[k][j][l])

    # Calculated H(n) based on 'dist' matrix
    for i in range(num_nodes):
        node = node_list[i]
        if not node.is_goal_node:
            st = [INF,INF]
            for e in goal_nodes:
                e_id = ord(e)-ord('A')
                st = [min(st[0],dist[i][e_id][0]),min(st[1],dist[i][e_id][1])]
            h_n_list[i] = st

    return h_n_list

# Tells if n1 dominates n2
def dominates(n1, n2):
    for e in n2.f_n:
        for e1 in n1.f_n:
            if is_strictly_dominating(e1[1],e[1]):
                return True
    return False

# Prints the path from source to goal node
def print_paths(best_cost_list, node_list):
    for e in best_cost_list:
        print('^^^^^^^^^')
        path = []
        node_id = e[0]
        g_n = e[1]
        path.append(node_id)
        while node_id != -1:
            node = node_list[node_id]
            parent_ = -1
            for e1 in node.g_n:
                if g_n[0] == e1[1][0] and g_n[1] == e1[1][1]:
                    parent_ = e1[0]
                    break
            if parent_ != -1:
                path.append(parent_)
                for e1 in node_list[parent_].adj_list:
                    if e1[0] == node_id:
                        g_n = [g_n[0]-e1[1][0],g_n[1]-e1[1][1]]
                        break
                node_id = parent_
            else:
                break
        for i in range(len(path)):
            print(path[len(path)-i-1],)

edge_list = []
h_n_list = []
num_nodes = 0
goal_nodes = []


# Various graphs for testing
SELECT = 5

# Source: AIFA slides
if SELECT == 1:
    edge_list = [
        ['AC',[12,12]],
        ['AB',[10,10]],
        ['AD',[5,5]],
        ['BE',[8,8]],
        ['BF',[7,7]],
        ['BG',[6,6]],
        ['CF',[7,7]],
        ['CG',[5,5]],
        ['DG',[6,6]],
        ['DH',[4,4]],
        ['HC',[2,2]],
        ['EI',[18,18]],
        ['EJ',[14,14]],
        ['FI',[18,18]],
        ['FE',[9,9]],
        ['FG',[1,1]],
        ['GJ',[3,3]],
        ['GH',[13,13]]
    ]
    h_n_list = [[6,6],[5,5],[4,4],[3,3],[9,9],[3,3],[1,1],[3,3],[0,0],[0,0]]
    num_nodes = len(h_n_list)
    goal_nodes = ['I','J']

# Source: Self made
if SELECT == 2:
    edge_list = [
        ['AB',[1,1]],
        ['BD',[2,2]],
        ['AC',[2,2]],
        ['CD',[2,2]]
    ]
    h_n_list = [[0,0],[0,0],[0,0],[0,0]]
    num_nodes = len(h_n_list)
    goal_nodes = ['D']

# Source: Self made
if SELECT == 3:
    edge_list = [
        ['AB',[5,6]],
        ['BD',[3,4]],
        ['AC',[4,5]],
        ['CD',[3,4]]
    ]
    h_n_list = [[3,4],[2,3],[2,3],[0,0]]
    num_nodes = len(h_n_list)
    goal_nodes = ['D']

# Source: self made
if SELECT == 4:
    edge_list = [
        ['AB',[1,6]],
        ['BD',[1,4]],
        ['AC',[4,1]],
        ['CD',[3,1]]
    ]
    h_n_list = [[2,2],[1,1],[1,1],[0,0]]
    num_nodes = len(h_n_list)
    goal_nodes = ['D']

# Source: Computational Intelligence Book
if SELECT == 5:
    edge_list = [
        ['AB',[2,3]],
        ['AC',[4,2]],
        ['BD',[3,2]],
        ['BE',[2,3]],
        ['CE',[3,2]],
        ['CF',[3,2]],
        ['DG',[5,7]],
        ['EG',[1,3]],
        ['EH',[1,1]],
        ['FH',[5,2]],
        ['GI',[2,5]],
        ['GJ',[9,8]],
        ['HJ',[6,8]],
        ['HK',[4,3]]
    ]
    h_n_list = []
    num_nodes = 11
    goal_nodes = ['I','J','K']

# Applies Multiobjective DFBB algorithm on the given Graph
def multiobjective_dfbb(num_nodes, edge_list, goal_nodes, h_n_list = []):
    open_ = set() # OPEN list
    closed_ = set() # CLOSED list
    best_cost_list = [[-1,[INF,INF]]] # [[goal_node1,cost1],[goal_node2,cost2],....]
    node_list = []

    for i in range(num_nodes):
        new_node  = Node(i)
        node_list.append(new_node)

    # Convert edge list to adjacency list 
    for e in edge_list:
        i1 = ord(e[0][0])-ord('A')
        i2 = ord(e[0][1])-ord('A')
        node_list[i1].adj_list.append([i2,e[1]])

    start_node = node_list[0]
    start_node.g_n = [[-1,[0,0]]]
    start_node.f_n = get_f_n(start_node.g_n,start_node.h_n)
    open_.add(start_node.id) # Add start node to open list

    # Mark the goal nodes
    for e in goal_nodes:
        node_list[ord(e)-ord('A')].is_goal_node = True
    
    # Calculated H(n) if not already supplied
    if len(h_n_list) == 0:
        h_n_list = floyd_warshall(node_list,goal_nodes)
        print("h_n_list =",h_n_list)

    for i in range(num_nodes):
        node_list[i].h_n = h_n_list[i]

    for e in node_list:
        e.print_()
        print('-----------')

    # DFBB loop
    while len(open_):
        # POP one node to process it
        node_id = open_.pop()
        node_now = node_list[node_id]
        closed_.add(node_id)
        node_now.status = IN_CLOSED

        # consider or reject this node based on solutions obtained till now
        if not consider_this_node(node_now,best_cost_list):
            print(node_id)
            continue

        print('###########')
        print(node_id)

        # Iterate through all the children
        for e in node_now.adj_list:

            print(e)
            ch_id = e[0]
            ch_node = node_list[ch_id]
            path_cost = e[1]
            new_costs = []      
            for e1 in node_now.g_n:
                cost = e1[1]
                new_costs.append([node_id,cost_sum(cost,path_cost)])     
            for e1 in ch_node.g_n:
                new_costs.append(e1)    
            print(new_costs)
            ch_node.g_n = get_non_dominating_subset(new_costs)
            print(ch_node.g_n)
            ch_node.f_n = get_f_n(ch_node.g_n,ch_node.h_n)  
            print(ch_node.f_n)

            # consider or reject this node based on solutions obtained till now
            if not consider_this_node(ch_node,best_cost_list):
                ch_node.status = IN_CLOSED
                closed_.add(ch_id)
                if ch_id in open_:
                    open_.remove(ch_id)
                print('h1')
                continue

            # If accepted add it to OPEN list
            if ch_node.status == UNDISCOVERED:
                ch_node.status = IN_OPEN
                open_.add(ch_id)
            elif ch_node == IN_CLOSED:
                ch_node.status = IN_OPEN
                closed_.remove(ch_id)
                open_.add(ch_id)

            print(open_)
            print(closed_)

            # If node is a goal node then update the best cost list
            if ch_node.is_goal_node:
                new_best_cost_list = []
                for e1 in best_cost_list:
                    new_best_cost_list.append(e1)
                for e1 in ch_node.f_n:
                    new_best_cost_list.append([ch_id,e1[1]])
                print(new_best_cost_list)
                best_cost_list = get_non_dominating_subset(new_best_cost_list)

            print(best_cost_list)

    print('-------')
    print(best_cost_list)
    return best_cost_list, node_list 

# Applies Multiobjective A Star algorithm on the given Graph
def multiobjective_a_star(num_nodes, edge_list, goal_nodes, h_n_list = []):
    open_ = set() # OPEN Set
    closed_ = set() # CLOSED Set
    open_nd_ = set() # Non-dominating Subset of the Open set
    best_cost_list = [[-1,[INF,INF]]] # [[goal_node1,cost1],[goal_node2,cost2],....]
    node_list = []

    for i in range(num_nodes):
        new_node  = Node(i)
        node_list.append(new_node)

    # Convert edge list to adjacency list 
    for e in edge_list:
        i1 = ord(e[0][0])-ord('A')
        i2 = ord(e[0][1])-ord('A')
        node_list[i1].adj_list.append([i2,e[1]])

    start_node = node_list[0]
    start_node.g_n = [[-1,[0,0]]]
    start_node.f_n = get_f_n(start_node.g_n,start_node.h_n)
    open_.add(start_node.id) # Add start node to open list

    # Mark the goal nodes
    for e in goal_nodes:
        node_list[ord(e)-ord('A')].is_goal_node = True

    # Calculated H(n) if not already supplied
    if len(h_n_list) == 0:
        h_n_list = floyd_warshall(node_list,goal_nodes)
        print("h_n_list =",h_n_list)

    for i in range(num_nodes):
        node_list[i].h_n = h_n_list[i]

    for e in node_list:
        e.print_()
        print('-----------')

    # A Star loop
    while True:
        # Calculated the non-dominating subset OPEN_ND_ of OPEN set
        temp = []
        open_nd_.clear()
        for e in open_:
            if consider_this_node(node_list[e],best_cost_list):
                temp.append(e)
        for i in range(len(temp)):
            add_to_open_nd = True
            for j in range(len(temp)):
                if i != j:
                    if dominates(node_list[temp[j]],node_list[temp[i]]):
                        add_to_open_nd = False
                        break
            if add_to_open_nd:
                open_nd_.add(temp[i])

        # If no non-dominating subet of OPEN set 
        # exists then terminate
        if len(open_nd_) == 0:
            break

        # POP one node to process it
        print(open_nd_)
        node_id = open_nd_.pop()
        open_.remove(node_id)
        node_now = node_list[node_id]
        closed_.add(node_id)
        node_now.status = IN_CLOSED

        print('###########')
        print(node_id)

        # Iterate through all the children
        for e in node_now.adj_list:
            print(e)
            ch_id = e[0]
            ch_node = node_list[ch_id]
            path_cost = e[1]
            new_costs = []      
            for e1 in node_now.g_n:
                cost = e1[1]
                new_costs.append([node_id,cost_sum(cost,path_cost)])     
            for e1 in ch_node.g_n:
                new_costs.append(e1)    
            print(new_costs)
            ch_node.g_n = get_non_dominating_subset(new_costs)
            print(ch_node.g_n)
            ch_node.f_n = get_f_n(ch_node.g_n,ch_node.h_n)  
            print(ch_node.f_n)
            open_.add(ch_id)
            if ch_id in closed_:
                closed_.remove(ch_id)
            ch_node.status = IN_OPEN

            print(open_)
            print(closed_)

            # If node is a goal node then update the best cost list
            if ch_node.is_goal_node:
                new_best_cost_list = []
                for e1 in best_cost_list:
                    new_best_cost_list.append(e1)
                for e1 in ch_node.f_n:
                    new_best_cost_list.append([ch_id,e1[1]])
                print(new_best_cost_list)
                best_cost_list = get_non_dominating_subset(new_best_cost_list)
            print(best_cost_list)

    print('-------')
    print(best_cost_list)
    return best_cost_list, node_list 


# best_cost_list, node_list = multiobjective_dfbb(num_nodes, edge_list, goal_nodes, h_n_list = [])
best_cost_list, node_list = multiobjective_a_star(num_nodes, edge_list, goal_nodes, h_n_list = [])
print_paths(best_cost_list, node_list)