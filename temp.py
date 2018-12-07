import numpy as np
from gurobipy import *
import random
import time

N_CONSTANT = 12
P_CONSTANT = .5
ITERATIONS_CONSTANT = 1
flag_res = False
flag_unres = True
flag_turn = False


def connected(vertex,adjacency,n,visited):
    visited[vertex] = 1
    row = adjacency[vertex]
    
    for node in range(n):
        if row[node] > 0 and visited[node] < 1:
            visited = connected(node,adjacency,n,visited)
    
    return visited

def rand_DAG(N,p):
    go = True
    while go:
        dag = np.zeros([N,N])
        for i in np.arange(N):
            for j in np.arange(i+1,N):
                if (random.random() < p):
                    dag[i][j] = 1
        visited = np.zeros(N)
        visited = connected(0,dag,N,visited)
        if(visited[N-1] == 1):
            go = False
            for i in range(N-1):
                if visited[i] == 0:
                    dag[i] = np.zeros(N)
    for i in range(N):
        visited = np.zeros(N)
        visited = connected(i,dag,N,visited)
        if(visited[N-1] == 0):
            dag[:,i] = np.zeros(N)
    return dag

def solve_TOP(adjacency,n,k,prizes,numPaths,s):
    nodes = []
    for i in range(n):
        nodes.append(i)
    m = Model("OP")
    
    m.setParam(GRB.Param.PoolSearchMode,2)
    m.setParam(GRB.Param.PoolSolutions,numPaths)
    
    
    inflow = np.zeros(n)
    inflow[s] = k
    inflow[n-1] = -k
    
    visited = m.addVars(nodes, vtype=GRB.BINARY,  name="visited")
    flows = m.addVars(nodes,nodes, vtype = GRB.INTEGER, name="flows")
    
    # Flow conservation
    m.addConstrs(
    (flows.sum('*',j) + inflow[j] == flows.sum(j,'*') for j in nodes), "Conservation")
    
    # Flow only on existing arcs
    m.addConstrs((flows[i,j] <= k*adjacency[i,j] for i in nodes for j in nodes),"EdgesPresent")
    
    # Visited nodes
    m.addConstrs((visited[j] - quicksum(flows[i,j] for i in nodes) <= 0 for j in nodes), "WhetherVisited")
    
#     obj = quicksum(visited[j] * prizes[j])
    
    # Objective Function
    m.setObjective((quicksum(visited[j]*prizes[j] for j in nodes)), GRB.MAXIMIZE)
    
    m.optimize()
    if m.status == GRB.Status.OPTIMAL:
        print(m.status)
        solution = m.getAttr('x',flows)
        for i in nodes:
            for j in nodes:
                if solution[i,j] > 0:
                    print('%s -> %s: %g' % (i, j, solution[i,j]))
    
    #save model
    
    return m,flows

# Assumes k = 1, constructs path from solution
def construct_path(m,n,flows):
    solution = m.getAttr('x',flows)
    print(solution)
    path = [0]
    nodes = range(n)
    for i in nodes:
        for j in nodes:
            if solution[i,j] > 0.5:
                path.append(j)
    return path

def score_paths(paths,k,prizes):
    scores = np.zeros(k,dtype=int)
    prizes_copy = prizes.copy()
    print("Point 1")
    go = True
    step = -1
    while go:
        go = False
        step += 1
        for player in range(k):
            path = paths[player]
            if len(path) > step:
                scores[player] += prizes_copy[path[step]]
                prizes_copy[path[step]] = 0
                go = True
    return scores

def unreserved_best_response(adjacency,n,prizes,first_path):
    nodes = []
    for i in range(n):
        nodes.append(i)
    m = Model("OP2")
    
    inflow = np.zeros(n)
    inflow[0] = 1
    inflow[n-1] = -1
    
    visited = m.addVars(nodes, vtype=GRB.BINARY,  name="visited")
    flows = m.addVars(nodes,nodes, vtype = GRB.INTEGER, name="flows")
    arrivals2 = m.addVars(nodes, vtype=GRB.INTEGER, name="arrivals2")
#     arrivals1 = m.addVars(nodes, name="arrivals1")
    arrivals1 = m.addVars(nodes, lb = -1, name = "arrivals1")
    second_prizes = m.addVars(nodes, vtype = GRB.BINARY, name = "second_prizes")
    first_prizes = m.addVars(nodes, vtype = GRB.BINARY, name = "first_prizes")
    
    full_arrivals1 = np.zeros(n)
    #full_arrivals1[0] = 0
#     print(first_path)
    for i in range(len(first_path)):
#         print(i)
        temp = first_path[i]
#         print(temp)
        full_arrivals1[temp] = i -.5
    print(full_arrivals1)
    
        # Flow conservation
    m.addConstrs((flows.sum('*',j) + inflow[j] == flows.sum(j,'*') for j in nodes), "Conservation")
    
    m.addConstrs((arrivals1[i] == full_arrivals1[i] for i in nodes),"ReallyStupid")
	
    #m.addConstrs((arrivals1[i] <= full_arrivals1[i]+. for i in nodes),"ReallyStupider")
	
    # Flow only on existing arcs
    m.addConstrs((flows[i,j] <= adjacency[i,j] for i in nodes for j in nodes),"EdgesPresent")
    
    # Prizes uncollected if not visited
    m.addConstrs((second_prizes[j] <= quicksum(flows[i,j] for i in nodes) for j in nodes), "WhetherVisited")
    
    # Prizes uncollected if not visited FIRST
    m.addConstrs((second_prizes[i] >= (arrivals1[i] - arrivals2[i])/n for i in nodes), "visitedFirst")
	
    m.addConstrs((first_prizes[i] >= (arrivals2[i] - arrivals1[i])/n for i in nodes), "visitedSecond")
	
    m.addConstrs((first_prizes[i] + second_prizes[i] <= 1 for i in nodes), "SingleDip")
	
    
    # Enforce correct arrival times
    #m.addConstrs((arrivals2[i] >= (arrivals2[j] + 1) * flows(j,i) for i >=1 in nodes for j < i in nodes), "arrivalEnforced")
    m.addConstrs((arrivals2[i] >= arrivals2[j] + (flows[j,i] - 1)* n + flows[j,i] for i in nodes for j in nodes), "arrivalEnforced")
#     m.addConstrs((arrivals2[i] >= arrivals2[j] -1 + 2*flows[j,i] for i in nodes for j in nodes), "arrivalEnforced")
    
    # Set starting point
    m.addConstr(arrivals2[0] == 0, "startEnforced")
    
    m.setObjective((quicksum(second_prizes[j] for j in nodes)), GRB.MAXIMIZE)
    
    m.optimize()
    print(m.objval)
    if m.status == GRB.Status.OPTIMAL:
        print(m.status)
        solution = m.getAttr('x',flows)
        for i in nodes:
            print("second prize",i,second_prizes[i].x)
            print("first_prize",i,first_prizes[i].x)
            print("first arrival",arrivals1[i].x)
            print("second arrival",arrivals2[i].x)
            print(arrivals1[i].x - arrivals2[i].x)
            for j in nodes:
                if solution[i,j] > 0:
                    print('%s -> %s: %g' % (i, j, solution[i,j]))
                    print(arrivals2[j].x)
                    print(arrivals2[i].x + (flows[i,j].x -1)*n + flows[i,j].x)
    
    #save model
    print(first_path)
    print(construct_path(m,n,flows))
	
	
    return m,flows
    
# Analysis for k best paths (2 players)
# def k_paths(adj,n,prizes,num_best):
#     model, flows = solve_TOP


def unreserved_solver(adjacency, n, prizes, k):
    model1, flows1 = solve_TOP(adjacency,n,1,prizes,1,0);
    print("point 1")
    path1 = construct_path(model1,n,flows1)
    solution = 0
    first_best_score = 0
    first_best_path = []
    second_best_score = 0
    second_best_path = []
    go = True
    while go:
        print("point 2")
        path1 = construct_path(model1,n,flows1)
        print(path1)
        path_val = model1.objval
        model2, flows2 = unreserved_best_response(adjacency,n,prizes,path1)
        print("\n\n Check here \n\n")
        print(model2.objval)
#         print(model2.getAttr('x'))
        path2 = construct_path(model2,n,flows2);
        print(path2)
        scores = score_paths([path1,path2],2,prizes)
#         print(scores)
        if scores[0] > first_best_score:
            first_best_score = scores[0]
            first_best_path = path1
            second_best_score = scores[1]
            second_best_path = path2
        if first_best_score >= path_val:
            go = False
        else:
            l = len(path1)-1
            model1.addConstr((quicksum(flows1[path1[j],path1[j+1]] for j in range(l)) <= l-2))
            model1.optimize()
    return [first_best_score,second_best_score],[first_best_path,second_best_path]

class State:
    def __init__(self,node1,node2,prizes_between):
        self.n1 = node1
        self.n2 = node2
        self.prizes = prizes_between
    
state_space = {}

def construct_turn_paths(state_space,initial_state,n):
    def recursive_path_constructor(state,n):
        nonlocal turn
        nonlocal path1
        nonlocal path2
        next_state = state_space[state]
        if state.n1 < n-1:
            if turn == 0:
                path1.append(next_state['next_step'])
            else:
                path2.append(next_state['next_step'])
        turn = (turn + 1) % 2
        if state.n2 < n-1:
            print(next_state)
            recursive_path_constructor(next_state['next_state'],n)
    
    turn = 0
    path1 = [0]
    path2 = [0]
    recursive_path_constructor(initial_state,n)
    return path1,path2

def turn_wrapper(adj,n,prizes):
    def turn_recurser(adj,n,prizes,state):
        nonlocal state_space
        next_step = -1
        value = 0
        n1 = state.n1
        n2 = state.n2
        if state in state_space:
            return state_space[state]
        else:
            # Check if next player is at the final node
            if n2 == n-1:
                # Check if current player is at the final node
                if n1 == n-1:
                    state_space[state] = {'val1':0,'val2':0,'next_state':state,'next_step': n-1}

                # Current player remains
                else:
                    new_prizes = prizes
                    n1 = state.n1
                    for i in range(n1+1,n):
                        new_prizes[i] = state.prizes[i-n1-1]
                    model,flows = solve_TOP(adj,n,1,new_prizes,1,n1)
                    state_space[state] = {'val1':model.objval,'val2':0, 'next_state':State(n-1,n-1,[]),'next_step': construct_path(model,n,flows)}

            # Check if current player is at final node
            elif n1 == n-1:
                new_state = State(state.n2,n-1,state.prizes)
                turn_recurser(adj,n,prizes,new_state)
                state_space[state] = {'val1':0, 'val2': state_space[new_state]['val1'],'next_state':new_state,'next_step':n-1}

            # If neither player is at final node
            else:
                row = adj[n1]
                next_moves = []
                next_step = n-1
                for i in range(n):
                    if row[i] == 1:
                        next_moves.append(i)

                # Assumes non-negative prizes
                state_space[state] = {'val1':-1,'val2':-1,'next_state':State(-1,-1,[]),'next_step':n-1}
                for new_step in next_moves:
                    step_prize = 0
                    new_prizes = []
                    if new_step == n2:
                        new_prizes = []
                    elif new_step > n2:
                        if n2 >= n1:
                            print([n1,n2])
                            new_prizes = prizes[n2+1:new_step+1]
                            step_prize = prizes[new_step]
                            print(new_prizes)
                            new_prizes[new_step-(n2+1)] = 0
                            print(new_prizes)
                        else: 
                            temp = prizes[n1+1:new_step+1]
                            step_prize = prizes[new_step]
                            temp[new_step-(n1+1)] = 0
                            new_prizes = state.prizes + temp
                    else:
                        new_prizes = state.prizes[new_step-(n1+1):]
                        step_prize = state.prizes[new_step-n1]
                    new_state = State(state.n2,new_step,new_prizes)
                    turn_recurser(adj,n,prizes,new_state)
                    if step_prize+state_space[new_state]['val2'] > state_space[state]['val1']:
                        state_space[state] = {'val1':step_prize+state_space[new_state]['val2'],
                                              'val2':state_space[new_state]['val1'], 'next_state':new_state,'next_step':new_step}
          
    state_space = {}
    state = State(0,0,[])
    turn_recurser(adj,n,prizes,state)
    path1,path2 = construct_turn_paths(state_space,state,n)
    p1 = []
    p2 = []
    for i in path1:
        if type(i) is list:
            for j in i:
                p1.append(j)
        else:
            p1.append(i)
    for i in path2:
        if type(i) is list:
            for j in i:
                p2.append(j)
        else:
            p2.append(i)
    scores = score_paths([p1,p2],2,prizes)
    return scores,[p1,p2]

def reserved_path_solver(adjacency,n,k,prizes):
    paths = []
    scores = []
    new_prizes =prizes.copy()
    for i in range(k):
        model, flows = solve_TOP(adjacency,n,1,new_prizes,1,0)
        path = construct_path(model,n,flows)
        score = score_paths([path],1,new_prizes)
        print(path,score)
        for j in path:
            new_prizes[j] = 0
        paths.append(path)
        scores.append(score[0])
    return scores,paths

def effective_n(adjacency, n):
    num = 1
    for i in range(n):
        row = adjacency[i]
        if np.sum(row) >=1:
            num +=1
    return num

def run_trials(n,p,iterations):

	central = []
	reserved = []
	turn = []
	unreserved = []
	nodes = []
	comp_central = []
	comp_reserved = []
	comp_unreserved = []
	comp_turn = []
	start_time = time.time()
	elapsed_time = time.time() - start_time
	for i in range(iterations):
		#n = random.randint(5,25)
		#p = .1 + .6*random.random()
		k = 2
		prizes = []
		for i in range(n):
			prizes.append(random.randint(1,1))
		prizes[0] = 0
		prizes[n-1] = 0
		adjacency = rand_DAG(n,p)
		n_prime = effective_n(adjacency,n)
		nodes.append(n_prime)
		print(n_prime)
		print(p)
		print(adjacency)
		
		# Compute central solution
		start_time = time.time()
		m1, flows1 = solve_TOP(adjacency,n,k,prizes,1,0)
		comp_central.append(time.time()-start_time)
		central.append(m1.objval)
		holder = 0;
		
		print("\n\n Finished Central \n\n")
		print(prizes)
		
		# Compute reserved solution
		if flag_res:
			start_time = time.time()
			scores, paths = reserved_path_solver(adjacency,n,k,prizes)
			comp_reserved.append(time.time()-start_time)
			reserved.append(np.sum(scores))
		
		# Compute unreserved solution
		if flag_unres:
			start_time = time.time()
			scores, holder = unreserved_solver(adjacency,n,prizes,k)
			comp_unreserved.append(time.time()-start_time)
			unreserved.append(np.sum(scores))
		
		# Compute turn-based solution
		if flag_turn:
			start_time = time.time()
			scores, paths = turn_wrapper(adjacency,n,prizes)
			comp_turn.append(time.time()-start_time)
			turn.append(np.sum(scores))


	print(adjacency)
	print(holder)
	print(prizes)
	print(central,reserved,unreserved,turn)
	print(comp_central,comp_reserved,comp_unreserved,comp_turn)
	
if __name__ == "__main__":
	run_trials(N_CONSTANT,P_CONSTANT,ITERATIONS_CONSTANT)