import numpy as np
from gurobipy import *
import random
import time
import csv
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

N_CONSTANT = 45
P_CONSTANT = .25
ITERATIONS_CONSTANT = 8
PRIZE_MAX_CONSTANT = 5
flag_res = True
flag_unres = True
flag_turn = True

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

def rand_Delaunay_graph(n):
    x = np.zeros(n)
    y = np.zeros(n)
    for i in range(n):
        x[i] = random.random()
        y[i] = random.random()
    x = np.sort(x)
    points = np.zeros([n,2])
    points[:,0] = x
    points[:,1] = y
    tri = Delaunay(points)
    # plt.triplot(points[:,0], points[:,1], tri.simplices.copy())
    # plt.plot(points[:,0], points[:,1], 'o')
    # plt.show()
    adj = np.zeros([n,n])
    for item in tri.vertices:
        adj[min(item[0],item[1]),max(item[0],item[1])] = 1
        adj[min(item[0],item[2]),max(item[0],item[2])] = 1
        adj[min(item[1],item[2]),max(item[1],item[2])] = 1
    return adj


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
    # print(solution)
    path = [0]
    nodes = range(n)
    for i in nodes:
        for j in nodes:
            if solution[i,j] > 0.5:
                path.append(j)
    return path

def score_paths(paths,k,prizes):
    scores = np.zeros(k,dtype=float)
    prizes_copy = prizes.copy()
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

    visited2 = m.addVars(nodes, vtype=GRB.BINARY,  name="visited1")
    flows = m.addVars(nodes,nodes, vtype = GRB.INTEGER, name="flows")
    arrivals2 = m.addVars(nodes, vtype=GRB.INTEGER, name="arrivals2")
#     arrivals1 = m.addVars(nodes, name="arrivals1")
    arrivals1 = m.addVars(nodes, lb = -1, name = "arrivals1")
    second_prizes = m.addVars(nodes, vtype = GRB.BINARY, name = "second_prizes")
    first_prizes = m.addVars(nodes, vtype = GRB.BINARY, name = "first_prizes")

    arrivals1 = np.zeros(n)
    #full_arrivals1[0] = 0
#     print(first_path)
    visited1 = np.zeros(n)
    for i in range(len(first_path)):
#         print(i)
        temp = first_path[i]
#         print(temp)
        arrivals1[temp] = i -.5
        visited1[temp] = 1
    # print(full_arrivals1)

        # Flow conservation
    m.addConstrs((flows.sum('*',j) + inflow[j] == flows.sum(j,'*') for j in nodes), "Conservation")

    #m.addConstrs((arrivals1[i] <= full_arrivals1[i]+. for i in nodes),"ReallyStupider")

    # Flow only on existing arcs
    m.addConstrs((flows[i,j] <= adjacency[i,j] for i in nodes for j in nodes),"EdgesPresent")

    # Prizes uncollected if not visited
    m.addConstrs((visited2[j] <= quicksum(flows[i,j] for i in nodes) for j in nodes), "WhetherVisited1")

    m.addConstrs((visited2[j] >= 1/n* quicksum(flows[i,j] for i in nodes) for j in nodes), "WhetherVisited2")

    m.addConstrs((visited2[j] >= second_prizes[j] for j in nodes), "VisitToCollect")

    # Prizes uncollected if not visited FIRST
    m.addConstrs((second_prizes[i] >= (visited2[i]*arrivals1[i] - arrivals2[i])/n for i in nodes), "visitedFirst")

    m.addConstrs((first_prizes[i] >= (visited1[i]*arrivals2[i] - arrivals1[i])/n for i in nodes), "visitedSecond")

    m.addConstrs((first_prizes[i] + second_prizes[i] <= 1 for i in nodes), "SingleDip")


    # Enforce correct arrival times
    m.addConstrs((arrivals2[i] >= arrivals2[j] + (flows[j,i] - 1)* n + flows[j,i] for i in nodes for j in nodes), "arrivalEnforced")

    # Set starting point
    m.addConstr(arrivals2[0] == 0, "startEnforced")

    m.setObjective((quicksum(second_prizes[j] for j in nodes)), GRB.MAXIMIZE)

    m.optimize()
    print(m.objval)
    if m.status == GRB.Status.OPTIMAL:
        print(m.status)
        solution = m.getAttr('x',flows)
        for i in nodes:
            # print("second prize",i,second_prizes[i].x)
            # print("first_prize",i,first_prizes[i].x)
            # print("first arrival",arrivals1[i])
            # print("second arrival",arrivals2[i].x)
            # print(arrivals1[i] - arrivals2[i].x)
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
    def __init__(self,node1,node2,visited_nodes):
        self.n1 = node1
        self.n2 = node2
        self.visited = visited_nodes

state_space = {}

def getKey(state):
    return (state.n1,state.n2, tuple(state.visited))

def construct_turn_paths(state_space,initial_state,n):
    def recursive_path_constructor(state,n):
        nonlocal turn
        nonlocal path1
        nonlocal path2
        next_state = state_space[getKey(state)]
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
    def turn_recurser(state):
        nonlocal state_space
        nonlocal adj
        nonlocal n
        nonlocal prizes

        value = 0
        n1 = state.n1
        n2 = state.n2
        key = getKey(state)
        if key in state_space:
            return
        elif n2 == n-1:
            current_prizes = prizes.copy()
            current_visited = []
            for i in state.visited:
                if i > n1:
                    current_visited.append(i)
            current_visited.append(n1)
            current_visited = list(set(current_visited))
            for i in current_visited:
                current_prizes[i] = 0

            # print("\nIn TOP Solver\n")
            # print("N1: ",n1)
            # print("N2: ",n2)
            # # print("New step: ",new_step)
            # # print("Prize at next step: ",step_prize)
            # print(current_prizes,current_visited)
            # print("\n\n")
            model,flows = solve_TOP(adj,n,1,current_prizes,1,n1)
            state_space[key] = {'val1':model.objval,'val2':0, 'next_state':State(n-1,n-1,[]),'next_step': construct_path(model,n,flows)}
            print(model.objval,"\n\n")
        else:
            row = adj[n1]
            next_moves = []
            for i in range(n):
                if row[i] == 1:
                    next_moves.append(i)
            # Assumes non-negative prizes
            state_space[key] = {'val1':-1,'val2':-1,'next_state':State(-1,-1,[]),'next_step':n-1}
            for new_step in next_moves:
                current_prizes = prizes.copy()
                current_visited = []
                lb = min(new_step,n2)
                step_prize = 0
                if n1 >= n2:
                    current_visited = state.visited.copy()
                    current_visited.append(new_step)
                    step_prize = prizes[new_step]
                elif new_step > n2:
                    current_visited = [new_step]
                    step_prize = prizes[new_step]
                else:
                    current_visited.append(n2)
                    for i in state.visited:
                        if i >= lb:
                            current_visited.append(i)
                    if new_step not in state.visited:
                        step_prize = prizes[new_step]

                current_visited = list(set(current_visited))
                for i in current_visited:
                    current_prizes[i] = 0
                current_prizes[new_step] = 0
                new_state = State(n2,new_step,current_visited)
                new_key = getKey(new_state)
                print("\n\n")
                print("N1: ",n1)
                print("N2: ",n2)
                print("New step: ",new_step)
                print("Prize at next step: ",step_prize)
                print(current_prizes,current_visited)

                print("\n\n")
                turn_recurser(new_state)
                if step_prize+state_space[new_key]['val2'] > state_space[key]['val1']:
                    state_space[key] = {'val1':step_prize+state_space[new_key]['val2'],
                                          'val2':state_space[new_key]['val1'], 'next_state':new_state,'next_step':new_step}

    state_space = {}
    state = State(0,0,[])
    turn_recurser(state)
    path1,path2 = construct_turn_paths(state_space,state,n)
    p1 = []
    p2 = []
    for i in path1:
        if type(i) is list:
            i.pop(0)
            for j in i:
                p1.append(j)
        else:
            p1.append(i)
    for i in path2:
        if type(i) is list:
            i.pop(0)
            for j in i:
                p2.append(j)
        else:
            p2.append(i)
    scores = score_paths([p1,p2],2,prizes)
    print(scores,[p1,p2])
    return scores,[p1,p2]

def reserved_path_solver(adjacency,n,k,prizes):
    paths = []
    scores = []
    new_prizes =prizes.copy()
    for i in range(k):
        print("\n\n",new_prizes,"\n\n")
        model, flows = solve_TOP(adjacency,n,1,new_prizes,1,0)
        path = construct_path(model,n,flows)
        score = score_paths([path],1,new_prizes)
        print(path,score)
        for j in path:
            new_prizes[j] = 0

        paths.append(path)
        print("\nScore: ",score[0],"\n")
        scores.append(score[0])
    return scores,paths

def effective_n(adjacency, n):
    num = 1
    for i in range(n):
        row = adjacency[i]
        if np.sum(row) >=1:
            num +=1
    return num

def run_test():
    adj = [[0,1,1,0,0,0],
            [0,0,0,1,1,0],
            [0,0,0,0,1,0],
            [0,0,0,0,0,1],
            [0,0,0,0,0,1],
            [0,0,0,0,0,0]]
    adjacency = np.zeros([6,6])
    for i in range(6):
        for j in range(6):
            adjacency[i,j] = adj[i][j]
    prizes = [0,1.1,1.0,1.0,1.1,0]
    k = 2
    n = 6
    scores_correct, paths = solve_TOP(adjacency,n,k,prizes,1,0)
    scores_res, paths = reserved_path_solver(adjacency,n,k,prizes)
    scores_unres, holder = unreserved_solver(adjacency,n,prizes,k)
    scores_turn, paths = turn_wrapper(adjacency,n,prizes)
    correct = [2.2,1]                           # All strategies should result in this
    print("TOP Scores: ", scores_correct.objval)       # Should be [2.1,2.1]
    print("Correct game scores: ",correct)
    print("Reserved scores: ",scores_res)
    print("Unreserved scores: ",scores_unres)
    print("Turn scores: ",scores_turn)

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
    # filename = ('ER_trial_results_{}_{}.csv'.format(n,p))
    filename = ('Delaunay_trial_results_{}.csv'.format(n))
    for i in range(iterations):
        #n = random.randint(5,25)
        #p = .1 + .6*random.random()
        k = 2
        prizes = []
        for i in range(n):
            prizes.append(random.randint(1,PRIZE_MAX_CONSTANT))
        prizes[0] = 0
        prizes[n-1] = 0
        # adjacency = rand_DAG(n,p)
        adjacency = rand_Delaunay_graph(n)
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
        with open(filename, mode='w') as trial_results:
            results = csv.writer(trial_results, delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)
            results.writerow(nodes)
            results.writerow(central)
            results.writerow(reserved)
            results.writerow(unreserved)
            results.writerow(turn)
            results.writerow(comp_central)
            results.writerow(comp_reserved)
            results.writerow(comp_unreserved)
            results.writerow(comp_turn)

    print(adjacency)
    print(holder)
    print(prizes)
    print(central,reserved,unreserved,turn)
    print(comp_central,comp_reserved,comp_unreserved,comp_turn)
    with open(filename, mode='w') as trial_results:
        results = csv.writer(trial_results, delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)
        results.writerow(nodes)
        results.writerow(central)
        results.writerow(reserved)
        results.writerow(unreserved)
        results.writerow(turn)
        results.writerow(comp_central)
        results.writerow(comp_reserved)
        results.writerow(comp_unreserved)
        results.writerow(comp_turn)


if __name__ == "__main__":
    run_trials(N_CONSTANT,P_CONSTANT,ITERATIONS_CONSTANT)
    # run_test()
