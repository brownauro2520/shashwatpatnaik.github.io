def calcneigh(E2N):
    node_neighbors = {}
    # iterate through the elements
    for element in E2N:
        # iterate through the nodes in the current element
        for node in element:
            # check if the node is already in the dictionary
            if node not in node_neighbors:
                node_neighbors[node] = []
            # add the two neighboring nodes to the set of neighbors
            for neighbors in element:
                if neighbors != node and not neighbors in node_neighbors[node]:
                    node_neighbors[node].append(neighbors)
    return node_neighbors

def smooth(neigh,N,C,w):
    # create an empty dictionary to store the node neighbors
    xi = C[N]

    eq1 = [(1-w)*x for x in xi]
    eq2 = w / len(neigh[N+1])
    temp = [0,0]
    for i in neigh[N+1]:
        temp = [x+y for x, y in zip(temp,C[i-1])]

    eq3 = [eq2*x for x in temp]
    return [x+y for x, y in zip(eq1,eq3)]

def smoothing(E2N,C,bcords,N,w):
    neigh = calcneigh(E2N) #smoothng
    for i in range(N):
        for j in range(len(C)):
            if j+1 not in bcords:
                C[j] = smooth(neigh,j,C,w)