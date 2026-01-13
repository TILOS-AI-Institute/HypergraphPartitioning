import os


"""
Evaluator:  Calculate cutsize and blocks_balance for a partitioning solution
"""
def Evaluator(hypergraph_file, solution_file, Nparts, UBfactor):
    if(os.path.exists(solution_file)):
        pass
    else:
        return 1e9, None

    with open(hypergraph_file) as f:
        content = f.read().splitlines()
    f.close()

    if(len(content) == 0):
        return 1e9, None

    items = content[0].split()
    num_hyperedges = int(items[0])
    num_vertices = int(items[1])
    flag = 0
    if(len(items) == 3):
        flag = int(items[2])

    hyperedges = []
    vertices_weight = [1 for i in range(num_vertices)]
    hyperedges_weight = [1 for i in range(num_hyperedges)]

    if((flag % 10) == 1):
        for i in range(1, num_hyperedges + 1):
            items = [int(item) for item in content[i].split()]
            hyperedges_weight[i - 1] = items[0]
            hyperedge = [items[j] - 1 for j in range(1, len(items))]
            hyperedges.append(hyperedge)
    else:
        for i in range(1, num_hyperedges + 1):
            items = [int(item) - 1 for item in content[i].split()]
            hyperedges.append(items)

    if(flag >= 10):
        for i in range(num_vertices):
            idx = i + num_hyperedges + 1
            vertices_weight[i] = int(content[idx])

    part = [-1 for i in range(num_vertices)]
    with open(solution_file) as f:
        content = f.read().splitlines()
    f.close()


    if(len(content) != num_vertices):
        return 1e9, None


    blocks_balance = [0.0 for i in range(Nparts)]
    total_weight = 0.0
    for i in range(len(content)):
        part_id = int(content[i])
        if(part_id == -1):
            return 1e9, None
        part[i] = part_id
        blocks_balance[part_id] += vertices_weight[i]
        total_weight += vertices_weight[i]

    max_balance = (100.0 / Nparts  + UBfactor) * 0.01;
    max_balance = max_balance * total_weight
    min_balance = (100.0 / Nparts  - UBfactor) * 0.01;
    min_balance = min_balance * total_weight

    flag = True
    for block_balance in blocks_balance:
        if (block_balance > max_balance):
            flag = False
        if (block_balance < min_balance):
            flag = False

    for i in range(len(blocks_balance)):
        blocks_balance[i] = blocks_balance[i] / total_weight

    num_cut = 0
    for i in range(num_hyperedges):
        part_list = []
        for vertex in hyperedges[i]:
            if(part[vertex] not in part_list):
                part_list.append(part[vertex])

        if(len(part_list) > 1):
            num_cut += hyperedges_weight[i]

    if (flag == False):
        print("Error !!!")
        print("cutsize   :   ", num_cut)
        print("blocks_balance :  ",  blocks_balance)
        return 1e9, blocks_balance, num_vertices, num_hyperedges

    return num_cut, blocks_balance, num_vertices, num_hyperedges
