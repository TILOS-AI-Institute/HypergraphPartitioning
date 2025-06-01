import argparse
argparser = argparse.ArgumentParser()
argparser.add_argument("hgr")
argparser.add_argument("num_parts", type=int)
argparser.add_argument("ub_factor", type=float)
args = argparser.parse_args()

path_hgr = args.hgr
num_parts = args.num_parts
ub_factor = args.ub_factor / 100

import numpy as np

with open(path_hgr, "r") as f:
    for idx_line, line in enumerate(f):
        words = list(map(int, line.strip().split(" ")))
        if idx_line == 0:
            nE, nV = words[:2]
            mode = 0 if (len(words) == 2) else words[2]
            q, r = divmod(mode, 10)
            is_e_weighted = (r == 1)
            is_v_weighted = (q == 1)
            wvs = np.ones(nV)
            wes = np.ones(nE)
            es = []
            continue
        if idx_line < (nE + 1):
            eIdx = idx_line - 1
            if is_e_weighted:
                wes[eIdx] = words[0]
                e = words[1:]
            else:
                e = words
            e = np.array(e, dtype=int) - np.ones(len(e), dtype=int)
            es.append(e)
            continue
        vIdx = idx_line - (nE + 1)
        wvs[vIdx] = words[0]

import gurobipy as gp
from gurobipy import GRB

m = gp.Model()
X = np.array([[
    m.addVar(vtype=GRB.BINARY, name=f"X_{vIdx}_{i}")
    for i in range(num_parts)]
    for vIdx in range(nV)], dtype=object)
Y = np.array([[
    m.addVar(vtype=GRB.BINARY, name=f"Y_{eIdx}_{i}")
    for i in range(num_parts)]
    for eIdx in range(nE)], dtype=object)
W = np.sum(wvs)
Wmin = W * (1/num_parts - ub_factor)
Wmax = W * (1/num_parts + ub_factor)
for i in range(num_parts):
    sumwx = gp.LinExpr(wvs, X[:, i])
    m.addConstr(sumwx >= Wmin)
    m.addConstr(sumwx <= Wmax)
for vIdx in range(nV):
    m.addConstr(gp.LinExpr(np.ones(num_parts), X[vIdx, :]) == 1)
for i in range(num_parts):
    for eIdx, e in enumerate(es):
        for vIdx in e:
            m.addConstr(Y[eIdx, i] <= X[vIdx, i])
m.setObjective(gp.LinExpr(
    (-wes[eIdx], Y[eIdx, i])
    for i in range(num_parts)
    for eIdx in range(nE)
))
status = m.optimize()
phi = np.zeros(nV, dtype=int)
for vIdx in range(nV):
    for i in range(num_parts):
        if X[vIdx, i].X == 1:
            phi[vIdx] = i
            break

path_sol = f"{path_hgr}.part.{num_parts}"
with open(path_sol, "w") as f:
    for vIdx in range(nV):
        f.write(f"{phi[vIdx]}\n")
