from __future__ import print_function
import time
import pickle
import numpy as np

import gurobipy as gp

from ExtremeTraffic import load_ETM


WiringSaveDir = "Wiring-save"


def load_Wiring(M, algo):
    fullpath = WiringSaveDir + "/" + wiringFilename(M, algo)
    fp = open(fullpath, "r")
    wiring = pickle.load(fp)
    fp.close()
    return wiring


def save_Wiring(wiring, M, algo):
    fullpath = WiringSaveDir + "/" + wiringFilename(M, algo)
    fp = open(fullpath, "w")
    pickle.dump(wiring, fp)
    fp.close()
    print("Save result to", fullpath)    


def wiringFilename(M, algo):
    name = "Wiring-M"
    for m in M:
        name += "_" + str(m)
    return name + "-" + algo


def optimizeWiring(radix, M, algo):
    K = len(M)
    numEdge = radix
    numIntf = radix/2

    etms = load_ETM(M)

    model = gp.Model()

    assign = gp.tuplelist((e, k) for e in range(numEdge) for k in range(K))
    w = model.addVars(assign, lb=0, vtype=gp.GRB.INTEGER, name="w_")
    for (e, k) in assign:
        w[e, k].setAttr(gp.GRB.Attr.UB, min(M[k], numIntf))
    if algo == "Optimal":
        pass
    elif algo == "Heuristic":
        print("Heuristic speedup is enable.")
        utm = etms[0].copy()
        etms = [utm]
        for etm in etms:
            utm = np.maximum(utm, etm)        
        for (e, k) in assign:
            w[e, k].setAttr(gp.GRB.Attr.LB, np.floor(M[k]/float(numEdge)))
            w[e, k].setAttr(gp.GRB.Attr.UB, np.ceil(M[k]/float(numEdge)))
    else:
        assert False, "Unknown algorithm"

    aux = gp.tuplelist((i, j, e, t) for i in range(K) for j in range(K) if i != j for e in range(numEdge) for t in range(len(etms)))
    a = model.addVars(aux, lb=0, vtype=gp.GRB.CONTINUOUS, name="a_")
    b = model.addVar(lb=0, vtype=gp.GRB.CONTINUOUS, name="b")
    model.update()

    model.setObjective(b, gp.GRB.MINIMIZE)
    for k in range(K):
        model.addConstr(w.sum("*", k) == M[k], name="tcap_"+str(k))
    for e in range(numEdge):
        model.addConstr(w.sum(e, "*") == numIntf, name="icap_"+str(e))

    for (i, j, e, t) in aux:
        cname = "a_" + str(i) + "_" + str(j) + "_" + str(e) + "_" + str(t)
        etm = etms[t]
        model.addConstr(etm[i, j]*(w[e, i]/M[i] - w[e, j]/M[j]) <= a[i, j, e, t], name=cname)
    for t in range(len(etms)):
        model.addConstr(a.sum("*", "*", "*", t) <= b)

    print("Solver starts")
    model.optimize()
    
    status = model.status
    wiring = list()
    if status == gp.GRB.Status.OPTIMAL:
        print("Wiring is found.")
        for e in range(numEdge):
            we = [0]*K
            for k in range(K):
                we[k] = int(round(w[e, k].x))
            wiring.append(we)
    else:
        print("Fail to determine wiring.")
        assert False
    wiring.sort()
    return tuple(tuple(we) for we in wiring)


def getBaselineWiring(radix, M):
    K = len(M)
    numEdge = radix
    numIntf = radix/2
    wiring = list()
    for e in range(numEdge):
        wiring.append([0]*K)
    R = list(M)
    k = 0
    e = 0
    i = numIntf
    while True:
        while k < len(M) and R[k] == 0:
            k += 1
        if k == len(M):
            break
        if i > 0:
            wiring[e][k] += 1
            i -= 1
        else:
            e += 1
            wiring[e][k] += 1
            i = numIntf - 1
        R[k] -= 1            
    return tuple(tuple(w) for w in wiring)


def generateWiring(radix, M, algo):
    if algo in ("Optimal", "Heuristic"):
        wiring = optimizeWiring(radix, M, algo)
    elif algo == "Baseline":
        wiring = getBaselineWiring(radix, M)
    else:
        assert False, "unknown algorithm"

    save_Wiring(wiring, M, algo)
    

if __name__ == "__main__":
    radix = 4
    M = (2, 2, 4)
    algo = "Optimal"
    generateWiring(radix, M, algo)
