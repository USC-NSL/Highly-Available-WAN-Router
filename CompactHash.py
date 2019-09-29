from __future__ import print_function
import numpy as np
import fractions as frlib
import time
import pickle
import gurobipy as gp

from LinkFailurePattern import load_RepLinkFP
from ExtremeTraffic import load_ETM
from Wiring import load_Wiring


HashSaveDir = "Hash-save"


def hashResultFilenameSingle(M, lkfp_id, algo, algoHash):
    fname = "Hash-M"
    for m in M:
        fname += "_" + str(m)
    fname += "-FP_" + str(lkfp_id)
    fname += "-" + algo + "-" + algoHash
    return fname


def save_HashResultSingle(result, M, lkfp_id, algo, algoHash):
    name = hashResultFilenameSingle(M, lkfp_id, algo, algoHash)
    fullpath = HashSaveDir + "/" + name
    fp = open(fullpath, "w")
    pickle.dump(result, fp)
    fp.close()
    print("Save result to", fullpath)

    
def load_HashResultSingle(M, lkfp_id, algo, algoHash):
    name = defaultHashResultFilenameSingle(M, lkfp_id, algo, algoHash)
    fullpath = HashSaveDir + "/" + name
    fp = open(fullpath, "r")
    hashResult = pickle.load(fp)
    fp.close()
    return hashResult


def getFGCD(f1, f2):
    dlcm = (f1.denominator*f2.denominator)/frlib.gcd(f1.denominator, f2.denominator)
    ngcd = frlib.gcd(f1*dlcm, f2*dlcm)
    return frlib.Fraction(ngcd, dlcm)
    

def fractionalize(radix, M, wiring):
    K = len(M)
    numEdge = radix

    # Determine fractional part and direction
    fracs = dict()
    for i in range(K):
        for j in range(K):
            if i == j:
                continue
            fracs[(i, j)] = dict()
            for e in range(numEdge):
                fin = frlib.Fraction(wiring[e][i])/frlib.Fraction(M[i])
                fout = frlib.Fraction(wiring[e][j])/frlib.Fraction(M[j])
                if fin > fout:
                    fval = fin - fout
                    fup = True
                elif fin < fout:
                    fval = fout - fin
                    fup = False
                else:
                    fval = 0
                    fup = None
                fracs[(i, j)][e] = { "val": fval, "up": fup }

    # Calculate FGCD of fractions
    fgcd = dict()
    for i in range(K):
        for j in range(K):
            if i == j:
                continue
            fgcd[(i, j)] = 0
            for e in range(numEdge):
                fgcd[(i, j)] = getFGCD(fgcd[(i, j)], fracs[(i, j)][e]["val"])

    # Calculate Flow counts
    flowCount = dict()
    upL1s = dict()
    downL1s = dict()
    for i in range(K):
        for j in range(K):
            if i == j:
                continue
            flowCount[(i, j)] = dict()
            upL1s[(i, j)] = set()
            downL1s[(i, j)] = set() 
            for e in range(numEdge):
                if fracs[(i, j)][e]["val"] == 0:
                    flowCount[(i, j)][e] = 0
                else:
                    flowCount[(i, j)][e] = fracs[(i, j)][e]["val"] / fgcd[(i, j)]
                assert flowCount[(i, j)][e] == int(flowCount[(i, j)][e])
                flowCount[(i, j)][e] = int(flowCount[(i, j)][e])
                if fracs[(i, j)][e]["up"] == True:
                    upL1s[(i, j)].add(e)
                elif fracs[(i, j)][e]["up"] == False:
                    downL1s[(i, j)].add(e)

    return (flowCount, fgcd, upL1s, downL1s)


def solveCompactWeight(radix, M, flowSet, lkfp, effCap, algoHash):
    numEdge = radix
    numCore = radix/2
    K = len(M)
    (flowCount, fgcd, upL1s, downL1s) = flowSet
    allUpL1s = set()
    allDownL1s = set()
    for i in range(K):
        for j in range(K):
            if i == j:
                continue
            allUpL1s.update(upL1s[(i, j)])
            allDownL1s.update(downL1s[(i, j)])    
    etms = load_ETM(M)
    if algoHash == "Optimal":
        pass
    elif algoHash == "Heuristic":
        utm = etms[0]
        for etm in etms:
            utm = np.maximum(utm, etm)
        etms = [utm]
    else:
        assert False, "Invalid algorithm"    

    model = gp.Model()

    x = gp.tupledict()
    for i in range(K):
        for j in range(K):
            if i == j:
                continue
            for e in upL1s[(i, j)]:
                for c in range(numCore):
                    en = "e"+str(e)
                    cn = "c"+str(c)
                    x[en, cn, i, j] = model.addVar(lb = 0,
                                                   vtype = gp.GRB.INTEGER,
                                                   name = "x_{0}_{1}_{2}_{3}".format(en, cn, i, j))
            for e in downL1s[(i, j)]:
                for c in range(numCore):
                    en = "e"+str(e)
                    cn = "c"+str(c)                
                    x[cn, en, i, j] = model.addVar(lb = 0,
                                                   vtype = gp.GRB.INTEGER,
                                                   name = "x_{0}_{1}_{2}_{3}".format(cn, en, i, j))

    xmax = model.addVar(lb = 0, vtype = gp.GRB.CONTINUOUS, name = "xmax")                    
    mu = model.addVar(lb = 1, vtype = gp.GRB.INTEGER, name = "mu")

    model.setObjective(xmax, sense=gp.GRB.MINIMIZE)

    # Auxiliary total number of entry at L2 switch
    for c in range(numCore):
        model.addConstr(x.sum("c"+str(c), "*", "*", "*") <= xmax, name = "Entry_c" + str(c))
    for e in range(numEdge):
        model.addConstr(x.sum("e"+str(e), "*", "*", "*") <= xmax, name = "Entry_e" + str(e))
        

    # Flow conservation at L2 switch
    for i in range(K):
        for j in range(K):
            if i == j:
                continue    
            for c in range(numCore):
                fin = x.sum("*", "c"+str(c), i, j)
                fout = x.sum("c"+str(c), "*", i, j)
                model.addConstr(fin == fout, name = "Match_c{0}_{1}_{2}".format(c, i, j))

    # Flow conservation at L1 switch
    for i in range(K):
        for j in range(K):
            if i == j:
                continue
            for e in upL1s[(i, j)]:
                fin = mu * flowCount[(i, j)][e]
                fout = x.sum("e"+str(e), "*", i, j)
                model.addConstr(fin == fout, name = "Match_+e{0}_{1}_{2}".format(e, i, j))
            for e in downL1s[(i, j)]:
                fin = x.sum("*", "e"+str(e), i, j)
                fout = mu * flowCount[(i, j)][e]
                model.addConstr(fin == fout, name = "Match_-e{0}_{1}_{2}".format(e, i, j))

    # Link capacity constraint L1 -> L2
    tmc = -1
    for etm in etms:
        tmc += 1
        for e in allUpL1s:
            for c in range(numCore):
                en = "e" + str(e)
                cn = "c" + str(c)
                rate = 0
                for (_, _, i, j) in x.subset(en, cn, "*", "*"):
                    rate += x[en, cn, i, j] * float(fgcd[(i, j)]) * effCap * etm[i, j]
                if (en, cn) in lkfp:
                    cap = 0
                else:
                    cap = 1
                model.addConstr(rate <= mu * cap, name = "Cap_{0}_{1}_{2}".format(en, cn, tmc))

    # Link capacity constraint L2 -> L1
    tmc = -1
    for etm in etms:
        tmc += 1
        for e in allDownL1s:
            for c in range(numCore):
                en = "e" + str(e)
                cn = "c" + str(c)
                rate = 0
                for (_, _, i, j) in x.subset(cn, en, "*", "*"):
                    rate += x[cn, en, i, j] * float(fgcd[(i, j)]) * effCap * etm[i, j]
                if (en, cn) in lkfp:
                    cap = 0
                else:
                    cap = 1
                model.addConstr(rate <= mu * cap, name = "Cap_{0}_{1}_{2}".format(cn, en, tmc))

    # model.write("test.lp")
    model.setParam("LogToConsole", 0)        
    model.optimize()
    if model.status == gp.GRB.Status.OPTIMAL:
        pass
    elif model.status == gp.GRB.Status.INFEASIBLE:
        print(M, effCap)
        assert False, "Model is infeasible."
    else:
        assert False, "Solver Status " + str(model.status)

    # Format flow result
    flowRoute = dict()
    for i in range(K):
        for j in range(K):
            if i == j:
                continue
            flowRoute[(i, j)] = dict()
            for e in upL1s[(i, j)]:
                for c in range(numCore):
                    en = "e" + str(e)
                    cn = "c" + str(c)
                    flowRoute[(i, j)][(en, cn)] = int(x[en, cn, i, j].X)
            for e in downL1s[(i, j)]:
                for c in range(numCore):
                    en = "e" + str(e)
                    cn = "c" + str(c)                    
                    flowRoute[(i, j)][(cn, en)] = int(x[cn, en, i, j].X)
            
    return flowRoute


def resolveWeight(radix, M, wiring, flowSet, flowRoute):
    numEdge = radix
    numCore = radix / 2
    K = len(M)
    (flowCount, fgcd, upL1s, downL1s) = flowSet    

    weights = dict()    
    # L2 hash
    for i in range(K):
        for j in range(K):
            if i == j:
                continue
            weights[(i, j)] = dict()
            for c in range(numCore):
                cn = "c" + str(c)
                weights[(i, j)][cn] = dict()
                gcd = 0
                for e in range(numEdge):
                    en = "e" + str(e)
                    if not flowRoute[(i, j)].has_key((cn, en)) or flowRoute[(i, j)][(cn, en)] == 0:
                        continue
                    gcd = frlib.gcd(gcd, flowRoute[(i, j)][(cn, en)])
                for e in range(numEdge):
                    en = "e" + str(e)
                    if not flowRoute[(i, j)].has_key((cn, en)) or flowRoute[(i, j)][(cn, en)] == 0:
                        continue
                    weights[(i, j)][cn][en] = flowRoute[(i, j)][(cn, en)] / gcd
    # L1 Hash
    for i in range(K):
        for j in range(K):
            if i == j:
                continue    
            for e in range(numEdge):
                en = "e" + str(e)
                # Augment trunk wires and internal links
                weights[(i, j)][en] = dict()                
                for t in range(wiring[e][j]):
                    weights[(i, j)][en]["t{0}_{1}".format(j, t)] = frlib.Fraction(1, M[j])
                for c in range(numCore):
                    cn = "c" + str(c)
                    if flowRoute[(i, j)].has_key((en, cn)):
                        weights[(i, j)][en][cn] = flowRoute[(i, j)][(en, cn)] * fgcd[(i, j)]
                # find FGCD and resolve weights
                tfgcd = 0
                for nh in weights[(i, j)][en].keys():
                    tfgcd = getFGCD(tfgcd, weights[(i, j)][en][nh])
                for nh in weights[(i, j)][en].keys():
                    if tfgcd == 0:
                        weights[(i, j)][en][nh] = int(weights[(i, j)][en][nh])
                    else:
                        assert weights[(i, j)][en][nh] / tfgcd == int(weights[(i, j)][en][nh] / tfgcd)
                        weights[(i, j)][en][nh] = int(weights[(i, j)][en][nh] / tfgcd)
    return weights


def computeCompactHash(radix, M, algo, algoHash, lkfp, effCap):
    K = len(M)
    wiring = load_Wiring(M, algo)
    flowSet = fractionalize(radix, M, wiring)
    flowRoute = solveCompactWeight(radix, M, flowSet, lkfp, effCap, algoHash)
    weights = resolveWeight(radix, M, wiring, flowSet, flowRoute)
    hashResult = (flowSet, flowRoute, weights)
    return hashResult


def test_LKFailure():
    radix = 4
    M = (2, 2, 4)
    K = len(M)
    algo = "Optimal"
    wiring = load_Wiring(M, algo)
    flowSet = fractionalize(radix, M, wiring)
    (flowCount, fgcd, upL1s, downL1s) = flowSet

    # Impose link failure pattern
    lkfp = (("e0", "c0"),)
    lkfp_id = 0
    effCap = 1
    algoHash = "Optimal"
    flowRoute = solveCompactWeight(radix, M, flowSet, lkfp, effCap, algoHash)
    weights = resolveWeight(radix, M, wiring, flowSet, flowRoute)

    print_hash_result(radix, M, wiring, flowSet, flowRoute, weights)
    result = (flowSet, flowRoute, weights)
    save_HashResultSingle(result, M, lkfp_id, algo, algoHash)


def test_L1Failure():
    radix = 4
    M = (2, 2, 4)
    K = len(M)
    algo = "Optimal"
    wiring = load_Wiring(M, algo)

    # Impose L1 switch failure
    MRes = (2, 1, 3)
    wiringRes = ((0, 0, 0), wiring[1], wiring[2], wiring[3])
    print("MRes:", MRes)
    print("wiringRes:", wiringRes)

    # Impose link failure
    lkfp = tuple()
    
    flowSet = fractionalize(radix, MRes, wiringRes)
    (flowCount, fgcd, upL1s, downL1s) = flowSet

    effCap = 0.5
    algoHash = "Optimal"
    flowRoute = solveCompactWeight(radix, M, flowSet, lkfp, effCap, algoHash)
    
    weights = resolveWeight(radix, MRes, wiringRes, flowSet, flowRoute)

    print_hash_result(radix, M, wiringRes, flowSet, flowRoute, weights)
    
    result = (flowSet, flowRoute, weights)



def print_hash_result(radix, M, wiring, flowSet, flowRoute, weights):
    numEdge = radix
    numCore = radix / 2
    K = len(M)
    (flowCount, fgcd, upL1s, downL1s) = flowSet    

    # Print hash weight per trunk pair
    for i in range(K):
        for j in range(K):
            if i == j:
                continue
            print("Trunk pair", (i, j))
            print("\tWiring", wiring)
            print("\tFGCD:", fgcd[(i, j)])
            print("\tFlowCount:", flowCount[(i, j)])
            print("\tUp:\t", upL1s[(i, j)])
            print("\tDown:\t", downL1s[(i, j)])
            print("\tRoute Up", end=": ")
            for e in upL1s[(i, j)]:
                en = "e" + str(e)                
                for c in range(numCore):
                    cn = "c" + str(c)
                    print("({0}, {1})".format(en, cn), flowRoute[(i, j)][(en, cn)], end=" ")
                print("\t", end=" ")
            print()

            print("\tRoute Down", end=": ")
            for c in range(numCore):
                cn = "c" + str(c)            
                for e in downL1s[(i, j)]:
                    en = "e" + str(e)
                    print("({0}, {1})".format(cn, en), flowRoute[(i, j)][(cn, en)], end=" ")
                print("\t", end=" ")
            print()

            print("\tWeight L1:")
            for e in range(numEdge):
                en = "e" + str(e)
                print("\t\t", en, weights[(i, j)][en])
            print("\tWeight L2:")
            for c in range(numCore):
                cn = "c" + str(c)            
                print("\t\t", cn, weights[(i, j)][cn])

    # Calculate statistic of WCMP entries
    maxWeight = 0
    print("Number of entries at L1:")
    for e in range(numEdge):
        en = "e" + str(e)
        cntWeight = 0
        for i in range(K):
            for j in range(K):
                if i == j:
                    continue
                cntWeight += sum(weights[(i, j)][en].values())
        maxWeight = max(maxWeight, cntWeight)
        print("\t", en, cntWeight)

    print("Number of entries at L2:")
    for c in range(numCore):
        cn = "c" + str(c)
        cntWeight = 0
        for i in range(K):
            for j in range(K):
                if i == j:
                    continue
                cntWeight += sum(weights[(i, j)][cn].values())
        maxWeight = max(maxWeight, cntWeight)
        print("\t", cn, cntWeight)

    print("Max number of entries", maxWeight)


if __name__ == "__main__":
    test_LKFailure()
    # test_L1Failure()


