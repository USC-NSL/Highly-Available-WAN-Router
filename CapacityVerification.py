from __future__ import print_function
import time
import numpy as np
import gurobipy as gp

from Router import Router
from ExtremeTraffic import load_ETM
from Wiring import load_Wiring
from LinkFailurePattern import load_RepLinkFP
from L1FailurePattern import getL1FP



def calculateEffectiveCapacityRouter(router, wiring, tmat):
    numEdge = router.numEdge
    numCore = router.numCore
    numTrunk = len(wiring[0])
    numTrunkWire = [ sum( [ wiring[ei][ti] for ei in range(numEdge) ] ) for ti in range(numTrunk) ]

    edgeNodes = set(["e"+str(ei) for ei in range(numEdge)])
    edgeResNodes = edgeNodes.difference(router.failedL1s)
    coreNodes = set(["c"+str(ci) for ci in range(numCore)])
    unLinks = set([(en, cn) for en in edgeNodes for cn in coreNodes])
    diLinks = set(list(unLinks) + [(l[1], l[0]) for l in unLinks])
    trunks = set(["t"+str(ti) for ti in range(numTrunk)])
    trunkPairs = set([(tn, tm) for tn in trunks for tm in trunks])

    # formulate optimization
    model = gp.Model()
    
    xindex = gp.tuplelist((ni, nj, ti, tj) for (ni, nj) in diLinks for (ti, tj) in trunkPairs if ti != tj)
    x = model.addVars(xindex, lb=0, ub=1, vtype=gp.GRB.CONTINUOUS, name="x_")

    beta = model.addVar(lb=0, ub=1, vtype=gp.GRB.CONTINUOUS, name="beta")
    
    model.update()
    model.setParam("LogToConsole", 0)
    model.setObjective(beta, gp.GRB.MAXIMIZE)

    # Flow conservation at core switches
    c_cons = gp.tupledict()
    for cn in coreNodes:
        for (ti, tj) in trunkPairs:
            if ti == tj:
                continue
            c_cons[cn, ti, tj] = model.addConstr(x.sum("*", cn, ti, tj) - x.sum(cn, "*", ti, tj) == 0, name="c_{0}_{1}_{2}".format(cn, ti, tj))

    # Flow conservation at edge switches
    tf = gp.tupledict()
    for ti in range(numTrunk):
        tn = "t"+str(ti)
        for tj in range(numTrunk):
            tm = "t"+str(tj)
            if ti == tj:
                assert tmat[ti, tj] == 0
                continue
            tf[tn, tm] = np.float64(tmat[ti, tj])

    numEdgeWire = gp.tupledict()
    for en in edgeNodes:
        for tn in trunks:
            numEdgeWire[en, tn] = wiring[int(en[1:])][int(tn[1:])]
    numResWire = gp.tupledict()
    for tn in trunks:
        numResWire[tn] = 0
        for en in edgeResNodes:
            numResWire[tn] += numEdgeWire[en, tn]
        if numResWire[tn] == 0:
            return 0


    e_cons = gp.tupledict()
    for en in edgeResNodes:
        for (ti, tj) in trunkPairs:
            if ti == tj:
                continue
            inflow = x.sum("*", en, ti, tj) + beta * numEdgeWire[en, ti] * tf[ti, tj] / numResWire[ti]
            outflow = x.sum(en, "*", ti, tj) + beta * numEdgeWire[en, tj] * tf[ti, tj] / numResWire[tj]
            e_cons[en, ti, tj] = model.addConstr(inflow == outflow, name="e_{0}_{1}_{2}".format(en, ti, tj))

    # Trunk capacity constraints
    t_cons = gp.tupledict()
    for en in edgeResNodes:
        for (ti, tj) in trunkPairs:
            if ti == tj:
                continue
            t_cons[en, ti, tj] = model.addConstr(beta * tf[ti, tj] / numResWire[ti] <= 1, name="t_{0}_{1}_{2}".format(en, ti , tj))

    # Link capacity constraints
    l_cons = gp.tupledict()
    failedLinks = router.getAllFailedLinks()
    for l in unLinks:
        (en, cn) = l        
        if l in failedLinks:
            cap = 0
        else:
            cap = 1
        l_cons[en, cn] = model.addConstr(x.sum(en, cn, "*", "*") <= cap, name="l_{0}_{1}".format(en, cn))
        l_cons[cn, en] = model.addConstr(x.sum(cn, en, "*", "*") <= cap, name="l_{0}_{1}".format(cn, en))

    # call solver
    model.optimize()
    if model.status == gp.GRB.Status.OPTIMAL:
        return beta.X
    elif model.status == gp.GRB.Status.INFEASIBLE:
        assert False, "Infeasible solution"
        return -1
    else:
        assert False, "Solver Status " + str(model.status)
        return -1


def verificationRouter(scenario):
    router = scenario["router"]
    wiring = scenario["wiring"]
    M = scenario["M"]
    
    etms = load_ETM(M)
    beta = 1
    for etm in etms:
        beta = min(beta, calculateEffectiveCapacityRouter(router, wiring, etm))
    return beta

    
def simpleTest():
    radix = 4
    M = (2, 2, 4)
    algo = "Optimal"
    numL2Failure = 0
    numL1Failure = 0
    numLinkFailure = 1
    
    router = Router(radix)
    wiring = load_Wiring(M, algo)
    etms = load_ETM(M)

    l2fp = set(["c"+str(ci) for ci in range(numL2Failure)])
    print("L2 FPS:", l2fp)
    l1fps = getL1FP(wiring, numL1Failure)
    print("L1 FPS:", l1fps)    
    linkfps = load_RepLinkFP(M, numLinkFailure, algo)
    print("Link FPS:", linkfps)

    for l1fp in l1fps:
        for linkfp in linkfps:        
            router.setL2Failure(l2fp)
            router.setL1Failure(l1fp)
            router.setLinkFailure(linkfp)

            scenario = { "router": router,
                         "wiring": wiring,
                         "M": M }
            beta = verificationRouter(scenario)
            print("Beta =", beta, ": L2-FP =", l2fp, ", L1-FP =", l1fp, ", LK-FP =", linkfp)


if __name__ == "__main__":
    simpleTest()
