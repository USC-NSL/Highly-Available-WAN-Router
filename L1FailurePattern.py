from __future__ import print_function
import itertools
import pickle

from ExtremeTraffic import load_ETM
from Wiring import load_Wiring


def deriveSubgroup(wiring):
    groups = dict()
    gkeys = dict()
    for eid in xrange(len(wiring)):
        we = wiring[eid]
        if not groups.has_key(we):
            gkeys[len(groups)] = we
            groups[we] = list()
        groups[we].append(eid)
    return (groups, gkeys)


def getL1FailurePatterns(groups, gkeys, numFailure):
    numGroup = len(groups)
    comb = [[]]
    for gid in xrange(numGroup):
        tgid = list()
        gk = gkeys[gid]
        for c in comb:
            for n in xrange(len(groups[gk])+1):
                tc = c[:]
                tc.append(n)
                if sum(tc) <= numFailure:
                    tgid.append(tc)
        comb = tgid

    fps = set()
    for c in comb:
        if sum(c) == numFailure:
            fps.add(tuple(c))
    return fps


def getL1FP(wiring, numFailure):
    (groups, gkeys) = deriveSubgroup(wiring)
    dists = getL1FailurePatterns(groups, gkeys, numFailure)
    fps = set()
    for dist in dists:
        fp = list()
        for gid in xrange(len(groups)):
            gk = gkeys[gid]
            for i in range(dist[gid]):
                fp.append("e"+str(groups[gk][i]))
        fps.add(tuple(fp))
    return fps


if __name__ == "__main__":
    radix = 4
    M = (2, 2, 4)
    algo = "Optimal"
    wiring = load_Wiring(M, algo)
    for numFailure in range(radix+1):
        fps = getL1FP(wiring, numFailure)
        print(numFailure, fps)
