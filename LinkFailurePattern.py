from __future__ import print_function
import itertools
import copy
import time
import numpy as np
import pickle

from Wiring import load_Wiring, getBaselineWiring

LinkFPSaveDir = "LKFP-save"


def defaultRepLinkFPFilename(M, numFailure, algo):
    fname = "RLKFP"
    for m in M:
        fname += "_" + str(m)
    return fname + "-" + str(numFailure) + "-" + algo


def load_RepLinkFP(M, numFailure, algo):
    fname = defaultRepLinkFPFilename(M, numFailure, algo)
    fullpath = LinkFPSaveDir + "/" + fname
    fpsFile = open(fullpath, "r")
    fps = pickle.load(fpsFile)
    fpsFile.close()
    return fps


def save_RepLinkFP(rfps, M, numFailure, algo):
    fname = defaultRepLinkFPFilename(M, numFailure, algo)
    fullpath = LinkFPSaveDir + "/" + fname
    fpsFile = open(fullpath, "w")
    pickle.dump(rfps, fpsFile)
    fpsFile.close()
    print("Save result to", fullpath)


class FailurePattern():
    def __init__(self, radix, wiring):
        self.wiring = tuple([tuple(ew) for ew in wiring])
        self.numEdge = radix
        self.numIntf = radix/2
        self.numCore = radix/2

        self.L1gs = None
        self.L1gsKey = None
        self.L2gs = None
        self.numFailure = None        
        self.L1FailureGroups = None
        self.L1FailureSubPatternTemplate = None
        self.L1FailureSubPatterns = None
        self.L1FailurePatterms = None
        self.L2FailurePatterns = None
        self.FailurePatterns = None

        self.__decomposeComponent()

        
    def __decomposeComponent(self):
        # L1 groups
        wiring = self.wiring
        L1gs = dict()
        L1gsKey = list()
        for e in xrange(len(wiring)):
            we = wiring[e]
            if not L1gs.has_key(we):
                L1gs[we] = list()
                L1gsKey.append(we)
            edgeID = e
            L1gs[we].append(edgeID)
        self.L1gs = L1gs
        self.L1gsKey = L1gsKey

        # L2 groups
        self.L2gs = tuple([c for c in xrange(self.numCore)])


    def __recurCreateL1FailureGroups(self, remainFailure, fg):
        if len(fg) == len(self.L1gs):
            if remainFailure == 0:
                self.L1FailureGroups.append(tuple(fg))
            return

        curGroupKey = self.L1gsKey[len(fg)]
        for i in xrange(min(self.numIntf * len(self.L1gs[curGroupKey]), remainFailure) + 1):
            tfg = fg[:]
            tfg.append(i)
            self.__recurCreateL1FailureGroups(remainFailure-i, tfg)


    def __recurGenL1FailureSubPatternTemplate(self, totEdge, totFailure, remainFailure, fe):
        if len(fe) == totEdge:
            if remainFailure == 0:
                tempat = (totEdge, totFailure)
                self.L1FailureSubPatternTemplate[tempat].append(tuple(fe))
            return

        prvFailureCnt = float("inf")
        if len(fe) > 0:
            prvFailureCnt = fe[-1]
        for i in xrange(min([prvFailureCnt, self.numIntf, remainFailure]), -1, -1):
            tfe = fe[:]
            tfe.append(i)
            self.__recurGenL1FailureSubPatternTemplate(totEdge, totFailure, remainFailure-i, tfe)

            
    def __createL1FailureGroups(self):
        self.L1FailureGroups = list()
        self.__recurCreateL1FailureGroups(self.numFailure, list())
        self.L1FailureGroups = tuple(self.L1FailureGroups)


    def __genL1FailureSubPatternTemplate(self):
        self.L1FailureSubPatternTemplate = dict()
        for fg in self.L1FailureGroups:
            for i in xrange(len(fg)):
                curGroupKey = self.L1gsKey[i]
                totEdge = len(self.L1gs[curGroupKey])
                totFailure = fg[i]
                tempat = (totEdge, totFailure)
                if not self.L1FailureSubPatternTemplate.has_key(tempat):
                    self.L1FailureSubPatternTemplate[tempat] = list()
                    self.__recurGenL1FailureSubPatternTemplate(totEdge, totFailure, totFailure, list())
                    self.L1FailureSubPatternTemplate[tempat] = tuple(self.L1FailureSubPatternTemplate[tempat])


    def __createL1FailurePatterns(self):
        self.L1FailurePatterns = dict()
        for fg in self.L1FailureGroups:
            # print("Failure group", fg)
            numSubGroup = len(fg)

            # Get a list of L1 failed switches
            l1s = dict()
            l1subpats = dict()

            for si in xrange(numSubGroup):
                l1s[si] = dict()
                nodes = self.L1gs[self.L1gsKey[si]]
                numSubL1 = fg[si]                
                subpat = (len(nodes), numSubL1)                
                l1subpats[si] = self.L1FailureSubPatternTemplate[subpat]
                for l1sp in l1subpats[si]:
                    l1s[si][l1sp] = nodes[:len(l1sp)]

            self.L1FailureSubPatterns = l1subpats                    
            # print("\t L1 subpatterns", l1subpats)
            # print("\t L1 nodes", l1s)

            l1pats = { tuple(): list() }
            for si in xrange(numSubGroup):
                tl1pats = copy.deepcopy(l1pats)
                l1pats = dict()
                for pat, patnodes in tl1pats.iteritems():
                    for l1sp in l1subpats[si]:
                        p = list(pat)
                        p.append(l1sp)
                        nodes = l1s[si][l1sp]
                        t = list()
                        t.extend(patnodes)
                        for i in xrange(len(l1sp)):
                            for rep in xrange(l1sp[i]):
                                t.append([nodes[i], None])
                        l1pats[tuple(p)] = tuple(t)

            self.L1FailurePatterns[fg] = l1pats
            # print("\t L1 patterns count", len(l1pats))
            # print("\t L1 patterns", l1pats)


    def __genL2FailureSubPatternTemplate(self):
        self.L2FailureSubPatternTemplate = dict()
        for tempat, subpats in self.L1FailureSubPatternTemplate.iteritems():
            # print("Template", tempat, ":", subpats)
            (numSubNode, numFailure) = tempat
            self.L2FailureSubPatternTemplate[tempat] = dict()
            for subpat in subpats:
                # print("\t Subpattern", subpat)
                numActSubNode = len(subpat)
                bin = ((0,)*numActSubNode,)*self.numCore
                bins = set([bin])
                for pi in xrange(numActSubNode):
                    # print("\t\t Active index", pi)
                    tbins = bins.copy()
                    bins = set()
                    numLink = subpat[pi]
                    perit = itertools.combinations(range(self.numCore), numLink)
                    for per in perit:
                        # print("\t\t\t Combination", per)
                        for tbin in tbins:
                            ltbin = [list(b) for b in tbin]
                            for idx in per:
                                ltbin[idx][pi] += 1
                            ltbin.sort(reverse=True)
                            bins.add(tuple([tuple(b) for b in ltbin]))
                #     print("\t\t Active subnode", pi, bins)
                # print("\t Bins", bins)
                self.L2FailureSubPatternTemplate[tempat][subpat] = bins


    def __createL2FailurePatterns(self):
        self.L2FailurePatterns = dict()
        for fg in self.L1FailureGroups:
            # print("Failure group", fg)
            self.L2FailurePatterns[fg] = dict()
            l1dists = self.L1FailurePatterns[fg].keys()
            l2subpats = dict()
            for l1dist in l1dists:
                # print("\t L1 failure distribution", l1dist)
                pats = set([tuple([() for x in xrange(self.numCore)])])                
                for sgidx in xrange(len(l1dist)):
                    # print("\t\t L1 index", sgidx)
                    tpats = pats.copy()
                    pats = set()
                    numNode = len(self.L1gs[self.L1gsKey[sgidx]])
                    numFailure = fg[sgidx]
                    l1subpat = l1dist[sgidx]
                    l2subpats[sgidx] = list(self.L2FailureSubPatternTemplate[(numNode, numFailure)][l1subpat])
                    # print("\t\t L1 subpat", l1subpat)
                    # print("\t\t L2 subpats", l2subpats[sgidx])
                    for l2subpat in l2subpats[sgidx]:
                        # print("\t\t\t L2 subpat", l2subpat)

                        # Create type and index to reduce duplicate
                        l2type = dict()
                        for t in l2subpat:
                            tt = list(t)
                            tt.sort(reverse=True)
                            tt = tuple(tt)
                            if not l2type.has_key(tt):
                                l2type[tt] = len(l2type)
                        used = set()
                        
                        iter = itertools.permutations(range(self.numCore))
                        for it in iter:

                            # Detect duplicate
                            v = list()
                            for i in xrange(len(it)):
                                y = list(l2subpat[it[i]])
                                y.sort(reverse=True)
                                y = tuple(y)
                                v.append(l2type[y])
                            # print("\t\t\t\t Permutation index", it, "-->", v)
                            if tuple(v) in used:
                                # print("\t\t\t\t\t Skip")
                                continue
                            else:
                                used.add(tuple(v))
                            
                            # print("\t\t\t\t Permutation index", it)
                            for tpat in tpats:
                                # print("\t\t\t\t\t init pattern", tpat)
                                npat = [list(t) for t in tpat]
                                for x in xrange(self.numCore):
                                    npat[x].extend(l2subpat[it[x]])
                                # print("\t\t\t\t\t npat", npat)
                                npat.sort(reverse=True)
                                pats.add(tuple([tuple(p) for p in npat]))
                self.L2FailurePatterns[fg][l1dist] = pats
                # print("\t\t L2 patterns", pats)


    def __constructRepresentativePatterns(self):
        self.FailurePatterns = set()
        for fg in self.L1FailureGroups:
            # print("Failure group", fg)
            l1dists = self.L1FailurePatterns[fg].keys()            
            for l1dist in l1dists:
                # print("\t L1 failure distribution", l1dist)
                # print("\t\t L1 failure patterns", self.L1FailurePatterns[fg][l1dist])
                # print("\t\t L2 failure patterns", self.L2FailurePatterns[fg][l1dist])
                for l2pat in self.L2FailurePatterns[fg][l1dist]:
                    fp = list()
                    for c in xrange(self.numCore):
                        for e in xrange(self.numEdge):
                            if l2pat[c][e] == 1:
                                fp.append(("e" + str(e), "c" + str(c)))
                    # print("\t\t\t", tuple(fp))
                    self.FailurePatterns.add(tuple(fp))
                

    def setNumFailure(self, numFailure):
        self.numFailure = numFailure
        
        self.__createL1FailureGroups()
        self.__genL1FailureSubPatternTemplate()
        self.__createL1FailurePatterns()

        self.__genL2FailureSubPatternTemplate()        
        self.__createL2FailurePatterns()
        self.__constructRepresentativePatterns()
        return self.FailurePatterns

    
    def getRepresentative(self, fp):
        numGroup = len(self.L1gsKey)        
        ifp = StrFPtoIntFP(fp)
        # print("Input faiure pattern", fp)
        # print("\t Integer failure pattern", ifp)

        re = list()
        for gid in xrange(numGroup):
            k = self.L1gsKey[gid]
            subgrouplen = len(self.L1gs[k])
            re.append([0]*subgrouplen)
        rc = [[0]*self.numEdge for x in xrange(self.numCore)]
        
        for l in ifp:
            eid = l[0]
            cid = l[1]
            gid = self.L1gsKey.index(self.wiring[eid])
            egid = self.L1gs[self.L1gsKey[gid]].index(eid)
            re[gid][egid] += 1
            rc[cid][eid] += 1
        # print("\t Edge groups", re)
        # print("\t Core group", rc)

        # Sort layer-1 subgroup
        offset = 0
        newcg = [ [ [] for gid in xrange(len(re)) ] for x in xrange(self.numCore)]
        neweg = [None]*len(re)
        newfp = list()
        for gid in xrange(numGroup):
            sg = re[gid]
            lensg = len(sg)            
            # print("subgroup:", sg)
            sortidx = list(np.argsort(sg))
            # print("\t sort index", sortidx)

            newsg = [None] * lensg
            for idx in xrange(lensg):
                ridx = lensg - idx - 1
                pidx = sortidx[ridx]
                newsg[idx] = sg[pidx]
                for l in ifp:
                    if l[0] != offset + pidx:
                        continue
                    cidx = l[1]
                    newcg[cidx][gid].append(sg[pidx])
                    newfp.append((offset + idx, cidx))
            neweg[gid] = newsg
            offset += lensg

        # Debug
        # print("neweg:", neweg)
        # print("newcg:", newcg)
        # print("newfp:", newfp)       

        # Sort layer-2 subgroup
        neweg2 = copy.deepcopy(neweg)        
        newcg2 = copy.deepcopy(newcg)
        newcg2.sort(reverse=True)
        newfp2 = list()
        # create mapping
        map2 = dict() # key = coreid, new coreid
        dupcheck = dict() # key = newcg[x], value = last index
        for cid in xrange(self.numCore):
            code = newcg[cid]
            tcode = tuple(tuple(x) for x in newcg[cid])            
            if dupcheck.has_key(tcode):
                lastindex = dupcheck[tcode]
            else:
                lastindex = -1
            newindex = newcg2.index(code, lastindex+1)
            map2[cid] = newindex
            dupcheck[tcode] = newindex
        # map newfp to newfp2
        newfp2 = [(l[0], map2[l[1]]) for l in newfp]

        # Debug
        # print("neweg2:", neweg2)
        # print("newcg2:", newcg2)
        # print("newfp2:", newfp2)

        # Sort layer-1 subgroup's group
        offset = 0
        newcg3 = [[0]*self.numEdge for x in xrange(self.numCore)]
        neweg3 = list()
        newfp3 = list()
        for gid in xrange(numGroup):
            sg = neweg[gid]
            # print("gid:", gid, sg)
            begidx = 0
            # Iterate by subsubgroup
            newsg = list()
            while begidx != len(sg):
                endidx = begidx + 1
                while endidx < len(sg) and sg[begidx] == sg[endidx]:
                    endidx += 1
                # print("\t subsubgroup index range", (begidx, endidx))
                lenssg = endidx - begidx
                ssgcode = list()

                # Create code of level-1 switch 
                for idx in xrange(begidx, endidx):
                    code = [0]*self.numCore
                    for l in newfp2:
                        if l[0] != offset + idx:
                            continue
                        cidx = l[1]
                        code[cidx] += 1
                    ssgcode.append(code)
                # print("\t subsubgroup code", ssgcode)

                tssgcode = copy.deepcopy(ssgcode)
                tssgcode.sort(reverse=True)
                # print("\t sorted code", tssgcode)

                # Update fp
                dupcheck = dict()
                # create map3
                map3 = dict() # key = edgeid, value = new edgeid
                dupcheck = dict() # key = ssgcode[x], value = last index
                for ssgid in xrange(len(ssgcode)):
                    code = ssgcode[ssgid]
                    tcode = tuple(code)
                    if dupcheck.has_key(tcode):
                        lastindex = dupcheck[tcode]
                    else:
                        lastindex = -1
                    newindex = tssgcode.index(code, lastindex+1)
                    map3[ssgid] = newindex
                    dupcheck[tcode] = newindex
                # map newfp to newfp2
                for l in newfp2:
                    if l[0] < begidx + offset or l[0] >= endidx + offset:
                        continue
                    newfp3.append( (offset + begidx + map3[l[0] - offset - begidx], l[1]) )

                # print("\t newfp", newfp3)

                # Update l2 code
                
                newsg.extend(tssgcode)
                begidx = endidx                

            neweg3.append(newsg)
            offset += len(sg)
        newfp3.sort()
        return tuple([("e" + str(l[0]), "c" + str(l[1])) for l in newfp3])


def StrFPtoIntFP(sfp):
    ifp = list()
    for l in sfp:
        if l[0][0] == "e" and l[1][0] == "c":
            pass
        elif l[0][0] == "c" and l[1][0] == "e":
            l = (l[1], l[0])
        else:
            assert False, "Invalid failure link {0}".format(str(l))
        ifp.append((int(l[0][1:]), int(l[1][1:])))
    return tuple(ifp)


def generateRepLinkFailurePatterns(radix, M, algo, minFailedLink, maxFailedLink):
    wiring = load_Wiring(M, algo)
    FP = FailurePattern(radix, wiring)

    assert minFailedLink >= 0
    upperboundFailure = radix**2 / 2
    assert maxFailedLink <= upperboundFailure

    for numFailure in range(minFailedLink, maxFailedLink + 1):
        t = time.time()
        fps = FP.setNumFailure(numFailure)
        repfps = set()
        for fp in fps:
            repfps.add(FP.getRepresentative(fp))

        save_RepLinkFP(repfps, M, numFailure, algo)


if __name__ == "__main__":
    radix = 4
    M = (2, 2, 4)
    algo = "Optimal"
    minFailedLink = 1
    maxFailedLink = 1
    generateRepLinkFailurePatterns(radix, M, algo, minFailedLink, maxFailedLink)

