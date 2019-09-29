from __future__ import print_function
import numpy as np

class FailureGraph():
    def __init__(self, radix, wiring):
        self.radix = radix
        self.numEdge = radix
        self.numCore = radix / 2
        self.numIntf = radix / 2

        lwiring = list(wiring)
        lwiring.sort()
        self.wiring = tuple(tuple(ew) for ew in lwiring)

        self.lkfp = set()
        self.replkfp = set()

        self.__createL1Groups()

        
    def __createL1Groups(self):
        self.L1GrpMap = dict()
        for ew in self.wiring:
            if self.L1GrpMap.has_key(ew):
                continue
            grpId = len(self.L1GrpMap)
            self.L1GrpMap[ew] = grpId

        self.L1GrpInfo = dict()
        for ewidx in range(self.numEdge):
            ew = self.wiring[ewidx]
            grpId = self.L1GrpMap[ew]
            if not self.L1GrpInfo.has_key(grpId):
                self.L1GrpInfo[grpId] = { "L1Index": [ewidx], "count": 1 }
            else:
                self.L1GrpInfo[grpId]["L1Index"].append(ewidx)
                self.L1GrpInfo[grpId]["count"] += 1

                
    def setFailedLinks(self, str_links):
        for l in str_links:
            assert set([l[0][0], l[1][0]]) == set(["e", "c"])
            if l[0][0] == "c":
                l = (l[1], l[0])
            self.lkfp.add( (int(l[0][1:]), int(l[1][1:])) )
        self.replkfp = None


    def getRepFailedLinks(self):
        if self.replkfp != None:
            return self.replkfp.copy()
        self.__L1sort()
        print("L1 sort", self.lkfp_l1sort)
        self.__L2sort()
        self.__L1resort()
        # return self.replkfp.copy()


    def __L1sort(self):
        # sort L1 per ngroup
        newlkfp = set()        
        for gid in self.L1GrpMap.values():
            l1Indexs = self.L1GrpInfo[gid]["L1Index"][:]
            l1Cards = dict( zip(l1Indexs, [0]*self.L1GrpInfo[gid]["count"]) )
            for eid in l1Indexs:
                for l in self.lkfp:
                    if l[0] == eid:
                        l1Cards[eid] += 1

            l1CardList = list()
            for eid in l1Indexs:
                l1CardList.append(l1Cards[eid])
            sortidx = np.argsort(l1CardList)
            sortmapping = dict()            
            for i in range(len(sortidx)-1, -1, -1):
                prvidx = sortidx[i]
                prvNodeId = l1Indexs[prvidx]
                newidx = len(sortidx)-1 - i
                newNodeId = l1Indexs[newidx]
                sortmapping[prvNodeId] = newNodeId
            for l in self.lkfp:
                if l[0] in sortmapping.keys():
                    newl = (int(sortmapping[l[0]]), int(l[1]))
                    newlkfp.add(newl)
        self.lkfp_l1sort = newlkfp
    
    def __L2sort(self):
        l2tuple = tuple([tuple([[0]*self.L1GrpInfo[gid]["count"] for gid in self.L1GrpMap.values()]) for cid in range(self.numCore)])
        for cid in range(self.numCore):
            for gid in self.L1GrpMap.values():
                for e in range(self.L1GrpInfo[gid]["count"]):
                    eid = self.L1GrpInfo[gid]["L1Index"][e]
                    fplink = (eid, cid)
                    if fplink not in self.lkfp_l1sort:
                        continue
                    for l in self.lkfp_l1sort:
                        if l[0] == eid:
                            l2tuple[cid][gid][e] += 1

        print(l2tuple)
        
        
    
    def __L1resort(self):
        pass



def test_FailureGraph():
    radix = 4
    wiring = ((1, 0, 1), (0, 1, 1), (1, 0, 1), (0, 1, 1))
    fg = FailureGraph(radix, wiring)

    assert fg.L1GrpMap == { (0, 1, 1): 0,
                            (1, 0, 1): 1 }
    assert fg.L1GrpInfo == { 0: {'L1Index': [0, 1], 'count': 2},
                             1: {'L1Index': [2, 3], 'count': 2} }

    lkfp = [("e1", "c1"), ("e3", "c0")]
    fg.setFailedLinks(lkfp)
    assert fg.lkfp == set([(1, 1), (3, 0)])

    fg.getRepFailedLinks()
    assert fg.lkfp_l1sort == set([(0, 1), (2, 0)])

    


    return fg

if __name__ == "__main__":
    fg = test_FailureGraph()
