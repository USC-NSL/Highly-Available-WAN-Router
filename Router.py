class Router:
    def __init__(self, radix):
        self.radix = radix
        self.numEdge = radix
        self.numCore = radix/2
        self.failedLinks = set()
        self.failedL1s = set()
        self.failedL2s = set()

    def setLinkFailure(self, failedLinks):
        links = set()
        for l in failedLinks:
            if set([l[0][0], l[1][0]]) != set(['e', 'c']):
                print("Invalid failed link", l)
                return
            if l[0][0] == 'e':
                tl = (l[0], l[1])
            else:
                tl = (l[1], l[0])
            if int(tl[0][1:]) < 0 or int(tl[0][1:]) > self.numEdge:
                print("Invalid failed link", l)
                return
            if int(tl[1][1:]) < 0 or int(tl[1][1:]) > self.numCore:
                print("Invalid failed link", l)
                return
            links.add(tl)
        self.failedLinks = links

    def setL1Failure(self, failedL1s):
        l1s = set()
        for e in failedL1s:
            if e[0] != 'e':
                print("Invalid L1 node", e)
                return
            if int(e[1:]) < 0 or int(e[1:]) >= self.numEdge:
                print("Invalid L1 node", e)
                return
            l1s.add(e)
        self.failedL1s = l1s

    def setL2Failure(self, failedL2s):
        l2s = set()
        for c in failedL2s:
            if c[0] != 'c':
                print("Invalid L2 node", c)
                return
            if int(c[1:]) < 0 or int(c[1:]) >= self.numCore:
                print("Invalid L2 node", c)
                return
            l2s.add(c)
        self.failedL2s = l2s

    def getAllFailedLinks(self):
        links = self.failedLinks.copy()
        for c in self.failedL2s:
            for ei in range(self.numEdge):
                l = ("e" + str(ei), c)
                links.add(l)
        for e in self.failedL1s:
            for ci in range(self.numCore):
                l = (e, "c" + str(ci))
                links.add(l)
        return links
