from __future__ import print_function
import time
import numpy as np
import pickle
import cdd

ETMSaveDir = "ETM-save"

def ETMFilename(M):
    fname = 'ETM'
    for m in M:
        fname += '_' + str(m)
    return fname


def load_ETM(M):
    fname = ETMFilename(M)
    fullpath = ETMSaveDir + "/" + fname
    fp = open(fullpath, 'r')
    etms = pickle.load(fp)
    fp.close()
    return etms


def save_ETM(wiring, M):
    fname = ETMFilename(M)
    fullpath = ETMSaveDir + "/" + fname
    fp = open(fullpath, 'w')
    pickle.dump(wiring, fp)
    fp.close()
    print("Save result to", fullpath)
    

def generate_ETM(M):
    K = len(M)
    nEle = K**2

    dt = time.time()
    rows = list()

    # row inequality
    for r in range(K):
        hfmt = [0]*(nEle + 1)
        # 1) b + ax >= 0 ==> [b, a1, a2, ...]
        # 2) ax <= b     ==> [b, -a1, -a2, ...]
        hfmt[0] = M[r]
        for c in range(K):
            hfmt[1+r*K+c] = -1
        rows.append(hfmt)

    # column inequality
    for c in range(K):
        hfmt = [0]*(nEle + 1)
        hfmt[0] = M[c]
        for r in range(K):
            hfmt[1+r*K+c] = -1
        rows.append(hfmt)

    # non-negativity
    for e in range(nEle):
        hfmt = [0]*(nEle + 1)
        hfmt[1+e] = 1
        rows.append(hfmt)

    # no reverse traffic
    for k in range(K):
        hfmt = [0]*(nEle + 1)
        hfmt[1+K*k+k] = -1
        rows.append(hfmt)

    mat = cdd.Matrix(rows, number_type='fraction')
    mat.rep_type = cdd.RepType.INEQUALITY
    poly = cdd.Polyhedron(mat)
    vs = poly.get_generators()
    dt = time.time() - dt
    print('Execution duration', dt, 'sec.')

    etms = list()
    for v in vs:
        etm = np.array(v[1:], dtype='float64')
        etm.resize((K, K))
        etms.append(etm)
    save_ETM(etms, M)

    
if __name__ == '__main__':
    M = (2, 2, 4)
    generate_ETM(M)
