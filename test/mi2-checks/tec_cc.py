#!/usr/bin/env python3
"""Read a hydra/ICE Tecplot BLOCK file with NODAL x,y,z + CELLCENTERED variables -> dict name -> (I-1,J-1,K-1) array.
Usage as a module: read_cc(path) ; as a script: tec_cc.py <file> <var> [i-index (0-based)] -> mean over j,k of var at that i."""
import re, sys
import numpy as np

def read_cc(path):
    txt = open(path).read()
    names = re.findall(r'"([^"]+)"', re.search(r'VARIABLES\s*=\s*(.*)', txt).group(1))
    z = re.search(r'ZONE.*?I\s*=\s*(\d+).*?J\s*=\s*(\d+).*?K\s*=\s*(\d+)', txt, re.S)
    I, J, K = map(int, z.groups())
    body = txt[z.end():]; body = body[body.index('\n') + 1:]
    vals = np.array([float(v) for v in re.findall(r'[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?', body.replace('D', 'E'))])
    nn, nc = I * J * K, (I - 1) * (J - 1) * (K - 1)
    out, p = {}, 0
    for n in names:
        cnt = nn if n in ('x', 'y', 'z') else nc
        a = vals[p:p + cnt]; p += cnt
        out[n] = a.reshape((K, J, I), order='C') if cnt == nn else a.reshape((K - 1, J - 1, I - 1), order='C')
    return out

if __name__ == '__main__':
    d = read_cc(sys.argv[1]); v = sys.argv[2]; i = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    a = d[v]
    print(f'{v}[i={i}] mean={a[:, :, i].mean():.6e} min={a[:, :, i].min():.6e} max={a[:, :, i].max():.6e}')
