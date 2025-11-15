#!/usr/bin/env python3
import sys
import os

def load_pot(path):
    xs, vs = [], []
    with open(path, 'r') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#') or s.startswith(';'):
                continue
            parts = s.split()
            if len(parts) >= 2:
                try:
                    xs.append(float(parts[0]))
                    vs.append(float(parts[1]))
                except:
                    pass
    return xs, vs

def deriv(xs, vs):
    n = len(xs)
    fs = [0.0]*n
    for i in range(n):
        if i == 0:
            dx = xs[i+1] - xs[i]
            dv = vs[i+1] - vs[i]
        elif i == n-1:
            dx = xs[i] - xs[i-1]
            dv = vs[i] - vs[i-1]
        else:
            dx = xs[i+1] - xs[i-1]
            dv = vs[i+1] - vs[i-1]
        dVdr = dv/dx
        fs[i] = -dVdr
    return fs

def write_xvg(path, xs, vs, fs):
    with open(path, 'w') as f:
        f.write("# GROMACS table: r  V(r)  F(r)=-dV/dr\n")
        for r, v, force in zip(xs, vs, fs):
            f.write(f"{r:.6f} {v:.8f} {force:.8f}\n")

def main():
    if len(sys.argv) < 3:
        print("usage: pot_to_table.py <input.pot> <output.xvg>")
        sys.exit(1)
    inp = sys.argv[1]
    outp = sys.argv[2]
    xs, vs = load_pot(inp)
    if not xs:
        print("no data in pot")
        sys.exit(2)
    fs = deriv(xs, vs)
    write_xvg(outp, xs, vs, fs)

if __name__ == '__main__':
    main()