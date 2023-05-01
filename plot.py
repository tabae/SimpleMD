#coding: utf-8

import matplotlib.pyplot as plt

iters = []
ks = []
vs = []
kvs = []

with open('results.txt') as f:
    line = f.readline()
    while(line):
        iter, k, v, kv = map(float, line.split())
        iter = int(iter)
        iters.append(iter)
        ks.append(k)
        vs.append(v)
        line = f.readline()

plt.plot(iters, ks, label='k')
plt.plot(iters, vs, label='v')
plt.legend(loc='best')
plt.savefig('results.png')