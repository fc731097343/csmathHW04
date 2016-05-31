import math
import numpy as np
import random

#find one extremum for function y = a *sin(b*x)
a = 2
b = np.matrix([1,2])

x = np.matrix([[3],[2.5]])

gk_norm = 100
k = 0
epson = 1e-13
mu = 0.1
fk = a * math.sin(b * x)
maxk = 1000
qk = 1
while(gk_norm > epson and k < maxk):
    gk = a * b * math.cos(b * x)
    Gk = -1 * a * b.T * b * math.sin(b * x)
    gk_norm = gk * gk.T
    [S, V, D] = np.linalg.svd(Gk+mu*np.eye(2), full_matrices=True)
    while(0 in V):
        mu = 4*mu
        [S, V, D] = np.linalg.svd(Gk+mu*np.eye(2), full_matrices=True)
    sk = np.array((np.linalg.inv(Gk+mu*np.eye(2)) * (-gk.T)).T)[0]

    fk_old = fk
    qk_old = qk

    fk = a * math.sin(b * (np.matrix(sk).T + x))

    qk = (fk + gk * np.matrix(sk).T + 0.5 * np.matrix(sk) * Gk * np.matrix(sk).T)[0,0]

    rk = (fk-fk_old) / (qk - qk_old)
    if rk < 0.25:
        mu = 4 * mu
    elif rk > 0.75:
        mu = 0.5 * mu

    if rk > 0:
        x = x + np.matrix(sk).T

    k += 1
print "iterations : %d, gk norm = %e, x = [%f, %f]" % (k,gk_norm,x[0,0],x[1,0])

print "extremum = 2 or -2, estimate = a*sin(b*x) = " , (a * math.sin(b * x))
