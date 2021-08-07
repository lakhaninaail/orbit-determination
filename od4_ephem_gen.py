from odlib import *
import numpy as np
import math

a = semimajorAxis("LakhaniInput.txt", "00:00:00.0000")
e = eccentricity("LakhaniInput.txt", "00:00:00.0000")
i = inclination("LakhaniInput.txt", "00:00:00.0000") * math.pi / 180
loan = LoAN("LakhaniInput.txt", "00:00:00.0000") * math.pi / 180
aop = AoP("LakhaniInput.txt", "00:00:00.0000") * math.pi / 180
#M = meanAnomaly("LakhaniInput.txt", "00:00:00.0000")
T = LPP("LakhaniInput.txt", "00:00:00.0000")
R = getSunVector("LakhaniInput1.txt", "00:00:00.0000")
#print("T",T)
#print("M",M)

t = 2458333.5 #JD of observation
k = .01720209894
n  = k * (1 / a ** 3) ** 0.5
M = n * (t - T)
#print("M",M)


def newtonRaphson(M, e, tol):
    E_predict = M
    f = M - (E_predict - e * math.sin(E_predict))
    iterations = 0
    while(abs(f) >= tol):
        iterations += 1
        #f = M - (E_predict - e * math.sin(E_predict))
        E_predict = E_predict - (M - (E_predict - e * math.sin(E_predict))) / (e * math.cos(E_predict) - 1)
        f = M - (E_predict - e * math.sin(E_predict))
    return E_predict

E = newtonRaphson(M, e, 1E-12)
#print("E",E)

def ephem_gen(a, R):
    obliquity = 23.4 * math.pi / 180
    x = a * math.cos(E) - a * E
    y = a * (1 - e ** 2) ** 0.5 * math.sin(E)
    z = 0
    r = [[x], [y], [z]]
    rnew = [x, y, z]
    rnew = np.array(rnew)
    rot_mat_aop = [[math.cos(aop), -math.sin(aop), 0],
                   [math.sin(aop), math.cos(aop), 0], 
                   [0, 0, 1]]
    rot_mat_i = [[1, 0, 0],
                 [0, math.cos(i), -math.sin(i)], 
                 [0, math.sin(i), math.cos(i)]]
    rot_mat_loan = [[math.cos(loan), -math.sin(loan), 0],
                   [math.sin(loan), math.cos(loan), 0], 
                   [0, 0, 1]]
    r = np.array(r)
    rot_mat_aop = np.array(rot_mat_aop)
    rot_mat_i = np.array(rot_mat_i)
    rot_mat_loan = np.array(rot_mat_loan)

    #print(type(r))
    #print(type(rot_mat_aop))
    #print(np.shape(r))
    #print(np.shape(rot_mat_aop))

    a = np.dot(rot_mat_aop, r)
    b = np.dot(rot_mat_i, a)
    r1 = np.dot(rot_mat_loan, b)
    #print(np.shape(b))
    #print(b)
    #print(np.shape(r))
    #print(r)
    #print(r1)
    #print(r2)
    rot_mat_obl = [[1, 0, 0],
                 [0, math.cos(obliquity), -math.sin(obliquity)], 
                 [0, math.sin(obliquity), math.cos(obliquity)]]
    
    rot_mat_obl = np.array(rot_mat_obl)
    r2 = np.dot(rot_mat_obl, r1)
    rho = []
    for k in range(3):
        rho.append(R[k] + r2[k])
    rho = np.array(rho)
    rho_mag = (rho[0] ** 2 + rho[1] ** 2 + rho[2] ** 2) ** 0.5
    rho_hat = 1 / rho_mag * rho
    print(rho_hat)
    print(rho_hat[2])
    dec = math.asin(rho_hat[2])
    ra = math.acos(rho_hat[0] / math.cos(dec))
    #print(ra)

    return ra * 180 / math.pi, dec * 180 / math.pi

print(ephem_gen(a, R))