from odlib import *
import numpy as np
import math

mu = 1
a = semimajorAxis("LakhaniInput.txt", "00:00:00.0000")
e = eccentricity("LakhaniInput.txt", "00:00:00.0000")
n = (mu / a ** 3) ** 0.5
tau3 = 0.050840808143482484
tau1 = -0.32618569435308475
r2 = [0.26640998194891174, -1.382856212643199, -0.505199925482389]
r2dot = [0.8439832722802604, -0.39937767878456487, 0.14200790188593015]
tol = 1E-12

def magnitude(vector):
    return (vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2) ** 0.5

def f(x, tau, r2, r2dot):
    return x - (1 - magnitude(r2) / a) * math.sin(x) + np.dot(r2, r2dot) / (n * a ** 2) * (1 - math.cos(x)) - n * tau

def df(x, tau, r2, r2dot):
    return 1 - (1 - magnitude(r2) / a) * math.cos(x) + np.dot(r2, r2dot) / (n * a ** 2) * math.sin(x)

def newtonRaphson(tau, e, tol, a, n):
    sign = np.dot(r2, r2dot) / (n * a ** 2) * math.cos(n * tau - np.dot(r2, r2dot) / (n * a ** 2)) + (1 - magnitude(r2) / a) * math.sin(n * tau - np.dot(r2, r2dot) / (n * a ** 2))
    plus_minus = math.copysign(1, sign)
    x = n * tau + plus_minus * 0.85 * e - np.dot(r2, r2dot) / (n * a ** 2)
    diff = f(x, tau, r2, r2dot) / df(x, tau, r2, r2dot)
    while(abs(diff) >= tol):
        diff = f(x, tau, r2, r2dot) / df(x, tau, r2, r2dot)
        x = x - diff
    return x

def fg(tau, r2, r2dot, flag):
    if flag == 0: #third order
        f = 1 - mu / (2 * magnitude(r2) ** 3) * tau ** 2 + mu * np.dot(r2, r2dot) / (2 * magnitude(r2) ** 5) * tau ** 3
        g = tau - mu * tau ** 3 / (6 * magnitude(r2) ** 3)
    
    if flag == 1: #fourth order
        u = 1 / magnitude(r2) ** 3
        z = np.dot(r2, r2dot) / magnitude(r2) ** 2
        q = np.dot(r2, r2dot) / magnitude(r2) ** 2 - u
        f = 1 - mu / (2 * magnitude(r2) ** 3) * tau ** 2 + mu * np.dot(r2, r2dot) / (2 * magnitude(r2) ** 5) * tau ** 3 + (3 * u * q - 15 * u * z ** 2 + u ** 2) / 24 * tau ** 4
        g = tau - mu * tau ** 3 / (6 * magnitude(r2) ** 3) + 6 * u * z * tau ** 4 / 24

    if flag == 2: # function
        delta_E = newtonRaphson(tau, e, tol, a, n)
        f = 1 - a * (1 - math.cos(delta_E)) / magnitude(r2)
        g = tau + (math.sin(delta_E) - delta_E) / n

    return f, g


#print(newtonRaphson(tau1, e, tol, a, n))
flag = int(input("Enter flag (0: third order, 1: fourth order, 2: function): "))
print("(f1, g1): ",format(fg(tau1, r2, r2dot, flag)))
print("(f3, g3): ",format(fg(tau3, r2, r2dot, flag)))

