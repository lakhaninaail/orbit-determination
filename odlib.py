import math
import numpy as np

def magnitude(vector):
    return (vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2) ** 0.5

def quadrantCheck(cos_angle, sin_angle):
    if sin_angle > 0 and cos_angle > 0:
        return math.asin(sin_angle) * 180 / math.pi
    if sin_angle > 0 and cos_angle < 0:
        return 180 - math.asin(sin_angle) * 180 / math.pi
    if sin_angle < 0 and cos_angle < 0:
        return 180 - math.asin(sin_angle) * 180 / math.pi
    if sin_angle < 0 and cos_angle > 0:
        return 360 - math.acos(cos_angle) * 180 / math.pi

def getJulian(filename, time):
    f = open(filename, "r")
    arr = []
    numlines = 0
    h = []
    for line in f.readlines():
        numlines += 1
        arr.append(line.split())
    f.close()
    arr = np.array(arr)
    for k in range(numlines):
        if time == arr[k,0]:
            obs = arr[k,1:].astype(np.float)
            return obs[6]

def getExpectedLPP(filename, time):
    f = open(filename, "r")
    arr = []
    numlines = 0
    for line in f.readlines():
        numlines += 1
        arr.append(line.split())
    f.close()
    arr = np.array(arr)
    for k in range(numlines):
        if time == arr[k,0]:
            obs = arr[k,1:].astype(np.float)
            return obs[7]

def angularMomentum(filename, time):
    f = open(filename, "r")
    arr = []
    numlines = 0
    h = []
    for line in f.readlines():
        numlines += 1
        arr.append(line.split())
    f.close()
    arr = np.array(arr)
    for k in range(numlines):
        if time == arr[k,0]:
            obs = arr[k,1:].astype(np.float)
            r = [obs[0], obs[1], obs[2]]
            rdot = [obs[3], obs[4], obs[5]]
            h = np.cross(r, rdot) * 365.2563835 / 2 / math.pi
            for k in range(len(h)):
                h[k] = round(h[k], 6)
            break
    return h

def semimajorAxis(filename, time):
    f = open(filename, "r")
    arr = []
    numlines = 0
    k = 0.01720209895
    mu = 1
    h = []
    a = 0
    for line in f.readlines():
        numlines += 1
        arr.append(line.split())
    f.close()
    arr = np.array(arr)
    for k in range(numlines):
        if time == arr[k,0]:
            obs = arr[k,1:].astype(np.float) 
            r = [obs[0], obs[1], obs[2]]
            rdot = [obs[3], obs[4], obs[5]]
            rdot = np.array(rdot) * 365.2563835 / 2 / math.pi
            a = (2 / magnitude(r) - magnitude(rdot) ** 2 / mu) ** -1
            return a


def eccentricity(filename, time):
    f = open(filename, "r")
    arr = []
    numlines = 0
    k = 0.01720209895
    mu = 1
    h = []
    for line in f.readlines():
        numlines += 1
        arr.append(line.split())
    f.close()
    arr = np.array(arr)
    for k in range(numlines):
        if time == arr[k,0]:
            obs = arr[k,1:].astype(np.float) 
            r = [obs[0], obs[1], obs[2]]
            rdot = [obs[3], obs[4], obs[5]] 
            rdot = np.array(rdot) * 365.2563835 / 2 / math.pi
            a = semimajorAxis(filename, time)
            h = angularMomentum(filename, time)
            #print(magnitude(h))
            #print(k ** 2)
            #print("value of h", h)
            #print(a)
            #print(1 - (magnitude(h) ** 2) / (mu * a))
            e = (1 - (magnitude(h) ** 2) / (mu * a)) ** 0.5
            #print(e)
            break
    return e

def inclination(filename, time):
    z = [0, 0, 1]
    h = angularMomentum(filename, time)
    return math.acos(np.dot(h, z) / magnitude(h)) * 180 / math.pi

def LoAN(filename, time):
    sin_i = math.sin(inclination(filename, time) * math.pi / 180)
    #print(sin_i)
    x = [1, 0, 0]
    y = [0, 1, 0]
    h = angularMomentum(filename, time)
    sin_omega = np.dot(h, x) / (magnitude(h) * sin_i)
    #print(sin_omega)
    cos_omega = -1 * np.dot(h, y) / (magnitude(h) * sin_i)
    #print(cos_omega)
    return quadrantCheck(cos_omega, sin_omega)

def AoP(filename, time):
    f = open(filename, "r")
    arr = []
    numlines = 0
    for line in f.readlines():
        numlines += 1
        arr.append(line.split())
    f.close()
    arr = np.array(arr)
    for k in range(numlines):
        if time == arr[k,0]:
            obs = arr[k,1:].astype(np.float)
            r = [obs[0], obs[1], obs[2]]
            rdot = [obs[3], obs[4], obs[5]]
            rdot = np.array(rdot) * 365.2563835 / 2 / math.pi
            break
    a = semimajorAxis(filename, time)
    cos_U = (r[0] * math.cos(LoAN(filename, time) * math.pi / 180) + r[1] * math.sin(LoAN(filename, time) * math.pi / 180)) / (magnitude(r))
    sin_U = r[2] / magnitude(r) / math.sin(inclination(filename, time) * math.pi / 180)
    cos_nu = (a * (1 - eccentricity(filename, time) ** 2) / magnitude(r) - 1) / eccentricity(filename, time)
    sin_nu = (a * (1 - eccentricity(filename, time) ** 2) * np.dot(r, rdot) / magnitude(r) / magnitude(angularMomentum(filename, time))) / eccentricity(filename, time)
    #print(cos_U)
    #print(sin_U)
    return (quadrantCheck(cos_U, sin_U) - quadrantCheck(cos_nu, sin_nu)) % 360

def meanAnomaly(filename, time):
    e = eccentricity(filename, time)
    a = semimajorAxis(filename, time)
    f = open(filename, "r")
    arr = []
    numlines = 0
    for line in f.readlines():
        numlines += 1
        arr.append(line.split())
    f.close()
    arr = np.array(arr)
    for k in range(numlines):
        if time == arr[k,0]:
            obs = arr[k,1:].astype(np.float)
            r = [obs[0], obs[1], obs[2]]
    cos_E = (1 - magnitude(r) / a) / e
    E = math.acos(cos_E)
    M = (E - e * math.sin(E)) * 180 / math.pi
    return M

def LPP(filename, time):
    M = meanAnomaly(filename, time)
    a = semimajorAxis(filename, time)
    k = 0.01720209895
    j = getJulian(filename, time)
    T = 2 * math.pi * a ** 1.5 / k
    lpp = j - M / 360 * T
    return lpp

def errorLPP(filename, time):
    lpp = LPP(filename, time)
    expectedLPP = getExpectedLPP(filename, time)
    return 100 * (expectedLPP - lpp) / lpp

def getSunVector(filename, time):
    f = open(filename, "r")
    arr = []
    numlines = 0
    for line in f.readlines():
        numlines += 1
        arr.append(line.split())
    f.close()
    arr = np.array(arr)
    for k in range(numlines):
        if time == arr[k,0]:
            obs = arr[k,1:].astype(np.float)
            R = [obs[0], obs[1], obs[2]]
    return R