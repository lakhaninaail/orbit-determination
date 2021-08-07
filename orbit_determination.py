#from odlib import *
import numpy as np
import math

obs1 = 1 #ENTER WHICH OBSERVATIONS YOU WANT TO USE HERE -- USE LINE NUMBERS ON LakhaniInput.txt FILE, NOT INDEX IN ARRAY
obs2 = 2
obs3 = 3

num_root = 1 #ENTER NUMBER OF ROOT YOU WISH TO USE HERE (first, second, third, etc.)

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

k = 0.01720209895 #Gaussian gravitational constant
cAU = 173.144632674240 #speed of light in au/(mean solar)day
eps = math.radians(23.4366) #Earth's obliquity
mu = 1

def convert_angle_dec(string, radians):
    dec_val = string.split(":")
    for k in range(len(dec_val)):
        dec_val[k] = float(dec_val[k])
    degrees = dec_val[0]
    minutes = dec_val[1]
    seconds = dec_val[2]  
    negative = math.copysign(1.0, degrees)
    angle = degrees + negative * minutes / 60 + negative * seconds / 3600
    if radians:
        angle *= math.pi / 180
    return angle

def convert_angle_ra(string, radians):
    ra_val = string.split(":")
    for k in range(len(ra_val)):
        ra_val[k] = float(ra_val[k])
    hours = ra_val[0]
    minutes = ra_val[1]
    seconds = ra_val[2]  
    negative = math.copysign(1.0, hours)
    angle = 360 * (hours + negative * minutes / 60 + negative * seconds / 3600) / 24
    if radians:
        angle *= math.pi / 180
    return angle

def fileReader(filename):
    f = open(filename, "r")
    arr = []
    numlines = 0
    for line in f.readlines():
        numlines += 1
        arr.append(line.split())
    f.close()
    arr = np.array(arr)
    for m in range(numlines):
        arr[m,1] = convert_angle_ra(arr[m,1], True)
        arr[m,2] = convert_angle_dec(arr[m,2], True)
    for k in range(numlines):
        obs = arr[k,:].astype(np.float)
    data1 = arr[obs1 - 1]
    data2 = arr[obs2 - 1]
    data3 = arr[obs3 - 1]
    return data1, data2, data3

o1, o2, o3 = fileReader("LakhaniInput2020.txt")
o1 = o1[:].astype(np.float)
o2 = o2[:].astype(np.float)
o3 = o3[:].astype(np.float)

t1 = float(o1[0])
t2 = float(o2[0])
t3 = float(o3[0])
tau = [k * (t3 - t2), k * (t1 - t2), k * (t3 - t1)] #(tau1, tau3, tau)
sun1 = [o1[3], o1[4], o1[5]] #read from input file
sun2 = [o2[3], o2[4], o2[5]] #read from input file
sun3 = [o3[3], o3[4], o3[5]] #read from input file

def get_rhohat(ra, dec):
    rhohat = [math.cos(ra) * math.cos(dec), math.sin(ra) * math.cos(dec), math.sin(dec)]
    return rhohat

rhohat1 = get_rhohat(o1[1], o1[2])
rhohat2 = get_rhohat(o2[1], o2[2])
rhohat3 = get_rhohat(o3[1], o3[2])

def getD(observation_number):
    if observation_number == 1:
        sun = sun1
    elif observation_number == 2:
        sun = sun2
    else:
        sun = sun3
    D1 = np.dot(np.cross(sun, rhohat2), rhohat3)
    D2 = np.dot(np.cross(rhohat1, sun), rhohat3)
    D3 = np.dot(rhohat1, np.cross(rhohat2, sun))
    D = [D1, D2, D3]
    return D

D0 = np.dot(rhohat1, np.cross(rhohat2, rhohat3))
D1j = getD(1)
D2j = getD(2)
D3j = getD(3)

def SEL():
    rhos = [] #range values for each real, positive root
    Tau = tau[2]
    A1 = tau[1] / Tau
    B1 = A1 * (Tau ** 2 - tau[1] ** 2) / 6
    A3 = -1 * tau[0] / Tau
    B3 = A3 * (Tau ** 2 - tau[0] ** 2) / 6
    A = (A1 * D2j[0] - D2j[1] + A3 * D2j[2]) / (-1 * D0)
    B = (B1 * D2j[0] + B3 * D2j[2]) / (-1 * D0)
    E = -2 * (np.dot(rhohat2, sun2))
    F = (sun2[0] ** 2 + sun2[1] ** 2 + sun2[2] ** 2)
    a = -1 * (A ** 2 + A * E + F)
    b = -1 * mu * (2 * A * B + B * E)
    c = -1 * mu ** 2 * B ** 2
    coeff = [c,0,0,b,0,0,a,0,1]
    roots = np.real(np.polynomial.polynomial.polyroots(coeff))
    positive_roots = roots[roots>0] # Daniel told me about this built in function in Python
    for elem in positive_roots:
        Rho = A + mu * B / elem ** 3
        rhos.append(Rho)
    return(positive_roots, rhos)

roots, rhos = SEL()
using_root = roots[num_root - 1] 

def fg(tau, r2, r2dot, flag):
    if flag == 0: #function
        delta_E = newtonRaphson(tau, e, tol, a, n)
        f = 1 - a * (1 - math.cos(delta_E)) / magnitude(r2)
        g = tau + (math.sin(delta_E) - delta_E) / n
    
    if flag == 2: #second order
        f = 1 - mu / (2 * r2 ** 3) * tau ** 2
        g = tau - mu * tau ** 3 / (6 * r2 ** 3)

    if flag == 3: #third order
        f = 1 - mu / (2 * magnitude(r2) ** 3) * tau ** 2 + mu * np.dot(r2, r2dot) / (2 * magnitude(r2) ** 5) * tau ** 3
        g = tau - mu * tau ** 3 / (6 * magnitude(r2) ** 3)
    
    if flag == 4: #fourth order
        u = 1 / magnitude(r2) ** 3
        z = np.dot(r2, r2dot) / magnitude(r2) ** 2
        q = np.dot(r2, r2dot) / magnitude(r2) ** 2 - u
        f = 1 - mu / (2 * magnitude(r2) ** 3) * tau ** 2 + mu * np.dot(r2, r2dot) / (2 * magnitude(r2) ** 5) * tau ** 3 + (3 * u * q - 15 * u * z ** 2 + u ** 2) / 24 * tau ** 4
        g = tau - mu * tau ** 3 / (6 * magnitude(r2) ** 3) + 6 * u * z * tau ** 4 / 24

    return f, g

def c_calculator(r2, r2dot, flag):
    f1, g1 = fg(tau[0], r2, r2dot, flag)
    f3, g3 = fg(tau[1], r2, r2dot, flag)
    c1 = g3 / (f1 * g3 - g1 * f3)
    c2 = -1
    c3 = -g1 / (f1 * g3 - g1 * f3)
    return c1, c2, c3

c_initial = c_calculator(using_root, 0, 2)
#print("c:",c_initial)

def get_rho(c, obs_num):
    if obs_num == 1:
        D = D1j
    elif obs_num == 2:
        D = D2j
    else:
        D = D3j
    rho = (c[0] * D[0] + c[1] * D[1] + c[2] * D[2]) / (D0 * c[obs_num - 1])
    return rho

rho1 = get_rho(c_initial, 1)
rho2 = get_rho(c_initial, 2)
rho3 = get_rho(c_initial, 3)
#print("rho2:",rho2)

def get_r(obs_num):
    if obs_num == 1:
        rho = rho1
        rhohat = rhohat1
        sun = sun1
    elif obs_num == 2:
        rho = rho2
        rhohat = rhohat2
        sun = sun2
    else:
        rho = rho3
        rhohat = rhohat3
        sun = sun3
    r = []
    for n in range(3):
        rho_vec = rho * rhohat[n]
        r.append(rho_vec - sun[n])
    return r

r1 = get_r(1)
r2 = get_r(2)
r3 = get_r(3)
#print("r2:",r2)

def central_velocity_vector(r2, r2dot, flag):
    f1, g1 = fg(tau[0], r2, r2dot, flag)
    f3, g3 = fg(tau[1], r2, r2dot, flag)
    r2dot = []
    d1 = -f3 / (f1 * g3 - f3 * g1)
    d3 = f1 / (f1 * g3 - f3 * g1)
    for k in range(3):
        r2dot.append(d1 * r1[k] + d3 * r3[k])
    return r2dot

#r1dot = central_velocity_vector(magnitude(r1), 0, 2)
r2dot = central_velocity_vector(magnitude(r2), 0, 2)
#r3dot = central_velocity_vector(magnitude(r3), 0, 2)

def new_time_obs(obs_num):
    if obs_num == 1:
        rho = rho1
        t_og = t1
    elif obs_num == 2:
        rho = rho2
        t_og = t2
    else:
        rho = rho3
        t_og = t3
    t = t_og - rho/ cAU
    return t

t1_updated = new_time_obs(1)
t2_updated = new_time_obs(2)
t3_updated = new_time_obs(3)

#a = (2 / magnitude(r2) - magnitude(r2dot) ** 2 / mu) ** -1
r2_mag = magnitude(r2)
r2dot_mag = magnitude(r2dot) / 365.2563835 * 2 * math.pi
a = (2 / r2_mag - r2dot_mag ** 2 / mu) ** -1
#a = 1.5
#print("a:",a)
n = (mu / a ** 3) ** 0.5
#print(n)

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
    return x #delta E

def get_r_new(obs_num, rho):
    if obs_num == 1:
        rhohat = rhohat1
        sun = sun1
    elif obs_num == 2:
        rhohat = rhohat2
        sun = sun2
    else:
        rhohat = rhohat3
        sun = sun3
    r = []
    for n in range(3):
        rho_vec = rho * rhohat[n]
        r.append(rho_vec - sun)
    return r

def new_time_obs_ITERATIVE(rho, t_og):
    t = t_og - rho/ cAU
    return t


def main():
    rho2_new = 0
    diff = rho2 - rho2_new
    #print(diff)
    r2_new = r2
    r2dot_new = r2dot
    t1_new = t1
    t2_new = t2
    t3_new = t3
    tau = [k * (t3_new - t2_new), k * (t1_new - t2_new), k * (t3_new - t1_new)] #(tau1, tau3, tau)
    diff = rho2 - rho2_new
    #print(abs(diff) >= tol)
    while(abs(diff.any()) >= 1E-12):
        old_rho = rho2_new
        print("PASS")
        #f, g = fg(tau, r2_new, r2dot_new, 4)
        f1, g1 = fg(tau[0], r2, r2dot, 3)
        f3, g3 = fg(tau[1], r2, r2dot, 3)
        #delta_E = newtonRaphson(tau, e, tol, a, n)
        c_new = c_calculator(r2_new, r2dot_new, 3)
        rho1_new = get_rho(c_new, 1)
        rho2_new = get_rho(c_new, 2)
        rho3_new = get_rho(c_new, 3)
        r1_new = get_r_new(1, rho1_new)
        r2_new = get_r_new(2, rho2_new)
        r3_new = get_r_new(3, rho3_new)
        #r1dot_new = central_velocity_vector(r1_new, 0, 4)
        r2dot_new = central_velocity_vector(r2_new, 0, 3)
        #r3dot_new = central_velocity_vector(r3_new, 0, 4)
        t1_new = new_time_obs_ITERATIVE(rho1_new, t1_new)
        t2_new = new_time_obs_ITERATIVE(rho2_new, t2_new)
        t3_new = new_time_obs_ITERATIVE(rho3_new, t3_new)
        tau = [k * (t3_new - t2_new), k * (t1_new - t2_new), k * (t3_new - t1_new)] #(tau1, tau3, tau)
        diff = old_rho - rho2_new
    return rho2_new, r2_new, r2dot_new

range_2, x2, v2 = main()
#print(main())

def ecliptic(vec):
    new_vec = [[vec[0]], [vec[1]], [vec[2]]]
    new_vec = np.array(new_vec)
    arr = [[1, 0, 0], [0, math.cos(eps), math.sin(eps)], [0, -math.sin(eps), math.cos(eps)]]
    arr = np.array(arr)
    ecliptic_coord = np.dot(arr, vec)
    return ecliptic_coord[0], ecliptic_coord[1], ecliptic_coord[2]

def angularMomentum():
    h = np.cross(x2, v2) * 365.2563835 / 2 / math.pi
    return h

h = angularMomentum()

def eccentricity():
    return (1 - magnitude(h) ** 2 / (mu * a)) ** 0.5

e = eccentricity()

def inclination():
    z = [0, 0, 1]
    return math.acos(np.dot(h, z) / magnitude(h)) * 180 / math.pi

i = inclination()

def LoAN():
    sin_i = math.sin(i * math.pi / 180)
    #print(sin_i)
    x = [1, 0, 0]
    y = [0, 1, 0]
    sin_omega = np.dot(h, x) / (magnitude(h) * sin_i)
    #print(sin_omega)
    cos_omega = -1 * np.dot(h, y) / (magnitude(h) * sin_i)
    #print(cos_omega)
    return quadrantCheck(cos_omega, sin_omega)

loan = LoAN()

def AoP():
    cos_U = (x2[0] * math.cos(loan * math.pi / 180) + x2[1] * math.sin(loan * math.pi / 180)) / (magnitude(x2))
    sin_U = x2[2] / magnitude(r) / math.sin(i * math.pi / 180)
    cos_nu = (a * (1 - e ** 2) / magnitude(x2) - 1) / e
    sin_nu = (a * (1 - e ** 2) * np.dot(x2, v2) / magnitude(x2) / magnitude(h)) / e
    #print(cos_U)
    #print(sin_U)
    return (quadrantCheck(cos_U, sin_U) - quadrantCheck(cos_nu, sin_nu)) % 360

aop = AoP()

def meanAnomaly():
    cos_E = (1 - magnitude(x2) / a) / e
    E = math.acos(cos_E)
    M = (E - e * math.sin(E)) * 180 / math.pi
    return M

M = meanAnomaly()

ecliptic_r2 = ecliptic(x2)
ecliptic_r2dot = ecliptic(v2)
print("range: ",range_2)
print("position: ",ecliptic_r2)
print("velocity: ",ecliptic_r2dot)
print("semimajor axis: ",a)
print("eccentricity: ",e)
print("orbit inclination: ", i)
print("LoAN: ",loan)
print("AoP: ",aop)
print("mean anomaly: ",M)