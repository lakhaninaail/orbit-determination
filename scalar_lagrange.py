import numpy as np

tau = [-0.1720209895, 0.1720209895, 0.344041979]
sun2 = [-0.2683394448727136, 0.8997620807182745, 0.3900022331276332]
rhohat2 =  [0.052719013914983195, -0.9121555187306237, 0.40643943610469035]
D = [0.0011797570297704812, 0.052586743761143424, 0.05848153743706686, 0.06274019190783499]

def SEL():
    rhos = [] #range values for each real, positive root
    mu = 1
    Tau = tau[2]
    A1 = tau[1] / Tau
    B1 = A1 * (Tau ** 2 - tau[1] ** 2) / 6
    A3 = -1 * tau[0] / Tau
    B3 = A3 * (Tau ** 2 - tau[0] ** 2) / 6
    A = (A1 * D[1] - D[2] + A3 * D[3]) / (-1 * D[0])
    B = (B1 * D[1] + B3 * D[3]) / (-1 * D[0])
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

#roots, rhos = SEL(tau,sun2,rhohat2,D)
roots, rhos = SEL()

print(roots)
print(rhos)