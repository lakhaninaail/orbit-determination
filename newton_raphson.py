import math

def newtonRaphson(M, e, tol):
    E_predict = M
    f = M - (E_predict - e * math.sin(E_predict))
    iterations = 0
    while(abs(f) >= tol):
        iterations += 1
        #f = M - (E_predict - e * math.sin(E_predict))
        E_predict = E_predict - (M - (E_predict - e * math.sin(E_predict))) / (e * math.cos(E_predict) - 1)
        f = M - (E_predict - e * math.sin(E_predict))
    return iterations, E_predict

M = float(input("Enter Mean Anomaly: "))
e = float(input("Enter Eccentricity: "))
tol = float(input("Enter Tolerance: "))

iterations, E = newtonRaphson(M, e, tol)

print("E = ",format(E))
print("Convergence Parameter = ",format(tol))
print("Number of Iterations = ",format(iterations))