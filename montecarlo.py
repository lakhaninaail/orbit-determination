# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 17:15:49 2020

@author: lei_m
"""

import math
import numpy as np
from LeiOD1SEL import SEL
from fgfunction import FGfunctions
from LeiOD2 import orbitalElements
from sspleastplatesquaresreduction import convert_angleRA, convert_angleDEC
#from LeiOD3 import meanAnamoly1 
#from NewtonRaphson import newtonRas
np.set_printoptions(precision=12)
#take in raw inputs from file w/RA, dec in hr;min;sec and changes into arrays of floats
#makes everything into arrays to be passed into OD function

def initialize(name):
    file = open(name, "r")  #name should be input file name, "2020test.txt"
    count = 1
    t = np.array([0.,0.,0.])
    ra = np.array([0.,0.,0.])
    dec = np.array([0.,0.,0.])
    R = np.array([np.array([0., 0., 0.]),
                  np.array([0., 0., 0.]),
                  np.array([0., 0., 0.])])
    
        #the following code works, but by default uses the first 3 lines in input file
    for line in file:
        if(count<4): #REMOVE THIS LATER!!!!! WHEN YOU CAN CHOOSE TO USE 3 DIF OBSERV
            data = line.split()
            t[count-1] = data[0]
            tempRA = data[1]  #puts the hr:min:sec form for Ra and dec into a temp variable to be parsed and split
            tempDEC = data[2]
            R[count-1][0] = data[3] 
            R[count-1][1] = data[4]
            R[count-1][2] = data[5]
            colonRA = tempRA.split(":")
            colonDEC = tempDEC.split(":")
            ra[count-1] = convert_angleRA(float(colonRA[0]),float(colonRA[1]),float(colonRA[2]))
            dec[count-1] = convert_angleDEC(float(colonDEC[0]),float(colonDEC[1]),float(colonDEC[2]))
            count = count+1
        
    """
    #the following code works, but by default uses the first 3 lines in input file
    for line in file:
        if(count<4): #REMOVE THIS LATER!!!!! WHEN YOU CAN CHOOSE TO USE 3 DIF OBSERV
            data = line.split()
            t[count-1] = data[0]
            tempRA = data[1]  #puts the hr:min:sec form for Ra and dec into a temp variable to be parsed and split
            tempDEC = data[2]
            R[count-1][0] = data[3] 
            R[count-1][1] = data[4]
            R[count-1][2] = data[5]
            colonRA = tempRA.split(":")
            colonDEC = tempDEC.split(":")
            ra[count-1] = convert_angleRA(float(colonRA[0]),float(colonRA[1]),float(colonRA[2]))
            dec[count-1] = convert_angleDEC(float(colonDEC[0]),float(colonDEC[1]),float(colonDEC[2]))
            count = count+1
            #all arrays should be formed when exiting this for loop
    """
    return ra, dec,t,R


#pass in arrays for ra, dec, t, and R
#p is rho, R is sun vector, tau is difference in time, D is constants
#tol is amount of tolerance allowed, flag is a string to say if you use the f&g functins, 3 term series, or 4 term series
def OD(ra,dec,t,R, tol):
    file1 = open("output.txt", "w")
    file2 = open("output2.txt", "w")
    phat = np.array([np.array([0., 0., 0.]),
                    np.array([0., 0., 0.]),
                    np.array([0., 0., 0.])]) #0th terms is for observation 1, 1th for obsv.2, 2th for obsv.3||each obsv. has 3 values of phat
    #constants
    t0 = t[0]
    t1 = t[1]
    t2 = t[2]
    for m in range(3):
        ra[m]=ra[m]*math.pi/180
        dec[m]=dec[m]*math.pi/180

    for i in range(3):
        phat[i][0] = math.cos(ra[i])*math.cos(dec[i])
        phat[i][1] = math.sin(ra[i])*math.cos(dec[i])
        phat[i][2] = math.sin(dec[i])
  
    D0 = np.dot(phat[0],np.cross(phat[1], phat[2]))
    D = np.array([np.array([D0,D0,D0]),
                  np.array([0., 0., 0.]),
                  np.array([0., 0., 0.]),
                  np.array([0., 0., 0.])])
                 
    #D0 is a constant for all observations
    for j in range(3):
        D[1][j] = np.dot(np.cross(R[j],phat[1]),phat[2])
        D[2][j] = np.dot(np.cross(phat[0],R[j]),phat[2])
        D[3][j] = np.dot(phat[0],np.cross(phat[1], R[j]))
    
    k = 0.01720209895 #Gaussian gravitational constant
    tau = np.array([0.,0.,0.]) #0th position item is tau, 1st position is tau_1, 2nd position is tau_3
    tau[1] = k*(t[0]-t[1]) #tau_1
    tau[2] = k*(t[2]-t[1]) #tau_3
    tau[0] = tau[2]-tau[1] #tau
    #print(tau)
    taus = np.array([tau[1],tau[2],tau[0]])
    Ds = np.array([D0, D[2][0],D[2][1],D[2][2]])
    r2array, pGuess = SEL(taus, R[1], phat[1], Ds)  #returns an estimate for r2 and pGuess vectors(position+rho)
    """
    #r2 is a magnitude
    for n in range(0,len(r2array)):
        A1 = tau[2]/tau[0]
        A3 = -tau[1]/tau[0]
        B1 = A1*(tau[0]**2-tau[2]**2)/6
        B3 = A3*(tau[0]**2-tau[1]**2)/6
        A = -(A1*D[2][0]-D[2][1]+A3*D[2][2])/D[0][0]
        B = -(B1*D[2][0]+B3*D[2][2])/D[0][0]
        p = A + B/(r2array[n]**3)
        if(p>0):
            r2 = r2array[n] #change later!! Should ask user to specify which specific root to use 
    r2mag = r2
    """
    #new attempt at taking in user input and can deal with multiple r2 values
    r2good = []
    for n in range(0,len(r2array)):
        A1 = tau[2]/tau[0]
        A3 = -tau[1]/tau[0]
        B1 = A1*(tau[0]**2-tau[2]**2)/6
        B3 = A3*(tau[0]**2-tau[1]**2)/6
        A = -(A1*D[2][0]-D[2][1]+A3*D[2][2])/D[0][0]
        B = -(B1*D[2][0]+B3*D[2][2])/D[0][0]
        p = A + B/(r2array[n]**3)
        if(p>0):
            r2good.append(r2array[n]) #change later!! Should ask user to specify which specific root to use 

    if(len(r2good)>1):
        r2 = r2good[-1]
        r2mag = r2
    else: 
        r2 = r2good[0]
        r2mag = r2
    
    
    #initial cycle using 2 term expression for f & g series 
    f1 = 1-u*(tau[1]**2)/(2*pow(r2mag,3))
    f3 = 1-u*(tau[2]**2)/(2*pow(r2mag,3))
    g1 = tau[1]
    g3 = tau[2]
    d1 = -f3/(f1*g3-f3*g1)
    d3 = f1/(f1*g3-f3*g1)
    c1 = g3/(f1*g3-f3*g1)
    c2 = -1
    c3 = -g1/(f1*g3-f3*g1)
    p = np.array([0., 0., 0.])
    p[0] = (c1*D[1][0]+c2*D[1][1]+c3*D[1][2])/(c1*D[0][0])
    p[1] = (c1*D[2][0]+c2*D[2][1]+c3*D[2][2])/(c2*D[0][0])
    p[2] = (c1*D[3][0]+c2*D[3][1]+c3*D[3][2])/(c3*D[0][0])
    #F AND G's CORRECT AT THIS POINT, EXACTLY CHECKED
    #INPUTs ALL CORRECT
    #p's all correct
    """
    print("HI")
    print(c1,c2,c3)
    print("Ds")
    print(D)
    print(p)
    """
    
    r1 = p[0]*phat[0] - R[0][:] # CHECK TO SEE IF r1,r2,r3 are VECTORS!!
    r2new = p[1]*phat[1] - R[1][:]
    r2magnew = pow(r2new[0]**2 + r2new[1]**2 + r2new[2]**2, 0.5)
    r3 = p[2]*phat[2] - R[2][:]
    #print("r2new")
    #print(p, phat, r2new,r1,r3)
    #print(r2new,r2magnew,r2,r2mag)
    r2dot = d1*r1 + d3*r3
    #correct for light travel time below!
    #print(t)
    for k in range(3):
        t[k] = t[k] - p[k]/cAU
        #print("LIGHT")
    #recalculuate taus to correct
    k = 0.01720209895
    tau[1] = k*(t[0]-t[1]) #tau_1
    tau[2] = k*(t[2]-t[1]) #tau_3
    tau[0] = tau[2]-tau[1] #tau
    
    variable = True
    #test for convergence?
    if(abs(r2magnew-r2mag)<tol):
        variable = False
    r2 = r2new
    r2mag = r2magnew
    iterations = 1
    """
    print("position vectors")
    print(r1,r2,r3,r2dot,r2magnew)
    #print("ÿeehaw")
    #print(r2mag, r2magnew, r2)
    print("tau")
    print(t)
    print(tau)
    print(f1,g1,f3,g3)
    #add in other taus
    #tau[0] = 0.3770495384812921
    #tau[1] = -0.10432058089523047
    #tau[2] = 0.2727289575860616
    """
    while(variable):
        """
        uSeries = u/r2mag**3
        #print(r2mag,uSeries)
        z = np.dot(r2,r2dot)/r2mag**2
        q = np.dot(r2dot,r2dot)/r2mag**2-uSeries
        print(uSeries,z,q)
        """
        #print(iterations)
        """
        f1 = 1 - u *(tau[1]**2) /(2*pow(r2mag,3)) + u * np.dot(r2, r2dot) * (tau[1]**3)/(2*pow(r2mag,5)) + (3*uSeries*q-15*uSeries*z*z+uSeries**2)*pow(tau[1],4)/24
        g1 = (tau[1]) - u * (tau[1]**3)/(6*pow(r2mag,3)) + 6*uSeries*z*pow(tau[1],4)/24
        f3 = 1 - u *(tau[2]**2) /(2*pow(r2mag,3)) + u * np.dot(r2, r2dot) * (tau[2]**3)/(2*pow(r2mag,5)) + (3*uSeries*q-15*uSeries*z*z+uSeries**2)*pow(tau[2],4)/24
        g3 = (tau[2]) - u * (tau[2]**3)/(6*pow(r2mag,3)) + 6*uSeries*z*pow(tau[2],4)/24
        """
        
        f1,g1,f3,g3 = FGfunctions(tau[1], tau[2], r2, r2dot, "4series")
        #FGfunctions(tau[1], tau[2], r2, r2dot, "4series")) #maybe try r2new?
        #print(tau[1],tau[2],r2,r2dot,"4series")
        #print(FGfunctions(-0.10431549 ,0.27271306 ,[0.08270936, -1.25067987, -0.52047478] , [ 0.8667031, -0.44634935, 0.09702106] , "4series"))
        d1 = -f3/(f1*g3-f3*g1)
        d3 = f1/(f1*g3-f3*g1)
        c1 = g3/(f1*g3-f3*g1)
        c2 = -1
        c3 = -g1/(f1*g3-f3*g1)
        #print("Cs")
        #print(c1,c2,c3)
        #print(d1,d3)
        pnew = np.array([0., 0., 0.])
        pnew[0] = (c1*D[1][0]+c2*D[1][1]+c3*D[1][2])/(c1*D[0][0])
        pnew[1] = (c1*D[2][0]+c2*D[2][1]+c3*D[2][2])/(c2*D[0][0])
        pnew[2] = (c1*D[3][0]+c2*D[3][1]+c3*D[3][2])/(c3*D[0][0])
        
        #pmagnew = pow(p[0]**2 + p[1]**2 + p[2]**2, 0.5)
        #print(pnew)
        r1 = pnew[0]*phat[0] - R[0][:] # CHECK TO SEE IF r1,r2,r3 are VECTORS!!
        r2new = pnew[1]*phat[1] - R[1][:]
        r2magnew = pow(r2new[0]**2 + r2new[1]**2 + r2new[2]**2, 0.5)
        r3 = pnew[2]*phat[2] - R[2][:]
        
        r2dot = d1*r1 + d3*r3
        #correct for light travel time below!
        for k in range(3):
            t[k] = t[k] - pnew[k]/cAU
        #print("here")
        #file1.write("next"+ str(p))
        #file1.write("\n")
        file2.write("next"+ str(pnew[1]-p[1]) +"split"+ str(pnew[1]*24*3600/cAU))
        file2.write("\n")
        
        #print(p[1]/cAU)
        #recalculuate taus to correct
        k = 0.01720209895
        tau[1] = k*(t[0]-t[1]) #tau_1
        tau[2] = k*(t[2]-t[1]) #tau_3
        tau[0] = tau[2]-tau[1]
        t[0]=t0
        t[1]=t1
        t[2]=t2
        #print("almost there")
        #print(r2, r2new, r2magnew, r2mag, r2dot)
        """
        if(abs(r2magnew-r2mag)<tol):
            break
        r2 = r2new
        r2mag = r2magnew
        iterations = iterations+1
        """
        if(abs(pnew[1]-p[1])<tol):
            break
        r2 = r2new
        r2mag = r2magnew
        p = pnew
        iterations = iterations+1
        
        #update all the variables???
    """
    r2 = pnew[1]*phat[1] - R[1][:]
    r2magnew = pow(r2[0]**2 + r2[1]**2 + r2[2]**2, 0.5)
    r1 = pnew[0]*phat[0] - R[0][:]
    r3 = pnew[2]*phat[2] - R[2][:]
    r2dot = d1*r1 + d3*r3
    """
    print("range to middle observation is " + str(pnew[1]))
    rotation = np.array([[1, 0, 0],[0, math.cos(eps), -math.sin(eps)], [0, math.sin(eps), math.cos(eps)]])
    r2 = np.dot(r2, rotation)
    r2dot = np.dot(r2dot, rotation)
    #file1.write(r2) 
    #file1.write(r2dot)     
    #r2vector = r2
    #h, myPosition, myVelocity = angularMomentum("output.txt")
    #posMag,a,e,i,longA,w,m_Actual = orbitalElements("LeiInput.txt", r2, r2dot)
    print(r2,r2dot, pnew)
    posMag,a,e,i,longA,w,m_Actual, v = orbitalElements("LeiInput.txt", r2,r2dot)
    #jDateLast = 2457852.08061
    #m_Actual, timePereh = meanAnamoly1(a, t[1], jDateLast)
    #r2vectormag = pow(r2vector[0]**2 + r2vector[1]**2 + r2vector[2]**2, 0.5)
    #FIND MEAN ANOMALY
    """
    E = math.cos((1-r2mag/a)/e)
    E = E * 180/math.pi #CONVERT E TO DEGREES AFTER CALLING NEWTON-RASPHON
    v = v * 180/math.pi
    if(v>=0 and v<180):
        if(E>=180 and E<360):
            E = 360 - E
    else:
        if(E<=180 and E>0):
            E = 360 - E
    E = E * math.pi/180 #convert back to radians for trig stuff
    m_Actual = E - e* math.sin(E)
    """
    v = v * 180/math.pi
    if(v>=0 and v<180):
        E = math.acos((1-r2mag/a)/e)
        m_Actual = E - e* math.sin(E)
    else:
        E = 2*math.pi - math.acos((1-r2mag/a)/e)
        m_Actual = E - e* math.sin(E)
    #print(a,e,i,longA,w,m_Actual, v,E)
    P = pow(a**3, 0.5)
    #timePereh = t[1] - m_Actual*P/2*math.pi + P * 365.256363004
    timePereh = t[1] + (2*math.pi-m_Actual)*P * 365.256363004/(2*math.pi)
    #convert to degrees:
    E = E * 180/math.pi
    m_Actual = m_Actual * 180/math.pi
    return a,e,i,longA,w,E, m_Actual, t[1], timePereh, P
    #return 1,2,3,4,5,6

    """"
    f1,g1,f3,g3 = FGfunctions(tau[1],tau[2],r2,r2dot,flag)
    d1 = -f3/(f1*g3-f3*g1)
    d3 = f1/(f1*g3-f3*g1)
    c1 = g3/(f1*g3-f3*g1)
    c3 = -g1/(f1*g3-f3*g1)
    """
    
r2vector = np.array([0.,0.,0.])
k = 0.01720209895 #Gaussian gravitational constant
cAU = 173.144632674240 #speed of light in au/(mean solar)day
eps = math.radians(23.4366) #Earth's obliquity
u = 1 #where u is equal to mu variable, which is 1 when in units of AU and Gaussian days
#ra, dec,t,R = initialize("2020test.txt")
ra, dec,t,R = initialize("montecarloInput.txt")
a,e,i,longA,w,E, m_Actual, t[1], timePereh, P = OD(ra,dec,t,R, 1e-12) #CHANGE THIS BACK TO 1E-12 LATER!!!!!!!!!!!!
#m_Actual = m_Actual*math.pi/180 #CONVERT M TO RADIANS BEFORE CALLING NEWTONRASPHON
#E = newtonRas(m_Actual,e,1e-12,m_Actual)
#jDateJuly25 = 2459056.16667
#m_July25, timePereh = meanAnamoly1(a, jDateJuly25, jDateLast)

print("The following is a,e,i,w,longA,mean anomaly, E, julian date, time of perehelion, P(period) in exactly that order:")
print(a,e,i,w,longA,m_Actual, E, t[1], timePereh, P)


def ODmontecarlo(name, uncertaintyRA, uncertaintyDEC, cycles):
    centerRA, centerDEC,t,R = initialize(name) # ra, dec from initialize function!!!
    a = np.zeros(cycles)
    e = np.zeros(cycles)
    i = np.zeros(cycles)
    longA = np.zeros(cycles)
    w = np.zeros(cycles)
    m_Actual = np.zeros(cycles)
    step = 0
    while(step<cycles):
        RA = np.array([0.,0.,0.])
        DEC = np.array([0.,0.,0.])
        for k in range(0,3):
            #print("JI")
            #print(np.random.normal(centerRA, uncertaintyRA))
            tempRA = np.random.normal(centerRA[k], uncertaintyRA[k])
            tempDEC = np.random.normal(centerDEC[k], uncertaintyDEC[k])
            #print(tempRA)
            RA[k] = tempRA
            DEC[k] = tempDEC
        print("f")
        a[step],e[step],i[step],longA[step],w[step],E,m_Actual[step], t[1], timePereh, P = OD(RA,DEC,t,R, 1e-12) #assuming t and R are accurate, no error
        step = step+1  
        #finds sum in prep to calculate averages
    print("A")
    suma,sume,sumi,sumlongA, sumw, summ_Actual = 0.,0.,0.,0.,0.,0.
    for y in range(0, cycles): 
        suma = suma+a[y]
        sume = sume+e[y]
        sumi = sumi+i[y]
        sumlongA = sumlongA + longA[y]
        sumw = sumw + w[y]
        summ_Actual = summ_Actual + m_Actual[y]
        #calculates averages
    averagea = suma/cycles
    averagee = sume/cycles
    averagei = sumi/cycles
    averagelongA = sumlongA/cycles
    averagew = sumw/cycles
    averagem_Actual = summ_Actual/cycles    
        #find sums for standard deviation
    sumdeviationa = 0
    sumdeviatione = 0
    sumdeviationi = 0
    sumdeviationlongA = 0
    sumdeviationw = 0
    sumdeviationm_Actual = 0
    print("B")
    for y in range(0, cycles):
        sumdeviationa = sumdeviationa + (a[y]-averagea)**2
        sumdeviatione = sumdeviatione + (e[y]-averagee)**2
        sumdeviationi = sumdeviationi + (i[y]-averagei)**2
        sumdeviationlongA = sumdeviationlongA + (longA[y]-averagelongA)**2
        sumdeviationw = sumdeviationw + (w[y]-averagew)**2
        sumdeviationm_Actual = sumdeviationm_Actual + (m_Actual[y]-averagem_Actual)**2
        #sd stands for standard deviation
    sda = (sumdeviationa/(cycles-1))**0.5
    sde= (sumdeviatione/(cycles-1))**0.5
    sdi = (sumdeviationi/(cycles-1))**0.5
    sdlongA = (sumdeviationlongA/(cycles-1))**0.5
    sdw = (sumdeviationw/(cycles-1))**0.5
    sdm_Actual = (sumdeviationm_Actual/(cycles-1))**0.5
    print("C")
    print(averagea, averagee, averagei, averagelongA, averagew,averagem_Actual, sda, sde, sdi, sdlongA, sdw, sdm_Actual)
    return averagea, averagee, averagei, averagelongA, averagew,averagem_Actual, sda, sde, sdi, sdlongA, sdw, sdm_Actual

uncertaintyRA = (0.00016125197856363312, 0.00010593071102531975, 0.0001582285458841291)
uncertaintyDEC = (9.657959724260052e-05, 0.00011440237895469723, 0.00015452302921633567)
print("LLLL")
print(ODmontecarlo("montecarloInput.txt",uncertaintyRA ,uncertaintyDEC, 1000))