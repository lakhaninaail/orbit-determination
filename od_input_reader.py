import numpy as np
import math

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
    obs1 = 1
    obs2 = 2
    obs3 = 3
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

o1, o2, o3 = fileReader("2020test (1).txt")
print(o1[0])
print(o2[1])
print(o3[2])