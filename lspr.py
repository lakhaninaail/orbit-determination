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

def ra_deg_hms(deg):
    deg *= 180 / math.pi
    negative = math.copysign(1, deg)
    deg = abs(deg)
    hours_decimal = deg * 24 / 360
    hours = int(hours_decimal // 1)
    minutes_decimal = 60 * (hours_decimal - hours)
    minutes = int(minutes_decimal // 1)
    seconds = round(60 * (minutes_decimal - minutes), 2)
    if deg < 0:
        pm = "-"
    else:
        pm = ""
    return "{}{}:{}:{}".format(pm, '%02d' % hours, '%02d' % minutes, seconds) # '%02d' % from Stack Overflow

def dec_deg_dms(deg):
    deg *= 180 / math.pi
    negative = math.copysign(1, deg)
    deg = abs(deg)
    degree = int(deg // 1)
    minute = int(60 * (deg - degree) // 1)
    second = round(60 * (60 * (deg - degree) - 60 * (deg - degree) // 1), 1)
    if deg < 0:
        pm = "-"
    else:
        pm = "+"
    return "{}{}:{}:{}".format(pm, '%02d' % degree, '%02d' % minute, second)

def lspr(filename, c_x, c_y):
    f = open(filename, "r")
    arr = []
    numlines = 0
    for line in f.readlines():
        numlines += 1
        arr.append(line.split())
    f.close()
    arr = np.array(arr)

    x = arr[:,0].astype(np.float)
    y = arr[:,1].astype(np.float)

    ra = []
    dec = []
    for k in range(numlines):
        ra.append(convert_angle_ra(arr[k,2], True))
        dec.append(convert_angle_dec(arr[k,3], True))

    coord_matrix = [[numlines, np.sum(x), np.sum(y)],
                    [np.sum(x), np.dot(x, x), np.dot(x, y)],
                    [np.sum(y), np.dot(x, y), np.dot(y, y)]]
    inv = np.linalg.inv(coord_matrix)
    ra_coords = [[np.sum(ra)], [np.dot(ra, x)], [np.dot(ra, y)]]
    ra_pc = np.dot(inv, ra_coords)

    dec_coords = [[np.sum(dec)], [np.dot(dec, x)], [np.dot(dec, y)]]
    dec_pc = np.dot(inv, dec_coords)


    sum_sig_ra = 0
    for k in range(numlines):
        sum_sig_ra += (ra[k] - ra_pc[0,0] - ra_pc[1,0] * x[k] - ra_pc[2,0] * y[k]) ** 2
    sig_ra = ((sum_sig_ra) / (numlines - 3)) ** 0.5

    sum_sig_dec = 0
    for k in range(numlines):
        sum_sig_dec += (dec[k] - dec_pc[0,0] - dec_pc[1,0] * x[k] - dec_pc[2,0] * y[k]) ** 2
    sig_dec = ((sum_sig_dec) / (numlines - 3)) ** 0.5

    print("*" * 15)
    print("plate constants")
    print("*" * 15)
    print("b1: {} deg".format(ra_pc[0,0] * 180 / math.pi))
    print("b2: {} deg".format(dec_pc[0,0] * 180 / math.pi))
    print("a11: {} deg/pix".format(ra_pc[1,0] * 180 / math.pi))
    print("a12: {} deg/pix".format(ra_pc[2,0] * 180 / math.pi))
    print("a21: {} deg/pix".format(dec_pc[1,0] * 180 / math.pi))
    print("a22: {} deg/pix".format(dec_pc[2,0] * 180 / math.pi))
    print("*" * 11)
    print("uncertainty")
    print("*" * 11)
    print(" RA: {} arcsec".format(round(sig_ra * 3600 * 180 / math.pi, 2)))
    print("DEC: {} arcsec".format(round(sig_dec * 3600 * 180 / math.pi, 2)))
    print("*********************************")
    print("astrometry for")
    print("(x,y) = ({},{})".format(c_x, c_y))
    print("*********************************")
    print(" RA = {}".format(ra_deg_hms(ra_pc[0,0] + ra_pc[1,0] * x[k] + ra_pc[2,0] * y[k])))
    print("Dec = {}".format(dec_deg_dms(dec_pc[0,0] + dec_pc[1,0] * x[k] + dec_pc[2,0] * y[k])))

file = input("Enter file name of input data: ")
x_c = input("Enter x-coordinate of centroid: ")
y_c = input("Enter y-coordinate of centroid: ")
lspr(file, x_c, y_c)

#lspr("LSPRtestinput2.txt", 484.35, 382.62)