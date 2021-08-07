# SSP 2020 Python pset 2
# Naail Lakhani
# 6/27/2020

#from math import asin, acos, pi
import math

# a function to determine the quadrant of an angle based on its sine and cosine (in radians)
# returns the angle in the correct quadrant (in radians)
def findQuadrant(sine, cosine):
    if cosine > 0 and sine > 0: #1
        return math.asin(sine)

    if cosine < 0 and sine > 0: #2
        return math.acos(cosine)

    if cosine < 0 and sine < 0: #3
        return math.pi - math.asin(sine)

    if cosine > 0 and sine < 0: #4
        return 2*math.pi + math.asin(sine)

# a function that given the values (in decimal degrees) of two sides and the included angle of a spheical triangle
# returns the values of the remaining side and two angles (in decimal degrees)
def SAS(a, B, c):
    a *= math.pi / 180
    B *= math.pi / 180
    c *= math.pi / 180
    b = math.acos(math.cos(a) * math.cos(c) + math.sin(a) * math.sin(c) * math.cos(B))
    #b = findQuadrant(math.sin(b), math.cos(b))
    #b *= 180.0 / math.pi
    A = math.asin(math.sin(a) * math.sin(B) / math.sin(b))
    C = math.asin(math.sin(c) * math.sin(B) / math.sin(b))
    side_b = findQuadrant(math.sin(b), math.cos(b))
    #side_b = math.atan2(math.sin(b), math.cos(b))
    angle_A = findQuadrant(math.sin(A), math.cos(A))
    angle_C = findQuadrant(math.sin(C), math.cos(C))

    side_b *= 180 / math.pi
    angle_A *= 180 / math.pi
    angle_C *= 180 / math.pi
    return side_b, angle_A, angle_C

# DO NOT REMOVE OR MODIFY THIS CODE
s3, a1, a2 = SAS(106, 114, 42)
if abs(s3 - 117.804) > 1e-3 or abs(a1 - 83.11) > 1e-3 or abs(a2 - 43.715) > 1e-3:
    print("SAS function INCORRECT, expected (117.804, 83.11, 43.715), but got", (s3, a1, a2))
else:
    print("SAS function CORRECT")
