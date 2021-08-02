# Programming Assignment #3 / SSP
# Naail Lakhani
# 6/30/2020

import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.io import fits

def findCentroid(fits_file, target_x, target_y, radius):
    counts_array = fits.getdata(fits_file)
    sumsx = []
    sumsy = []
    small_array = counts_array[target_y - radius : target_y + radius + 1, target_x - radius : target_x + radius + 1]
    for x in range(len(small_array)):
        sumsx.append(np.sum(small_array[:,x]))
        
    totalx = 0
    for k in range(len(sumsx)):
        totalx += (k + target_x - radius) * sumsx[k]
    c_x = totalx / np.sum(sumsx)

    sum_sd_x = 0
    for p in range(len(small_array)):
        sum_sd_x += ((p + target_x - radius - c_x) / np.sum(sumsx)) ** 2 * sumsx[p]
    sdx = sum_sd_x ** 0.5

    for y in range(len(small_array)):
        sumsy.append(np.sum(small_array[y,:]))
    totaly = 0
    for k in range(len(sumsy)):
        totaly += (k + target_y - radius) * sumsy[k]
    c_y = totaly / np.sum(sumsy)

    sum_sd_y = 0
    for i in range(len(small_array)):
        sum_sd_y += ((i + target_y - radius - c_y) / np.sum(sumsy)) ** 2 * sumsy[i]
    sdy = math.sqrt(sum_sd_y)

    return c_x, c_y, sdx, sdy  # <- replace with x centroid, y centroid, uncertainty in x, uncertainty in y

centroid_x, centroid_y, uncert_x, uncert_y = findCentroid("sampleimage.fits", 351, 154, 1)
if abs(centroid_x - 350.9958) < 1e-3 and abs(centroid_y - 153.9955) < 1e-3:
    print("centroid calculation CORRECT")
else:
    print(f"centroid calculation INCORRECT, expected (350.9958, 153.9955), got ({centroid_x}, {centroid_y})")
if abs(uncert_x - 0.005254018) < 1e-6 and abs(uncert_y - 0.005249733) < 1e-6:
    print("uncertainty calculation CORRECT")
else:
    print(f"uncertainty calculation INCORRECT, expected (0.005254018, 0.005249733), got ({uncert_x}, {uncert_y})")

print(findCentroid("sampleimage.fits", 351, 154, 1))