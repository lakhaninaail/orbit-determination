import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.io import fits

#constants
dark = 10 #e-/pixel
gain = 0.8 #e-/ADU
read = 11 #e-/pixel

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

    for y in range(len(small_array)):
        sumsy.append(np.sum(small_array[y,:]))
    totaly = 0
    for k in range(len(sumsy)):
        totaly += (k + target_y - radius) * sumsy[k]
    c_y = totaly / np.sum(sumsy)

    return c_x, c_y  # <- replace with x centroid, y centroid, uncertainty in x, uncertainty in y

def include_pix_in_circle(r, fits_file, target_x, target_y, radius):
    centroid_x, centroid_y = findCentroid(fits_file, target_x, target_y, radius)
    c_x = int(centroid_x)
    c_y = int(centroid_y)
    arr_x = []
    arr_y = []
    x_range_low = c_x - r
    x_range_high = c_x + r + 1
    y_range_low = c_y - r
    y_range_high = c_y + r + 1
    for x in range(x_range_low, x_range_high):
        for y in range(y_range_low, y_range_high):
            if (centroid_x - x) ** 2 + (centroid_y - y) ** 2 <= r ** 2:
                arr_x.append(x)
                arr_y.append(y)
    #print("x: ", len(arr_x))
    #print("y: ", len(arr_y))
    coords = []
    for a in range(len(arr_x)):
            elem = [arr_x[a], arr_y[a]]
            coords.append(elem)
    coords = np.array(coords)
    return coords

def include_pix_between_circles(rin, rout, fits_file, target_x, target_y, radius):
    centroid_x, centroid_y = findCentroid(fits_file, target_x, target_y, radius)
    c_x = int(centroid_x)
    c_y = int(centroid_y)
    arr_x = []
    arr_y = []
    x_range_low = c_x - rout
    x_range_high = c_x + rout + 1
    y_range_low = c_y - rout
    y_range_high = c_y + rout + 1
    for x in range(x_range_low, x_range_high):
        for y in range(y_range_low, y_range_high):
            if (centroid_x - x) ** 2 + (centroid_y - y) ** 2 <= rout ** 2 and (centroid_x - x) ** 2 + (centroid_y - y) ** 2 >= rin ** 2:
                arr_x.append(x)
                arr_y.append(y)
    #print("x: ", len(arr_x))
    #print("y: ", len(arr_y))
    coords = []
    for a in range(len(arr_x)):
            elem = [arr_x[a], arr_y[a]]
            coords.append(elem)
    coords = np.array(coords)
    return coords
    
def find_apperture(fits_file, target_x, target_y, radius):
    coords_in_ap = include_pix_in_circle(radius, fits_file, target_x, target_y, radius)
    #print("SIZE: ", len(coords_in_ap))
    return coords_in_ap

def find_annulus(fits_file, target_x, target_y, radius, inner_radius, outer_radius):
    an_coords = include_pix_between_circles(inner_radius, outer_radius, fits_file, target_x, target_y, radius)
    #print("LENGTH ANN: ", len(an_coords))
    return an_coords

def adu_signal(fits_file, target_x, target_y, radius, inner_radius, outer_radius):
    counts = fits.getdata(fits_file)
    ap_pix = find_apperture(fits_file, target_x, target_y, radius)
    an_pix = find_annulus(fits_file, target_x, target_y, radius, inner_radius, outer_radius)
    ap_adu = 0
    an_adu = 0
    n_ap = 0
    n_an = 0
    for a in ap_pix:
        n_ap += 1
        x, y = a[0], a[1]
        ap_adu += counts[x,y]
    for b in an_pix:
        n_an += 1
        x, y = b[0], b[1]
        an_adu += counts[x,y]
    #print("n_ap = ",n_ap)
    #print("n_an = ",n_an)
    #print("AP ADU: ", ap_adu)
    #print("AN ADU: ", an_adu)
    #print("N_AP: ",len(ap_pix))
    #print("ADU SIGNAL: ",ap_adu - an_adu * n_ap / n_an)
    return ap_adu - an_adu * n_ap / n_an

def snr(fits_file, target_x, target_y, radius, inner_radius, outer_radius):
    counts = fits.getdata(fits_file)
    adu = adu_signal(fits_file, target_x, target_y, radius, inner_radius, outer_radius)
    signal = gain * adu
    ap = find_apperture(fits_file, target_x, target_y, radius)
    an = find_annulus(fits_file, target_x, target_y, radius, inner_radius, outer_radius)
    n_ap, n_an = len(ap), len(an)
    an_adu = 0
    for b in an:
        x, y = b[0], b[1]
        an_adu += counts[x,y]
    sky_e = an_adu / n_an
    electronic_noise_squared = read ** 2 + gain ** 2 / 12
    snr = (signal ** 0.5) / (1 + n_ap * (1 + n_ap / n_an) * ((sky_e + dark + electronic_noise_squared) / (signal))) ** 0.5
    return snr


def instMag(fits_file, target_x, target_y, radius, inner_radius, outer_radius):
    adu = adu_signal(fits_file, target_x, target_y, radius, inner_radius, outer_radius)
    inst_mag = -2.5 * math.log(adu, 10)
    return inst_mag

def instMag_uncertainty(fits_file, target_x, target_y, radius, inner_radius, outer_radius):
    signalnoiseratio = snr(fits_file, target_x, target_y, radius, inner_radius, outer_radius)
    return 1.0875 / signalnoiseratio

file = input("Enter File Name of Input Data: ")
x_approx = int(input("Enter Approximate x-Coordinate of Centroid: "))
y_approx = int(input("Enter Approximate y-Coordinate of Centroid: "))
ap_r = int(input("Enter Apperture Radius: "))
an_r_low = int(input("Enter Annulus Inner Radius: "))
an_r_high = int(input("Enter Annulus Outer Radius: "))

# centroid = findCentroid("aptest.FIT", 490, 293, 5)
# #print(include_pix_in_circle(5, "aptest.FIT", 490, 293, 5))
# #print(find_apperture("aptest.FIT", 490, 293, 5))
# #find_annulus("aptest.FIT", 490, 293, 5, 8, 13)
signal = adu_signal(file, x_approx, y_approx, ap_r, an_r_low, an_r_high)
inst_mag = instMag(file, x_approx, y_approx, ap_r, an_r_low, an_r_high)
uncertainty = instMag_uncertainty(file, x_approx, y_approx, ap_r, an_r_low, an_r_high)
snr = snr(file, x_approx, y_approx, ap_r, an_r_low, an_r_high)

print("Test Input File: {}".format(file))
print("Object Near (x,y) = ({},{})".format(x_approx, y_approx))
print("Aperture Radius = {}".format(ap_r))
print("Annulus Inner Radius = {}".format(an_r_low))
print("Annulus Outer Radius = {}".format(an_r_high))
print("Signal = {}".format(signal))
print("SNR = {}".format(snr))
print("Instrumental Magnitude = {} +/- {}".format(inst_mag, uncertainty))