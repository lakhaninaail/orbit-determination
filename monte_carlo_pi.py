import math
import numpy as np

def pi_approx():
    counter = 0
    inside = 0
    pi_sum = 0
    while counter < 10000:
        x = 2 * np.random.random()
        y = 2 * np.random.random()
        if (x - 1) ** 2 + (y - 1) ** 2 <= 1:
            inside += 1
        counter += 1

    return 4 * inside / counter

def errors():
    pi_sum = 0
    sd_sum = 0
    n = 1000
    for k in range(n):
        approx = pi_approx()
        pi_sum += approx
        sd_sum += (math.pi - approx) ** 2
    mean = pi_sum / n
    sd = (sd_sum / (n-1)) ** 0.5
    sdom = sd / (n ** 0.5)
    percent_error = 100 * (mean - math.pi) / math.pi
    return percent_error, mean, sd, sdom

print(errors())