from odlib import *

a = semimajorAxis("LakhaniInput.txt", "00:00:00.0000")
e = eccentricity("LakhaniInput.txt", "00:00:00.0000")
i = inclination("LakhaniInput.txt", "00:00:00.0000")
loan = LoAN("LakhaniInput.txt", "00:00:00.0000")
aop = AoP("LakhaniInput.txt", "00:00:00.0000")
M = meanAnomaly("LakhaniInput.txt", "00:00:00.0000")
lpp = LPP("LakhaniInput.txt", "00:00:00.0000")

# all expected values and percent differences are for 00:00:00.0000
a_ = 1.056800392709198E+00
e_ = 3.442329555807205E-01
i_ = 2.515525486975122E+01
loan_ = 2.362380026346535E+02
aop_ = 2.555046578188225E+02
M_ = 1.404194630208429E+02
lppExpected = getExpectedLPP("LakhaniInput.txt", "00:00:00.0000")

errorLPP = errorLPP("LakhaniInput.txt", "00:00:00.0000")

print("Julian Date: {}".format(getJulian("LakhaniInput.txt", "00:00:00.0000")))
print("Semimajor Axis: Expected Value = {} Calculated Value = {} Percent Error = {}".format(a_, a, 100 * (a_ - a) / a_))
print("Eccentricity: Expected Value = {} Calculated Value = {} Percent Error = {}".format(e_, e, 100 * (e_ - e) / e_))
print("Inclination: Expected Value = {} Calculated Value = {} Percent Error = {}".format(i_, i, 100 * (i_ - i) / i_))
print("Longitude of Ascending Node: Expected Value = {} Calculated Value = {} Percent Error = {}".format(loan_, loan, 100 * (loan_ - loan) / loan_))
print("Argument of Perihelion: Expected Value = {} Calculated Value = {} Percent Error = {}".format(aop_, aop, 100 * (aop_ - aop) / aop_))
print("Mean Anomaly: Expected Value = {} Calculated Value = {} Percent Error = {}".format(M_, M, 100 * (M_ - M) / M_))
print("Last Perihelion Passage: Expected Value = {} Calculated Value = {} Percent Error = {}".format(lppExpected, lpp, errorLPP))