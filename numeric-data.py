# integers are type int
# floating point numbers (real numbers) are type float
# use int and float functions to convert from text to numbers

# represent 1/3 -> 0.3, 0.33, 0.333, ...
print((2.5*0.1)*1.5 == 2.5*(0.1*1.5))
print(0.1*1.5*2.5 == 2.5*0.1*1.5)
print(0.1 + 0.1 + 0.1 == 0.3)

print(format(0.1, ".30f"))
print(format(0.3, ".30f"))
print(format(0.1 + 0.1 + 0.1, ".30f"))
# be very careful when comparing floats for equality
# better to compare an absolute difference to a threshold
# instead of a == b, abs(a - b) < 1e-9