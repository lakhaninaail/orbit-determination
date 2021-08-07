from astropy.io import fits
table = fits.open("corr5.fits")[1].data
print(table)