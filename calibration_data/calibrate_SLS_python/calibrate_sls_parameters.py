import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import griddata



# Import the data
raw_data = []
with open('exported_visco_properties_ALL_data.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0

    for row in csv_reader:

        # Skip first two header lines
        if ((line_count == 0) or (line_count == 1)):
            line_count += 1
            continue

        else:
            row_floats = [float(x) for x in row]
            raw_data.append(row_floats)
            line_count += 1

raw_data = np.array(raw_data)

# Organize the data
data_0da = raw_data[:,0:5]
data_1da = raw_data[:,5:10]
data_2da = raw_data[:,10:15]
data_3da = raw_data[:,15:20]

# Delete empty data rows denoted by -1 values
data_0da_list = []
for ii, row in enumerate(data_0da):
    if row[0] > 0:
        data_0da_list.append(row)

data_1da_list = []
for ii, row in enumerate(data_1da):
    if row[0] > 0:
        data_1da_list.append(row)

data_2da_list = []
for ii, row in enumerate(data_2da):
    if row[0] > 0:
        data_2da_list.append(row)

data_3da_list = []
for ii, row in enumerate(data_3da):
    if row[0] > 0:
        data_3da_list.append(row)

data_0da = np.array(data_0da_list)
data_1da = np.array(data_1da_list)
data_2da = np.array(data_2da_list)
data_3da = np.array(data_3da_list)

data_012da = np.vstack((data_0da,data_1da,data_2da))
data_0123da = np.vstack((data_0da,data_1da,data_2da, data_3da))


"""
# Interpolate the data onto a regular grid of temperatures and dosages
# for plotting purposes later

# Temperature range in Kelvin
xlin = np.linspace(259.0, 333.0, num=60, endpoint=True)

# Dosage range in Grays. 
ylin = np.linspace(0.0, 16200, num=10, endpoint=True)

# Make a regular 2D grid
xgrid, ygrid = np.meshgrid(xlin, ylin)

zgrid_gamma = griddata(data_012da[:,0:2], data_012da[:,3], (xgrid, ygrid), 
    method='cubic')

zgrid_tau = griddata(data_012da[:,0:2], data_012da[:,4], (xgrid, ygrid), 
    method='cubic')

# Plot
fig1, ax1 = plt.subplots()
c_obj = ax1.contourf(xgrid, ygrid, zgrid_gamma)
cbar = fig1.colorbar(c_obj)
plt.show()
"""


# Definitions of the polynomial surfaces to be fitted to the 1D and 2D data
def func_poly41(X, a40, a31, a30, a21, a20, a11, a10, a01, a00):
    x,y = X # Separate the tuple for clarity

    f4 = a40*(x**4) + a31*(x**3)*y
    f3 = a30*(x**3) + a21*(x**2)*y
    f2 = a20*(x**2) + a11*x*y
    f1 = a10*x + a01*y
    f0 = a00

    f41 = f4 + f3 + f2 + f1 + f0

    return f41

# Definitions of the polynomial surfaces to be fitted to the 1D and 2D data
def func_poly51(X, a50, a41, a40, a31, a30, a21, a20, a11, a10, a01, a00):
    x,y = X # Separate the tuple for clarity

    f5 = a50*(x**5) + a41*(x**4)*y
    f4 = a40*(x**4) + a31*(x**3)*y
    f3 = a30*(x**3) + a21*(x**2)*y
    f2 = a20*(x**2) + a11*x*y
    f1 = a10*x + a01*y
    f0 = a00

    f51 = f5 + f4 + f3 + f2 + f1 + f0

    return f51

#!!!
def func_poly71(X, b10, b11, b12, b13, b20, b21, b22, b23, a70, a61, a60, a51, a50, a41, a40, a31, a30, a21, a20, a11, a10, a01, a00):
    x,y = X # Separate the tuple for clarity

    f10 = b20*y*np.tanh(b21*x) + b22*x*np.tanh(b23*y**2)
    f9 = b20*x*np.tanh(b21*x) + b22*y*np.tanh(b23*y**2)
    f8 = b10*np.tanh(b11*x) + b12*np.tanh(b13*y)
    f7 = a70*(x**7) + a61*(x**6)*y
    f6 = a60*(x**6) + a51*(x**5)*y
    f5 = a50*(x**5) + a41*(x**4)*y
    f4 = a40*(x**4) + a31*(x**3)*y
    f3 = a30*(x**3) + a21*(x**2)*y
    f2 = a20*(x**2) + a11*x*y
    f1 = a10*x + a01*y
    f0 = a00

    f71 = f10 + f9 + f8 + f7 + f6 + f5 + f4 + f3 + f2 + f1 + f0

    return f71

def func_poly42(X, a40, a31, a22, a30, a21, a12, a20, a11, a02, a10, a01, a00):
    x,y = X # Separate the tuple for clarity

    f4 = a40*(x**4) + a31*(x**3)*y + a22*(x**2)*(y**2)
    f3 = a30*(x**3) + a21*(x**2)*y + a12*x*(y**2)
    f2 = a20*(x**2) + a11*x*y + a02*(y**2)
    f1 = a10*x + a01*y
    f0 = a00

    f42 = f4 + f3 + f2 + f1 + f0

    return f42


# ---------------- FIT GAMMA ---------------- 

# Fit polynomial surface to the gamma data
xdata = data_0123da[:,0] # temperature
ydata = data_0123da[:,1] # dosage
zdata = data_0123da[:,3] # gamma

popt, pcov = curve_fit(func_poly42, (xdata, ydata),
    zdata)

print(f"\n\nFitted poly42 gamma parameters:")
print(f"a40 = {popt[0]}")
print(f"a31 = {popt[1]}")
print(f"a22 = {popt[2]}")
print(f"a30 = {popt[3]}")
print(f"a21 = {popt[4]}")
print(f"a12 = {popt[5]}")
print(f"a20 = {popt[6]}")
print(f"a11 = {popt[7]}")
print(f"a02 = {popt[8]}")
print(f"a10 = {popt[9]}")
print(f"a01 = {popt[10]}")
print(f"a00 = {popt[11]}")

model_z = func_poly42((xdata, ydata), *popt) 

absError = model_z - zdata
SE = np.square(absError) # squared errors
MSE = np.mean(SE) # mean squared errors
RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
Rsquared = 1.0 - (np.var(absError) / np.var(zdata))
print('\nRMSE:', RMSE)
print('R-squared:', Rsquared)


#!!! Fit polynomial surface to the gamma data
xdata = data_012da[:,0] # temperature
ydata = data_012da[:,1] # dosage
zdata = data_012da[:,3] # gamma

popt, pcov = curve_fit(func_poly71, (xdata, ydata),
    zdata)

print(f"\n\nFitted poly71 gamma parameters:")
print(f"a70 = {popt[0]}")
print(f"a61 = {popt[1]}")
print(f"a60 = {popt[2]}")
print(f"a51 = {popt[3]}")
print(f"a50 = {popt[4]}")
print(f"a41 = {popt[5]}")
print(f"a40 = {popt[6]}")
print(f"a31 = {popt[7]}")
print(f"a30 = {popt[8]}")
print(f"a21 = {popt[9]}")
print(f"a20 = {popt[10]}")
print(f"a11 = {popt[11]}")
print(f"a10 = {popt[12]}")
print(f"a01 = {popt[13]}")
print(f"a00 = {popt[14]}")

model_z = func_poly71((xdata, ydata), *popt) 

absError = model_z - zdata
SE = np.square(absError) # squared errors
MSE = np.mean(SE) # mean squared errors
RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
Rsquared = 1.0 - (np.var(absError) / np.var(zdata))
print('\nRMSE:', RMSE)
print('R-squared:', Rsquared)


# ---------------- FIT TAU ---------------- 

# Fit polynomial surface to the tau (relaxation time) data
xdata = data_0123da[:,0] # temperature
ydata = data_0123da[:,1] # dosage
zdata = data_0123da[:,4] # tau

popt, pcov = curve_fit(func_poly42, (xdata, ydata),
    zdata)

print(f"\n\nFitted poly42 relaxation time (tau) parameters:")
print(f"a40 = {popt[0]}")
print(f"a31 = {popt[1]}")
print(f"a22 = {popt[2]}")
print(f"a30 = {popt[3]}")
print(f"a21 = {popt[4]}")
print(f"a12 = {popt[5]}")
print(f"a20 = {popt[6]}")
print(f"a11 = {popt[7]}")
print(f"a02 = {popt[8]}")
print(f"a10 = {popt[9]}")
print(f"a01 = {popt[10]}")
print(f"a00 = {popt[11]}")

model_z = func_poly42((xdata, ydata), *popt) 

absError = model_z - zdata
SE = np.square(absError) # squared errors
MSE = np.mean(SE) # mean squared errors
RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
Rsquared = 1.0 - (np.var(absError) / np.var(zdata))
print('\nRMSE:', RMSE)
print('R-squared:', Rsquared)


# Fit polynomial surface to the tau (relaxation time) data
xdata = data_012da[:,0] # temperature
ydata = data_012da[:,1] # dosage
zdata = data_012da[:,4] # tau

popt, pcov = curve_fit(func_poly51, (xdata, ydata),
    zdata)

print(f"\n\nFitted poly41 relaxation time (tau) parameters:")
print(f"a50 = {popt[0]}")
print(f"a41 = {popt[1]}")
print(f"a40 = {popt[2]}")
print(f"a31 = {popt[3]}")
print(f"a30 = {popt[4]}")
print(f"a21 = {popt[5]}")
print(f"a20 = {popt[6]}")
print(f"a11 = {popt[7]}")
print(f"a10 = {popt[8]}")
print(f"a01 = {popt[9]}")
print(f"a00 = {popt[10]}")

model_z = func_poly51((xdata, ydata), *popt) 

absError = model_z - zdata
SE = np.square(absError) # squared errors
MSE = np.mean(SE) # mean squared errors
RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
Rsquared = 1.0 - (np.var(absError) / np.var(zdata))
print('\nRMSE:', RMSE)
print('R-squared:', Rsquared)


# ---------------- FIT YOUNG'S MODULUS ---------------- 

# Fit polynomial line to the base spring (Young's modulus) data. Note, this
# value only changes with dosage -- it is not temperature dependent.
xdata = data_0123da[:,1] # dosage
zdata = data_0123da[:,2] # young's modulus

# Fit a polynomial p(x) = p[0] * x**deg + ... + p[deg]
# For deg = 2, p(x) = p[0] * x**2 + p[1] * x + p[2]
popt = np.polyfit(xdata, zdata, deg=2)

print(f"\n\nFitted poly2 Young's modulus parameters:")
print(f"a2 * x**2 + a1 * x + a0")
print(f"a2 = {popt[0]}")
print(f"a1 = {popt[1]}")
print(f"a0 = {popt[2]}")

model_z = np.polyval(popt, xdata) 

absError = model_z - zdata
SE = np.square(absError) # squared errors
MSE = np.mean(SE) # mean squared errors
RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
Rsquared = 1.0 - (np.var(absError) / np.var(zdata))
print('\nRMSE:', RMSE)
print('R-squared:', Rsquared)


# Fit polynomial line to the base spring (Young's modulus) data. Note, this
# value only changes with dosage -- it is not temperature dependent.
xdata = data_012da[:,1] # dosage
zdata = data_012da[:,2] # young's modulus

# Fit a polynomial p(x) = p[0] * x**deg + ... + p[deg]
# For deg = 1, p(x) = p[0] * x + p[1]
popt = np.polyfit(xdata, zdata, deg=1)

print(f"\n\nFitted poly1 Young's modulus parameters:")
print(f"a1 * x + a0")
print(f"a1 = {popt[0]}")
print(f"a0 = {popt[1]}")

model_z = np.polyval(popt, xdata) 

absError = model_z - zdata
SE = np.square(absError) # squared errors
MSE = np.mean(SE) # mean squared errors
RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
Rsquared = 1.0 - (np.var(absError) / np.var(zdata))
print('\nRMSE:', RMSE)
print('R-squared:', Rsquared)