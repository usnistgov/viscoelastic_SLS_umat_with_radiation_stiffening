import csv
import glob
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

# Find filenames to be imported from "dat2txt.py"
file_name_strain = glob.glob("strains elem*")
file_name_stress = glob.glob("stresses elem*")

print(f"\nImporting '{file_name_strain[0]}'")
df_strain = pd.read_csv(file_name_strain[0], sep='\s+', header=None,
    skip_blank_lines=True, dtype=np.float64)
strain_arr = df_strain.to_numpy()

print(f"Importing '{file_name_stress[0]}'")
df_stress = pd.read_csv(file_name_stress[0], sep='\s+', header=None,
    skip_blank_lines=True, dtype=np.float64)
stress_arr = df_stress.to_numpy()

del df_strain, df_stress

num_time_pnts = stress_arr.shape[0]
num_val_pnts = stress_arr.shape[1]

if num_val_pnts > 9:
    print("\nAssuming 2x2x2 Gaussian integration was used. Averaging integration points...")
    strain_arr_temp = np.zeros((num_time_pnts,9), dtype=float)
    stress_arr_temp = np.zeros((num_time_pnts,9), dtype=float)
    
    for m, cur_row in enumerate(strain_arr):
        time_i = cur_row[0]

        el_num1_i = cur_row[1]
        pnt1_i = cur_row[2]
        vals1_i = cur_row[3:9]

        el_num2_i = cur_row[9]
        pnt2_i = cur_row[10]
        vals2_i = cur_row[11:17]

        el_num3_i = cur_row[17]
        pnt3_i = cur_row[18]
        vals3_i = cur_row[19:25]

        el_num4_i = cur_row[25]
        pnt4_i = cur_row[26]
        vals4_i = cur_row[27:33]

        el_num5_i = cur_row[33]
        pnt5_i = cur_row[34]
        vals5_i = cur_row[35:41]

        el_num6_i = cur_row[41]
        pnt6_i = cur_row[42]
        vals6_i = cur_row[43:49]

        el_num7_i = cur_row[49]
        pnt7_i = cur_row[50]
        vals7_i = cur_row[51:57]

        el_num8_i = cur_row[57]
        pnt8_i = cur_row[58]
        vals8_i = cur_row[59:65]

        vals_sum = vals1_i + vals2_i + vals3_i + vals4_i + vals5_i + \
            vals6_i + vals7_i + vals8_i

        vals_avg = vals_sum/8

        strain_arr_temp[m,0] = time_i
        strain_arr_temp[m,1] = el_num1_i
        strain_arr_temp[m,2] = 1
        strain_arr_temp[m,3:9] = vals_avg

    for m, cur_row in enumerate(stress_arr):
        time_i = cur_row[0]

        el_num1_i = cur_row[1]
        pnt1_i = cur_row[2]
        vals1_i = cur_row[3:9]

        el_num2_i = cur_row[9]
        pnt2_i = cur_row[10]
        vals2_i = cur_row[11:17]

        el_num3_i = cur_row[17]
        pnt3_i = cur_row[18]
        vals3_i = cur_row[19:25]

        el_num4_i = cur_row[25]
        pnt4_i = cur_row[26]
        vals4_i = cur_row[27:33]

        el_num5_i = cur_row[33]
        pnt5_i = cur_row[34]
        vals5_i = cur_row[35:41]

        el_num6_i = cur_row[41]
        pnt6_i = cur_row[42]
        vals6_i = cur_row[43:49]

        el_num7_i = cur_row[49]
        pnt7_i = cur_row[50]
        vals7_i = cur_row[51:57]

        el_num8_i = cur_row[57]
        pnt8_i = cur_row[58]
        vals8_i = cur_row[59:65]

        vals_sum = vals1_i + vals2_i + vals3_i + vals4_i + vals5_i + \
            vals6_i + vals7_i + vals8_i

        vals_avg = vals_sum/8

        stress_arr_temp[m,0] = time_i
        stress_arr_temp[m,1] = el_num1_i
        stress_arr_temp[m,2] = 1
        stress_arr_temp[m,3:9] = vals_avg

    strain_arr = strain_arr_temp.copy()
    stress_arr = stress_arr_temp.copy()

    del strain_arr_temp, stress_arr_temp

# Now strains exist in strain_arr with columns of:
# [time, elem ID, integ pnt ID, Exx, Eyy, Ezz, Exy, Exz, Eyz]
#
# And stresses exist in stress_arr with columns of:
# [time, elem ID, integ pnt ID, Sxx, Syy, Szz, Sxy, Sxz, Syz]

# Assuming that loadings occur in the Y-direction
time_arr = strain_arr[:,0]
eps_arr = strain_arr[:,4]
sig_arr = stress_arr[:,4]

def func_sin(time, amplitude, frequency, phase):
    return amplitude*(np.sin(frequency*time + phase))

# Use only the latter portion of the data for fitting, after
# everything has settled out.
t_fit_start = int(np.round(num_time_pnts/2))

time_arr_fit = time_arr[t_fit_start:]
eps_arr_fit = eps_arr[t_fit_start:]
sig_arr_fit = sig_arr[t_fit_start:]

print(f"\nFitting sin() function from: {time_arr_fit[0]} <= time <= {time_arr_fit[-1]}")

# Curve fit for strain data
eps_init = [np.amax(eps_arr)/2.0, 10.0, 0.0]
eps_opt, eps_cov = curve_fit(func_sin, time_arr_fit, eps_arr_fit, p0=eps_init)

# Curve fit for strain data
sig_init = [np.amax(sig_arr)/2.0, 10.0, 0.0]
sig_opt, sig_cov = curve_fit(func_sin, time_arr_fit, sig_arr_fit, p0=sig_init)

print(f"\nFitted parameters for strain data:")
print(f"  strain(t) = a*sin(w*t + d)")
print(f"  a = {eps_opt[0]}")
print(f"  w = {eps_opt[1]}")
print(f"  d = {eps_opt[2]}  <-- Should be close to zero")

print(f"\nFitted parameters for stress data:")
print(f"  stress(t) = a*sin(w*t + d)")
print(f"  a = {sig_opt[0]}")
print(f"  w = {sig_opt[1]}")
print(f"  d = {sig_opt[2]}")

phase_angle = np.absolute(sig_opt[2] - eps_opt[2]) # [radians]
loss_tangent = np.tan(phase_angle) 
storage_mod = (sig_opt[0]/eps_opt[0])*np.cos(phase_angle)
loss_mod = (sig_opt[0]/eps_opt[0])*np.sin(phase_angle)

print(f"\nAssuming stresses were provided in [Pa]. Summary of results:")
print(f"  Loss tangent, tan(delta) = {loss_tangent}")
print(f"  Storage modulus (E') = {storage_mod} [Pa]")
print(f"  Loss modulus (E'') = {loss_mod} [Pa]")

eps_model = np.zeros(time_arr.shape, dtype=float)
sig_model = np.zeros(time_arr.shape, dtype=float)

for m, time in enumerate(time_arr):
    eps_model[m] = func_sin(time, *eps_opt)
    sig_model[m] = func_sin(time, *sig_opt)

visc_fitted_path_out = "./viscoelastic_fitted_data.csv"

print(f"\nSaving fitted data and cofficents: {visc_fitted_path_out}")

with open(visc_fitted_path_out, 'w') as f_obj:

    # create the csv writer
    writer = csv.writer(f_obj)

    head_str = ["time (s)", "strain_yy", "stress_yy (Pa)", 
        "fitted strain", "fitted_stress", "loss tangent, tan(delta)",
        "storage modulus, E' (Pa)", "loss modulus, E'' (Pa)", 
        "fitted strain amplitude", "fitted strain freq (rad/s)", 
        "fitted strain phase angle (rad)", "fitted stress amplitude (Pa)", 
        "fitted stress freq (rad/s)", "fitted stress phase angle (rad)"]

    # write a row to the csv file
    writer.writerow(head_str)

    for m, cur_t in enumerate(time_arr):

        temp_row = [cur_t, eps_arr[m], sig_arr[m], eps_model[m], sig_model[m],
            loss_tangent, storage_mod, loss_mod, eps_opt[0], eps_opt[1], 
            eps_opt[2], sig_opt[0], sig_opt[1], sig_opt[2]]

        temp_row_str = [str(x) for x in temp_row]

        writer.writerow(temp_row_str)

print(f"\nPost-processing of virtual DMA experiment complete!")