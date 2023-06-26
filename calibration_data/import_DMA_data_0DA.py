import os
import math
import statistics as stats
import pandas as pd
import numpy as np
import requests as rq
import h5py as hdf

# Directory where the .xls files are located
directory = './VN01_001_DMA_results/0DA'
output_file = './VN01_001_DMA_results/0DA/visco_properties_0DA_data.xlsx'

strain_total = 0.0924292 + 0.0000827615  #[mm/mm]
stress_t1600 = 22707.3 #Pa

# Define the mat. props. to use
E_0 = (stress_t1600/strain_total)/10e6 #MPa, for 0 day dose (stress relaxation experiment @ 1600 s)
freq = 10.0 #rad/s, test setup param

# Load data and compute poisson's ratio
data_location = 'https://data.materialsdatafacility.org/mdf_open/foam_db_v1.1/quasistatic_rate_data/VN01/VN01_001_003_QS06_00/VN01_001_003_QS06_00_dic.mat'
dest_dic = "./VN01_001_003_QS06_00_dic.mat"
r_dic_data = rq.get(data_location)
open(dest_dic , 'wb+').write(r_dic_data.content)
dic_data = hdf.File(dest_dic)
strain_e11 = (np.array(dic_data['complete_data']['E'][1, :]) - 1)
strain_e22_c = 1.0*(np.array(dic_data['complete_data']['E'][0, :]) - 1)
poisson_ratio = -stats.mean(strain_e11[1:20]/strain_e22_c[1:20])


# List to store the data from each file
data = []

# Iterate over each file in the directory
for filename in os.listdir(directory):
    if filename.endswith('.xls'):
        # Full path of the file
        file_path = os.path.join(directory, filename)
        
        # Read the Excel file into a pandas DataFrame
        df = pd.read_excel(file_path,sheet_name="Time - 1")
        # Append the DataFrame to the data list
        df.drop(index=df.index[0:2], axis=0, inplace=True)
        data.append(df)

        # Now try to do more sheets
        try: 
            df = pd.read_excel(file_path,sheet_name="Time - 2")
            df.drop(index=df.index[0:2], axis=0, inplace=True)
            data.append(df)
            df = pd.read_excel(file_path,sheet_name="Time - 3")
            df.drop(index=df.index[0:2], axis=0, inplace=True)
            data.append(df)
            df = pd.read_excel(file_path,sheet_name="Time - 4")
            df.drop(index=df.index[0:2], axis=0, inplace=True)
            data.append(df)
            df = pd.read_excel(file_path,sheet_name="Time - 5")
            df.drop(index=df.index[0:2], axis=0, inplace=True)
            data.append(df)
            df = pd.read_excel(file_path,sheet_name="Time - 6")
            df.drop(index=df.index[0:2], axis=0, inplace=True)
            data.append(df)
            df = pd.read_excel(file_path,sheet_name="Time - 7")
            df.drop(index=df.index[0:2], axis=0, inplace=True)
            data.append(df)
            df = pd.read_excel(file_path,sheet_name="Time - 8")
            df.drop(index=df.index[0:2], axis=0, inplace=True)
            data.append(df)
        except ValueError:
            pass


# Properties to compute
tau = []
gamma = []
temps = []

# Loop through all temperatures
for curr_temp in data:
    #compute the current relax time
    tau_i = stats.mean(1.0/(freq*curr_temp['Unnamed: 5']))
    tau.append(tau_i)

    #compute the current normalized stiffness 
    Ei = stats.mean(curr_temp['Unnamed: 6'])*(1.0 + freq**2*tau_i**2)/(freq**2*tau_i**2)
    gamma.append(Ei/E_0)
    
    #save the associated temp too
    temps.append(stats.mean(curr_temp['Unnamed: 2']))

# Combine into a dataframe
complete_data = pd.DataFrame({'Base spring modulus [MPa]':E_0,'Poisson ratio':poisson_ratio,'Temperature [degC]':temps,'Tau':tau,'Gamma':gamma})
#Write out the excel
with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
    complete_data.to_excel(writer, sheet_name='data', header=True)

print('Data written to ',output_file)
