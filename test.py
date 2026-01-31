
import pandas as pd
import numpy as np

from Python.parameters import * 

hcpb_u_rr_df   = pd.read_csv('./Figures/Data/HCPB_U238_rxns_900K.csv')
hcpb_th_rr_df  = pd.read_csv('./Figures/Data/HCPB_Th232_rxns_900K.csv')

df_u    = pd.DataFrame(HCPB_CONV_U_TBR, columns=['fertile_kgm3', 'ratio'])
df_th   = pd.DataFrame(HCPB_CONV_TH_TBR, columns=['fertile_kgm3', 'ratio'])

# NB. np.interp(x_values_to_calculate, original_x, original_y)
def get_u_ratio(x):
    return np.interp(x, df_u['fertile_kgm3'], df_u['ratio'])

def get_th_ratio(x):
    return np.interp(x, df_th['fertile_kgm3'], df_th['ratio'])

# Apply interpolation to the Uranium dataframe
hcpb_u_rr_df_corr = hcpb_u_rr_df.copy()
hcpb_u_rr_df_corr['tbr'] = hcpb_u_rr_df_corr['tbr'] * get_u_ratio(hcpb_u_rr_df_corr['fertile_kg/m3'])

# Apply interpolation to the Thorium dataframe
hcpb_th_rr_df_corr = hcpb_th_rr_df.copy()
hcpb_th_rr_df_corr['tbr'] = hcpb_th_rr_df_corr['tbr'] * get_th_ratio(hcpb_th_rr_df_corr['fertile_kg/m3'])

print(hcpb_u_rr_df['tbr'])
print(hcpb_u_rr_df_corr['tbr'])

print(get_u_ratio(hcpb_u_rr_df_corr['fertile_kg/m3']))