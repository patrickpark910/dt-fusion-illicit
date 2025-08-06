import math
import numpy as np
import pandas as pd
from scipy.interpolate import Akima1DInterpolator
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.font_manager import findSystemFonts
from matplotlib.font_manager import FontProperties


def reverse_and_append(x):
    x = list(x)
    reversed_part = x[-2::-1]  # Slicing to reverse excluding the last element
    return x + reversed_part

def project_forward(original_list, step_pattern, projections=2):
    extended_list = original_list.copy()  # To avoid modifying the original list
    pattern_length = len(step_pattern)
    
    for i in range(projections):
        step = step_pattern[i % pattern_length]
        next_element = extended_list[-1] + step
        list(extended_list).append(next_element)
    
    return extended_list



ZR_paths = ["..//MCNP//A-cladZr-Li000-bu10//elwr-coreA-cladZr-Li000-bu10.csv",
            "..//MCNP//A-cladZr-Li091-bu10//elwr-coreA-cladZr-Li091-bu10.csv",
            "..//MCNP//A-cladZr-Li139-bu10//elwr-coreA-cladZr-Li139-bu10.csv",
            "ZR"]

SS_paths = ["..//MCNP//A-cladSS-Li000-bu10//elwr-coreA-cladSS-Li000-bu10.csv",
            "..//MCNP//A-cladSS-Li057-bu10//elwr-coreA-cladSS-Li057-bu10.csv",
            "..//MCNP//A-cladSS-Li079-bu10//elwr-coreA-cladSS-Li079-bu10.csv",
            "SS"]

filepaths = [ZR_paths, SS_paths]

df1, df2, df3 = pd.read_csv(ZR_paths[0]), pd.read_csv(ZR_paths[1]), pd.read_csv(ZR_paths[2])
df4, df5, df6 = pd.read_csv(SS_paths[0]), pd.read_csv(SS_paths[1]), pd.read_csv(SS_paths[2])

b1 = np.array([x for x in df1['Burnup'] if not math.isnan(x)])
k1 = np.array([x for x in df1['k-eff'] if not math.isnan(x)])
d1 = np.array([x for x in df1['EFPD'] if not math.isnan(x)])

b2 = np.array(df2['Burnup'])
k2 = np.array(df2['k-eff'])
d2 = np.array(df2['EFPD'])

b3 = np.array(df3['Burnup'])
k3 = np.array(df3['k-eff'])
d3 = np.array(df3['EFPD'])

b4 = np.array([x for x in df4['Burnup'] if not math.isnan(x)])
k4 = np.array([x for x in df4['k-eff'] if not math.isnan(x)])
d5 = np.array([x for x in df4['EFPD'] if not math.isnan(x)])

b5 = np.array(df5['Burnup'])
k5 = np.array(df5['k-eff'])
d5 = np.array(df5['EFPD'])

b6 = np.array(df6['Burnup'])
k6 = np.array(df6['k-eff'])
d6 = np.array(df6['EFPD'])

b_fit = np.linspace(b1.min(), b1[:4].max(), 500)
coef11 = np.polyfit(b1[:4], k1[:4], deg=2)
quad11 = np.poly1d(coef11)(b_fit)
coef12 = np.polyfit(b1[3:], k1[3:], deg=2)
quad12 = np.poly1d(coef12)

# set up plot
font_path = 'C://Users//patri//Downloads//_Fonts//din//DIN-Regular.ttf'
fm.fontManager.addfont(font_path)
plt.rc('font',family='DIN',size=12) # ,stretch='normal')

plt.figure(figsize=(10, 4))

# zr clean
akima1 = Akima1DInterpolator(b1, k1) # Create Akima interpolator 
akima2 = Akima1DInterpolator(b2, k2)
akima3 = Akima1DInterpolator(b3, k3)
akima4 = Akima1DInterpolator(b4, k4)
akima5 = Akima1DInterpolator(b5, k5)
akima6 = Akima1DInterpolator(b6, k6)

x_fine = np.linspace(b1.min(), b1.max(), 300) # Evaluate on a fine grid
y_akima1 = akima1(x_fine)
y_akima2 = akima2(x_fine)
y_akima3 = akima3(x_fine)
y_akima4 = akima4(x_fine)
y_akima5 = akima5(x_fine)
y_akima6 = akima6(x_fine)

y_mono1 = np.copy(y_akima1) # Enforce monotonic decreasing post-processing
y_mono2 = np.copy(y_akima2)
y_mono3 = np.copy(y_akima3)
y_mono4 = np.copy(y_akima4)
y_mono5 = np.copy(y_akima5)
y_mono6 = np.copy(y_akima6)

for y_mono in [y_mono1,y_mono2,y_mono3,y_mono4,y_mono5,y_mono6]:
    for i in range(1, len(y_mono)):
        if y_mono[i] > y_mono[i-1]:
            # Adjust the current point to ensure it's not greater than the previous point
            y_mono[i] = y_mono[i-1]

# Plotting with specified colors and markers
plt.grid(True)
plt.axhline(y=1.03, color='#f8cc04', linewidth=2, linestyle='-')

plt.scatter(b1, k1, marker='.', color='black') # ZR clean
plt.plot(x_fine, y_akima1, linestyle='solid', color='black', linewidth=1.0)

plt.scatter(b2, k2, marker='+', color='#b41f24')  # ZR T(max)' 
plt.plot(x_fine, y_akima2, linestyle='solid',color='#b41f24', linewidth=1.0)

plt.scatter(b3, k3, marker='x', color='#0047ba')   # ZR T/WGPu
plt.plot(x_fine, y_akima3, linestyle='solid', color='#0047ba', linewidth=1.0)

plt.scatter(b4, k4, marker='.', color='black') # SS clean
plt.plot(x_fine, y_akima4, linestyle='dashed', color='black', linewidth=1.0)

plt.scatter(b5, k5, marker='+', color='#b41f24')  # SS T(max)
plt.plot(x_fine, y_akima5, linestyle='dashed', color='#b41f24', linewidth=1.0)

plt.scatter(b6, k6, marker='x', color='#0047ba')   # SS T/WGPu
plt.plot(x_fine, y_akima6, linestyle='dashed', color='#0047ba', linewidth=1.0)


# Dummy functions for legend
plt.plot([100,101], [100,101], label='SS-clad clean',  marker='.', linestyle='dashed', color='black', linewidth=1.5, alpha=1)
plt.plot([100,101], [100,101], label='SS-clad tritium', marker='+', linestyle='dashed', color='#b41f24', linewidth=1.5, alpha=1)
plt.plot([100,101], [100,101], label='SS-clad co-prod.', marker='x', linestyle='dashed', color='#0047ba', linewidth=1.5, alpha=1)
plt.plot([100,101], [100,101], label='ZR-clad clean', marker='.', linestyle='solid', color='black', linewidth=1.5, alpha=1)
plt.plot([100,101], [100,101], label='ZR-clad tritium', marker='+', linestyle='solid', color='#b41f24', linewidth=1.5, alpha=1)
plt.plot([100,101], [100,101], label='ZR-clad co-prod.', marker='x', linestyle='solid', color='#0047ba', linewidth=1.5, alpha=1)

 
# Setting font to Arial
# plt.title('Burnup vs. k-eff', fontname='Arial')
plt.xlabel('Burnup [MWd/kgU]')
plt.ylabel('k-eff')



# Setting axis limits
plt.ylim(0.99, 1.31)
plt.xlim(-0.5, 28.5)

plt.xticks(range(0, 29, 2))  # X-axis: intervals of 2 from 0 to 20
plt.yticks(np.arange(1.00, 1.31, 0.03))
plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(1))
plt.tick_params(which='minor', direction="in", length=3, width=1) #, color='gray')  # Same length as major ticks
plt.tick_params(which='major', direction="in", length=4, width=1)



# Adding legend, grid, and customizing appearance
plt.legend(fancybox=False, edgecolor='black', frameon=True, framealpha=.75, ncol=2)


ax_top = plt.twiny()
# top_ticks = [0,2*35.26,4*35.26,6*35.26,8*35.26,10*35.26,12*35.26,14*35.26,16*35.26,18*35.26,20*35.26,22*35.26,24*35.26,26*35.26,28*35.26]
top_ticks = [0,5.61*35.26,11.2*35.26,22.8*35.26,8.27*35.26,14.9*35.26]
ax_top.set_xlim(-18, 1005)
ax_top.set_xticks(top_ticks)
formatter = ticker.FuncFormatter(lambda x, pos: f"{round(x):.0f}")
ax_top.xaxis.set_major_formatter(formatter)
ax_top.set_xlabel("Effective Full-Power Days (EFPD)")
ax_top.tick_params(which='major', direction="in", length=4, width=1)


plt.savefig(f'fig4-keff-smooth-efpd.svg', format='svg', bbox_inches='tight')
plt.show()



