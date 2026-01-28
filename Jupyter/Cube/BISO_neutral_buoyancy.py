import numpy as np


density = {'UO2': 10.5,
           'FLiBe': 1.94, 'Pb-Li': 9.4, 
           'SiC': 3.2, 'PyC': 1.9, 'C': 1.05} # g/cc

radius  = {'UO2': 0.04} # cm 

c = 4/3*np.pi

def cube_root(x):
    return x**(1/3)


def main():
    
    for B in ['FLiBe', 'Pb-Li']:

        for C in ['PyC', 'SiC']:

            y = (density['UO2'] - density[C]) * c * radius['UO2']**3 / (density[B]*c - density[C]*c)

            r = cube_root(y)

            print(f"Breeder: {B} | Coating: {C} | Particle radius: {r:.4f} cm | diameter {2*r:.4f} cm")


if __name__ == '__main__':
    main()





