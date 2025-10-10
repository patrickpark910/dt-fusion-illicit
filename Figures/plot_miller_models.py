import numpy as np
import matplotlib.pyplot as plt


"""
This file is meant to be standalone (execute this by itself, independent of main.py)
"""

class Specs:
    def __init__(self, name:str, R0:float, a:float, kappa:float, delta:float, layers:list):
        self.name   = name
        self.R0, self.a, self.kappa, self.delta = R0, a, kappa, delta
        self.aspect = self.R0 / self.a
        self.layers = layers

        self.surf_area = miller_surface_area(R0, a, kappa, delta)
        # print(f"{self.name} plasma-facing surface area: {self.surf_area/1e4:.4f} m^2")

    def __repr__(self):
        return (f"TokamakReactor(R0= {self.R0:.4f} m, a={self.a:.4f} m, "
                f"kappa= {self.kappa:.4f}, delta={self.delta:.4f})")


def main():

    ARC_Sorbom     = Specs('ARC (Sorbom 15)',     
                           330, 113, 1.84, 0.45, [(0.3,'fw'), (1,'st1'), (2,'br1'), (3,'st2'), (100,'br2'), (3,'st3')])
    ARC_Ball_code  = Specs('ARC-class (Ball 25 code)',
                           400, 120, 1.5 , 0.5 , [(0.3,'fw'), (1,'st1'), (2,'br1'), (3,'st2'), (100,'br2'), (3,'st3')])
    ARC_Ball_paper = Specs('ARC-class (Ball 25 paper)',
                           400, 100, 1.6 , 0.5 , [(0.3,'fw'), (1,'st1'), (2,'br1'), (3,'st2'), (100,'br2'), (3,'st3')])
    Emma_FLiBe     = Specs('Emma FLiBe',
                           600, 200, 1.72, 0.4 , [(0.3,'fw'), (1,'st1'), (2,'br1'), (3,'st2'), (100,'br2'), (3,'st3')])
    
    Emma_LL        = Specs('Emma LL/PB',
                           620, 207, 1.72, 0.4 , [(0.2,'fw'), (0.4,'st1'), (2,'br1'), (0.4,'st2'), (22.5,'br2'), (3.2,'st3'), (21.0,'br3'), (3.2,'st4'), (21.0,'br4'), (8.0,'st5')])
    
    EU_DEMO        = Specs('EU DEMO',
                           907.2, 292.7, 1.59, 0.33, [(0.2,'fw'), (0.4,'st1'), (2,'br1'), (0.4,'st2'), (22.5,'br2'), (3.2,'st3'), (21.0,'br3'), (3.2,'st4'), (21.0,'br4'), (8.0,'st5')]) 

    # plot_separate([ARC_Sorbom, ARC_Ball_code, Emma_FLiBe, Emma_LL, EU_DEMO]) # , ITER_LL, EU_DEMO_PB]
    plot_together([ARC_Ball_code]) # , ITER_LL, EU_DEMO_PB]


def plot_separate(reactors_to_plot, n=10000):

    fig, ax = plt.subplots(1, len(reactors_to_plot), figsize=(7*len(reactors_to_plot), 6))

    # Create parameter array
    i, t = 0, np.linspace(0, 2*np.pi, n)

    # Define colors for different reactors
    colors = plt.cm.Set1(np.arange(len(reactors_to_plot))) 

    for reactor in reactors_to_plot:
        color = colors[i]

        ax[i].set_xlim(100, 1300)
        ax[i].set_ylim(-600, 600)
        ax[i].set_xlabel('radius (cm)', fontsize=10)
        ax[i].set_ylabel('height (cm)', fontsize=10)
        ax[i].set_title(f'{reactor.name}', fontsize=10)
        ax[i].grid(True, alpha=0.3)
        # ax[i].axis('equal')
        ax[i].legend(loc='best', fontsize=8)

        # Plasma shape
        R, Z, V = miller_model(reactor.R0, reactor.a, reactor.kappa, reactor.delta, 0)
        ax[i].plot(R, Z, '-', color=color, linewidth=1) # , label='Original Miller D-shape')
        
        offset = 0
        for layer in reactor.layers:
            offset += layer[0]
            R_offset, Z_offset, V = miller_model(reactor.R0, reactor.a, reactor.kappa, reactor.delta, offset)
            ax[i].plot(R_offset, Z_offset, '-', color=color, linewidth=1, label=f'{layer[1]}')

        i+=1

    plt.tight_layout()
    # plt.show()
    # plt.savefig('tokamaks_separate_lowres.png',dpi=300,bbox_inches='tight',transparent=False)  



def plot_together(reactors_to_plot, n=10000):
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))  # Single subplot
    
    # Create parameter array
    t = np.linspace(0, 2*np.pi, n)
    
    # Define colors for different reactors
    colors = plt.cm.Set1(np.arange(len(reactors_to_plot))) 
    
    for i, reactor in enumerate(reactors_to_plot):
        color = colors[i]
        
        # Plasma shape
        R, Z, V = miller_model(reactor.R0, reactor.a, reactor.kappa, reactor.delta, 0)
        ax.plot(R, Z, color=color, linewidth=1, label=f'{reactor.name}')
        
        offset, V0 = 0, 0
        breeding_vol = 0
        for j, layer in enumerate(reactor.layers):
            offset += layer[0]
            R, Z, V = miller_model(reactor.R0, reactor.a, reactor.kappa, reactor.delta, offset)

            print(f"{reactor.name} {layer[1]} vol = {(V - V0)/1e6:.4f} m^3") # for debugging :)
            
            if layer[1].startswith('br'):
                breeding_vol += (V - V0)

            V0 = V

            # Use same color but different linestyle for layers
            linestyle = '-' # '--' if j % 2 == 0 else ':'
            ax.plot(R, Z, color=color, linestyle=linestyle, 
                   linewidth=1) # , alpha=0.7) # , label=f'{reactor.name} {layer[1]}')

        print(f"{reactor.name} breeding vol = {breeding_vol/1e6:.4f} m^3")
    
    # Set common properties
    ax.set_xlim(100, 1300)
    ax.set_ylim(-600, 600)
    ax.set_xlabel('radius (cm)', fontsize=10)
    ax.set_ylabel('height (cm)', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=8)
    plt.tight_layout()
    # plt.show()
    # plt.savefig('tokamaks_together_lowres.png', dpi=300, bbox_inches='tight', transparent=False)  






def miller_model_simple(t, R0, a, kappa, delta):
    """
    FOR V&V ONLY -- SUGGEST USING 'miller_model()' IN PROD CODE --ppark
    Calculate parametric coordinates from the Miller local equilibrium model.
    cf. R. L. Miller et al., doi.org/10.1063/1.872666
    cf. (Justin) Ball et al., arxiv.org/pdf/1510.08923
    
    Args:
        t  : array : parameter from 0 to 2pi
        R0 : float : major radius (cm)
        a  : float : minor radius (cm)
        kappa : float : elongation
        delta : float : triangularity
    
    Returns:
        R, Z : arrays of R and Z coordinates 
    """
    R = R0 + a * np.cos(t + delta * np.sin(t))
    Z = kappa * a * np.sin(t)

    # Ensure the contours are closed (first point = last point) for OpenMC
    if not (np.isclose(R[0], R[-1]) and np.isclose(Z[0], Z[-1])):
        R = np.append(R, R[0])
        Z = np.append(Z, Z[0])

    return R, Z


def miller_surface_area(R0, a, kappa, delta, n=1000):
    """
    Surface area of the toroidal surface obtained by revolving the Miller contour
    about the z-axis.
   
    Returns: 
      A : float : Surface area in (length units)^2.
    """
    t = np.linspace(0.0, 2*np.pi, n, endpoint=False)
    phi = t + delta*np.sin(t)

    R = R0 + a*np.cos(phi)
    dRdt = -a*np.sin(phi)*(1.0 + delta*np.cos(t))

    Z = kappa*a*np.sin(t)
    dZdt = kappa*a*np.cos(t)

    integrand = R*np.sqrt(dRdt**2 + dZdt**2)
    A = 2*np.pi*np.trapezoid(integrand, t)
    return A


def miller_model(R0, a, kappa, delta, extrude=0, calc_vol=True, n=100):
    """
    Calculate parametric coordinates from the Miller local equilibrium model.
    cf. R. L. Miller et al., doi.org/10.1063/1.872666
    cf. (Justin) Ball et al., arxiv.org/pdf/1510.08923
    
    Args:
        t  : array : parameter from 0 to 2pi
        R0 : float : major radius (cm)
        a  : float : minor radius (cm)
        kappa : float : elongation
        delta : float : triangularity
        extrude : float : thickness by which to extrude the boundary (for blanket layers)
    
    Returns:
        R, Z : arrays of R and Z coordinates 
    """
    t = np.linspace(0, 2*np.pi, n)
    # Original Miller contour
    R = R0 + a * np.cos(t + delta * np.sin(t))
    Z = kappa * a * np.sin(t)

    # Derivatives
    dR_dt = -a * np.sin(t + delta * np.sin(t)) * (1 + delta * np.cos(t))
    dZ_dt = kappa * a * np.cos(t)

    # Outward normal
    N_R = dZ_dt
    N_Z = -dR_dt
    N_mag = np.sqrt(N_R**2 + N_Z**2)
    N_R_unit, N_Z_unit = N_R/N_mag, N_Z/N_mag

    # Offset contour
    R_offset = R + extrude * N_R_unit
    Z_offset = Z + extrude * N_Z_unit

    # Close polygon
    if not (np.isclose(R_offset[0], R_offset[-1]) and np.isclose(Z_offset[0], Z_offset[-1])):
        R_offset = np.append(R_offset, R_offset[0])
        Z_offset = np.append(Z_offset, Z_offset[0])

    if calc_vol:
        # Shoelace area
        A = 0.5 * np.sum(R_offset[:-1]*Z_offset[1:] - R_offset[1:]*Z_offset[:-1])

        # Centroid in R
        Cx = (1/(6*A)) * np.sum((R_offset[:-1] + R_offset[1:]) *
                                (R_offset[:-1]*Z_offset[1:] - R_offset[1:]*Z_offset[:-1]))
        # Torus volume
        volume = 2*np.pi*Cx*abs(A)

        return R_offset, Z_offset, volume

    return R_offset, Z_offset # use this for OpenMC: list(zip(R_offset, Z_offset))


# def miller_volume(t, R0, a, kappa, delta, d):




if __name__ == '__main__':
    main()

