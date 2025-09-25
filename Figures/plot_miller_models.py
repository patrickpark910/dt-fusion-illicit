import numpy as np
import matplotlib.pyplot as plt
matplotlib.use("Agg")  # non-GUI backend
import os
os.makedirs("figures", exist_ok=True)


"""
This file is meant to be standalone (execute this by itself, independent of main.py)
"""

class Specs:
    def __init__(self, name:str, R0:float, a:float, kappa:float, delta:float, outboardlayers:list, inboardlayers:list=None):
        self.name   = name
        self.R0, self.a, self.kappa, self.delta = R0, a, kappa, delta
        self.aspect = self.R0 / self.a
        self.outboardlayers = outboardlayers
        self.inboardlayers = inboardlayers if inboardlayers is not None else outboardlayers

        self.surf_area = miller_surface_area(R0, a, kappa, delta)
        print(f"{self.name} plasma-facing surface area: {self.surf_area/1e4:.4f} m^2")

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
                           620, 207, 1.72, 0.4 , [(0.2,'fw'), (0.4,'fwf'), (2,'fwc'), (0.4,'fwb'), (22.5,'br1'), (3.2,'d1'), (21.0,'br2'), (3.2,'d2'), (21.0,'br3'), (8.0,'im'), (1.5, 'bp')], 
                           [(0.2,'fw'), (0.4,'fwf'), (2,'fwc'), (0.4,'fwb'), (22.5,'br1'), (3.2,'d1'), (21.0,'br2'), (8.0,'im'), (1.5, 'bp')])
    
    EU_DEMO        = Specs('EU DEMO',
                           907.2, 292.7, 1.59, 0.33, [(0.2,'fw'), (0.4,'st1'), (2,'br1'), (0.4,'st2'), (22.5,'br2'), (3.2,'st3'), (21.0,'br3'), (3.2,'st4'), (21.0,'br4'), (8.0,'st5')]) 

    # plot_separate([ARC_Sorbom, ARC_Ball_code, Emma_FLiBe, Emma_LL, EU_DEMO]) # , ITER_LL, EU_DEMO_PB]
    plot_together([ARC_Sorbom, ARC_Ball_code, ARC_Ball_paper, Emma_FLiBe, Emma_LL, EU_DEMO]) # , ITER_LL, EU_DEMO_PB]


def plot_separate(reactors_to_plot, n=10000):

    fig, ax = plt.subplots(1, len(reactors_to_plot), figsize=(7*len(reactors_to_plot), 6))

    # Create parameter array
    i, t = 0, np.linspace(0, 2*np.pi, n)

    # Define colors for different reactors
    colors = plt.cm.Set1(np.arange(len(reactors_to_plot))) 

    for i, reactor in enumerate(reactors_to_plot):
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
        R, Z, V = miller_offset(t, reactor.R0, reactor.a, reactor.kappa, reactor.delta, 0)
        ax[i].plot(R, Z, '-', color=color, linewidth=1) # , label='Original Miller D-shape')
        
        # --- Case 1: asymmetric (separate inboard vs outboard) ---
        if reactor.inboardlayers is not None and reactor.inboardlayers != reactor.outboardlayers:
            # Outboard
            offset_out = 0
            for d, label in enumerate(reactor.outboardlayers):
                offset_out += label[0]
                (R_in, Z_in, _), (R_out, Z_out, _) = miller_offset_split(
                    t, reactor.R0, reactor.a, reactor.kappa, reactor.delta,
                    d_in=0, d_out=offset_out, calc_vol=False
                )
                ax[i].plot(R_out, Z_out, '-', color=color, linewidth=1, label=f'OB {label[1]}')

            # Inboard
            offset_in = 0
            for d, label in enumerate(reactor.inboardlayers):
                offset_in += label[0]
                (R_in, Z_in, _), (R_out, Z_out, _) = miller_offset_split(
                    t, reactor.R0, reactor.a, reactor.kappa, reactor.delta,
                    d_in=offset_in, d_out=0, calc_vol=False
                )
                ax[i].plot(R_in, Z_in, '--', color=color, linewidth=1, label=f'IB {label[1]}')

        # --- Case 2: symmetric (use old loop) ---
        else:
            offset = 0
            for d, label in enumerate(reactor.outboardlayers):
                offset += label[0]
                R_offset, Z_offset, _ = miller_offset(t, reactor.R0, reactor.a, reactor.kappa, reactor.delta, offset)
                ax[i].plot(R_offset, Z_offset, '-', color=color, linewidth=1, label=f'{label[1]}')

        ax[i].legend(loc='best', fontsize=8)

    plt.tight_layout()
    plt.show()
    plt.savefig('tokamaks_sym_asym.png', dpi=300, bbox_inches='tight', transparent=False)  


def plot_together(reactors_to_plot, n=10000):
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))  # Single subplot
    
    # Create parameter array
    t = np.linspace(0, 2*np.pi, n)
    
    # Define colors for different reactors
    colors = plt.cm.Set1(np.arange(len(reactors_to_plot))) 
    
    for i, reactor in enumerate(reactors_to_plot):
        color = colors[i]
        
        # Plasma shape
        R, Z, VPlasma = miller_offset(t, reactor.R0, reactor.a, reactor.kappa, reactor.delta, 0)
        ax.plot(R, Z, color=color, linewidth=1, label=f'{reactor.name}')
        
        breeding_vol = 0  # initialize here for both cases

        # Asymmetric inboard and outboard layers
        if reactor.inboardlayers is not None and reactor.inboardlayers != reactor.outboardlayers:
            # --- Outboard layers ---
            offset_out, V_out_prev = 0, VPlasma
            for d, tag in enumerate(reactor.outboardlayers):
                offset_out += tag[0]
                (R_in, Z_in, V_in), (R_out, Z_out, V_out) = miller_offset_split(
                    t, reactor.R0, reactor.a, reactor.kappa, reactor.delta,
                    d_in=0, d_out=offset_out, calc_vol=True
                )
                dV_out = V_out - V_out_prev
                print(f"{reactor.name}  OUT {tag[1]} vol = {dV_out/1e6:.4f} m^3")
                if tag[1].startswith('br'):
                    breeding_vol += dV_out
                V_out_prev = V_out
                ax.plot(R_out, Z_out, color=color, linestyle='-', linewidth=1)

            # --- Inboard layers ---
            offset_in, V_in_prev = 0, VPlasma
            for d, tag in enumerate(reactor.inboardlayers):
                offset_in += tag[0]
                (R_in, Z_in, V_in), (R_out, Z_out, V_out) = miller_offset_split(
                    t, reactor.R0, reactor.a, reactor.kappa, reactor.delta,
                    d_in=offset_in, d_out=0, calc_vol=True
                )
                dV_in = V_in - V_in_prev
                print(f"{reactor.name} IN {tag[1]} vol = {dV_in/1e6:.4f} m^3")
                if tag[1].startswith('br'):
                    breeding_vol += dV_in 
                V_in_prev = V_in
                ax.plot(R_in, Z_in, color=color, linestyle='-', linewidth=1)
        else:
            # Symmetric layers around full circumference
            offset, V0 = 0, 0
            breeding_vol = 0
            for j, layer in enumerate(reactor.outboardlayers):
                offset += layer[0]
                R, Z, V = miller_offset(t, reactor.R0, reactor.a, reactor.kappa, reactor.delta, offset)

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
    plt.show()
    plt.savefig('tokamaks_together_sym_asym.png', dpi=300, bbox_inches='tight', transparent=False)  






def miller_model(t, R0, a, kappa, delta):
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

def split_inboard_outboard(R, Z):
    # Find indices of top and bottom of our D
    i_top = np.argmax(Z)
    i_bot = np.argmin(Z)

    # Path 1: top to bottom 
    if i_top < i_bot:
        R1, Z1 = R[i_top:i_bot+1], Z[i_top:i_bot+1]
        R2, Z2 = np.concatenate([R[i_bot:], R[:i_top+1]]), np.concatenate([Z[i_bot:], Z[:i_top+1]])
    else: # Path 2: top to bottom
        R1, Z1 = R[i_top:], Z[i_top:]
        R1, Z1 = np.concatenate([R1, R[:i_bot+1]]), np.concatenate([Z1, Z[:i_bot+1]])
        R2, Z2 = R[i_bot:i_top+1], Z[i_bot:i_top+1]

    # Now decide which is inboard vs outboard
    # Outboard has larger R (to the right)
    if np.mean(R1) > np.mean(R2):
        R_out, Z_out = R1, Z1
        R_in, Z_in   = R2, Z2
    else:
        R_out, Z_out = R2, Z2
        R_in, Z_in   = R1, Z1

    return (R_in, Z_in), (R_out, Z_out)


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


def miller_offset(t, R0, a, kappa, delta, d, calc_vol=True):

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
    R_offset = R + d * N_R_unit
    Z_offset = Z + d * N_Z_unit

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

    return R_offset, Z_offset

def miller_offset_split(t, R0, a, kappa, delta, d_in, d_out, calc_vol=True):
    """
    Return separate offset curves for inboard and outboard sides.
    """
    # Base Miller contour
    R = R0 + a * np.cos(t + delta * np.sin(t))
    Z = kappa * a * np.sin(t)

    # Derivatives for normal vector
    dR_dt = -a * np.sin(t + delta * np.sin(t)) * (1 + delta * np.cos(t))
    dZ_dt = kappa * a * np.cos(t)

    N_R = dZ_dt
    N_Z = -dR_dt
    N_mag = np.sqrt(N_R**2 + N_Z**2)
    N_R_unit, N_Z_unit = N_R / N_mag, N_Z / N_mag

    # Split curve into inboard vs outboard
    i_top = np.argmax(Z)
    i_bot = np.argmin(Z)

    if i_top < i_bot:
        idx_out = np.arange(i_top, i_bot+1)
        idx_in  = np.concatenate([np.arange(i_bot, len(R)), np.arange(0, i_top+1)])
    else:
        idx_out = np.concatenate([np.arange(i_top, len(R)), np.arange(0, i_bot+1)])
        idx_in  = np.arange(i_bot, i_top+1)

    # Apply different offsets
    R_out = R[idx_out] + d_out * N_R_unit[idx_out]
    Z_out = Z[idx_out] + d_out * N_Z_unit[idx_out]

    R_in  = R[idx_in]  + d_in  * N_R_unit[idx_in]
    Z_in  = Z[idx_in]  + d_in  * N_Z_unit[idx_in]

    if calc_vol:
        def polygon_volume(Rc, Zc):
            # Close polygon by adding the cut line (straight vertical from bottom to top)
            Rc_closed = np.append(Rc, Rc[0])
            Zc_closed = np.append(Zc, Zc[0])

            # Shoelace area
            A = 0.5 * np.sum(Rc_closed[:-1]*Zc_closed[1:] - Rc_closed[1:]*Zc_closed[:-1])

            # Centroid in R
            Cx = (1/(6*A)) * np.sum((Rc_closed[:-1] + Rc_closed[1:]) *
                                    (Rc_closed[:-1]*Zc_closed[1:] - Rc_closed[1:]*Zc_closed[:-1]))

            # Toroidal volume
            return 2*np.pi*Cx*abs(A)

        V_out = polygon_volume(R_out, Z_out)
        V_in  = polygon_volume(R_in,  Z_in)

        return (R_in, Z_in, V_in), (R_out, Z_out, V_out)

    return (R_in, Z_in), (R_out, Z_out)




# def miller_volume(t, R0, a, kappa, delta, d):




if __name__ == '__main__':
    main()


