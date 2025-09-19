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

    def __repr__(self):
        return (f"TokamakReactor(R0= {self.R0:.4f} m, a={self.a:.4f} m, "
                f"kappa= {self.kappa:.4f}, delta={self.delta:.4f})")

    def calc_volumes(self):
        pass

def main():

    ARC_Sorbom     = Specs('ARC (Sorbom 15)',     
                           330, 113, 1.84, 0.45, [(0.3,'fw'), (1,'st1'), (2,'br1'), (3,'st2'), (100,'br2'), (3,'st3')])
    ARC_Ball_code  = Specs('ARC-class (Ball 25)',
                           400, 120, 1.5 , 0.5 , [(0.3,'fw'), (1,'st1'), (2,'br1'), (3,'st2'), (100,'br2'), (3,'st3')])
    # ARC_Ball_paper = Specs(400, 100, 1.6 , 0.5 , [(0.3,'fw'), (1,'st1'), (2,'br1'), (3,'st2'), (100,'br2'), (3,'st3')])
    Emma_FLiBe     = Specs('Emma FLiBe',
                           600, 200, 1.72, 0.4 , [(0.3,'fw'), (1,'st1'), (2,'br1'), (3,'st2'), (100,'br2'), (3,'st3')])
    Emma_LL        = Specs('Emma LL/PB',
                           620, 207, 1.72, 0.4 , [(0.2,'fw'), (0.4,'st1'), (2,'br1'), (0.4,'st2'), (22.5,'br2'), (3.2,'st3'), (21.0,'br3'), (3.2,'st4'), (21.0,'br4'), (8.0,'st5')])
    
    EU_DEMO        = Specs('EU DEMO',
                           907.2, 292.7, 1.59, 0.33, [(0.2,'fw'), (0.4,'st1'), (2,'br1'), (0.4,'st2'), (22.5,'br2'), (3.2,'st3'), (21.0,'br3'), (3.2,'st4'), (21.0,'br4'), (8.0,'st5')]) 

    plot_separate([ARC_Sorbom, ARC_Ball_code, Emma_FLiBe, Emma_LL, EU_DEMO]) # , ITER_LL, EU_DEMO_PB]
    plot_together([ARC_Sorbom, ARC_Ball_code, Emma_FLiBe, Emma_LL, EU_DEMO]) # , ITER_LL, EU_DEMO_PB]


def plot_separate(reactors_to_plot):

    fig, ax = plt.subplots(1, len(reactors_to_plot), figsize=(7*len(reactors_to_plot), 6))

    # Create parameter array
    i, t = 0, np.linspace(0, 2*np.pi, 100)

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
        R, Z = miller_model(t, reactor.R0, reactor.a, reactor.kappa, reactor.delta)
        ax[i].plot(R, Z, '-', color=color, linewidth=1) # , label='Original Miller D-shape')
        
        offset = 0
        for layer in reactor.layers:
            offset += layer[0]
            R_offset, Z_offset = miller_offset(t, reactor.R0, reactor.a, reactor.kappa, reactor.delta, offset)
            ax[i].plot(R_offset, Z_offset, '-', color=color, linewidth=1, label=f'{layer[1]}')

        i+=1

    plt.tight_layout()
    # plt.show()
    plt.savefig('tokamaks_separate_lowres.png', 
            dpi=300,            
            bbox_inches='tight',
            transparent=False)  



def plot_together(reactors_to_plot):
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))  # Single subplot
    
    # Create parameter array
    t = np.linspace(0, 2*np.pi, 100)
    
    # Define colors for different reactors
    colors = plt.cm.Set1(np.arange(len(reactors_to_plot))) 
    
    for i, reactor in enumerate(reactors_to_plot):
        color = colors[i]
        
        # Plasma shape
        R, Z = miller_model(t, reactor.R0, reactor.a, reactor.kappa, reactor.delta)
        ax.plot(R, Z, color=color, linewidth=1, label=f'{reactor.name}')
        
        offset = 0
        for j, layer in enumerate(reactor.layers):
            offset += layer[0]
            R_offset, Z_offset = miller_offset(t, reactor.R0, reactor.a, reactor.kappa, reactor.delta, offset)
            # Use same color but different linestyle for layers
            linestyle = '-' # '--' if j % 2 == 0 else ':'
            ax.plot(R_offset, Z_offset, color=color, linestyle=linestyle, 
                   linewidth=1) # , alpha=0.7) # , label=f'{reactor.name} {layer[1]}')
    
    # Set common properties
    ax.set_xlim(100, 1300)
    ax.set_ylim(-600, 600)
    ax.set_xlabel('radius (cm)', fontsize=10)
    ax.set_ylabel('height (cm)', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=8)
    plt.tight_layout()
    # plt.show()
    plt.savefig('tokamaks_together_lowres.png', 
            dpi=300,            
            bbox_inches='tight',
            transparent=False)  





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




def miller_offset(t, R0, a, kappa, delta, d):
    """
    Calculate a new contour that is precisely [d] meters extruded from a shape generated by miller_model().
    
    Args:
        t  : array : parameter from 0 to 2pi
        R0 : float : major radius (m)
        a  : float : minor radius (m)
        kappa : float : elongation
        delta : float : triangularity
        d  : float : offset distance (m)
    
    Returns:
        R_offset, Z_offset : arrays of offset R and Z coordinates
    """
    # Original shape
    R = R0 + a * np.cos(t + delta * np.sin(t))
    Z = kappa * a * np.sin(t)
    
    # Derivatives
    dR_dt = -a * np.sin(t + delta * np.sin(t)) * (1 + delta * np.cos(t))
    dZ_dt = kappa * a * np.cos(t)
    
    # Normal vector components (not normalized)
    N_R = dZ_dt   # kappa * a * cos(t)
    N_Z = -dR_dt  # a * sin(t + delta*sin(t)) * (1 + delta*cos(t))
    
    # Magnitude of normal vector
    N_mag = np.sqrt(N_R**2 + N_Z**2)
    
    # Unit normal components
    N_R_unit = N_R / N_mag
    N_Z_unit = N_Z / N_mag
    
    # Offset coordinates
    R_offset = R + d * N_R_unit
    Z_offset = Z + d * N_Z_unit
    
    # Ensure the contours are closed (first point = last point) for OpenMC
    if not (np.isclose(R_offset[0], R_offset[-1]) and np.isclose(Z_offset[0], Z_offset[-1])):
        R_offset = np.append(R_offset, R_offset[0])
        Z_offset = np.append(Z_offset, Z_offset[0])

    return R_offset, Z_offset




if __name__ == '__main__':
    main()




    # # Plot 2: Zoomed in view to see the offset clearly
    # # Find a region to zoom (e.g., the outer midplane)
    # zoom_center_R = R0 + a
    # zoom_center_Z = 0
    # zoom_window = 0.5
    # 
    # ax2.plot(R_orig, Z_orig, 'b-', linewidth=2, label='Original')
    # ax2.plot(R_offset, Z_offset, 'r-', linewidth=2, label=f'Offset ({d*100:.1f} cm)')
    # 
    # # Add some normal vectors for visualization
    # n_vectors = 20
    # t_vectors = np.linspace(0, 2*np.pi, n_vectors)
    # for t_val in t_vectors:
    #     R_pt, Z_pt = miller_d_shape(t_val, R0, a, kappa, delta)
    #     R_off_pt, Z_off_pt = miller_d_shape_offset(t_val, R0, a, kappa, delta, d)
    #     ax2.arrow(R_pt, Z_pt, R_off_pt - R_pt, Z_off_pt - Z_pt,
    #               head_width=0.02, head_length=0.01, fc='gray', ec='gray', alpha=0.5)
    # 
    # ax2.set_xlim([zoom_center_R - zoom_window, zoom_center_R + zoom_window])
    # ax2.set_ylim([zoom_center_Z - zoom_window, zoom_center_Z + zoom_window])
    # ax2.set_xlabel('R (m)', fontsize=12)
    # ax2.set_ylabel('Z (m)', fontsize=12)
    # ax2.set_title('Zoomed View with Normal Vectors', fontsize=14)
    # ax2.grid(True, alpha=0.3)
    # ax2.axis('equal')
    # ax2.legend(loc='best', fontsize=10)
    # 
    # plt.tight_layout()
    # 
    # # Print some information
    # print(f"Miller D-shape parameters:")
    # print(f"  Major radius R0 = {R0} m")
    # print(f"  Minor radius a = {a} m")
    # print(f"  Elongation κ = {kappa}")
    # print(f"  Triangularity δ = {delta}")
    # print(f"  Offset distance = {d*100} cm")
    # print(f"\nAspect ratio = {R0/a:.2f}")
    # print(f"Plasma height ≈ {2*kappa*a:.2f} m")
    # print(f"Plasma width ≈ {2*a:.2f} m")
    # 
    # # Calculate and print the perimeter (approximate)
    # dR = np.diff(R_orig)
    # dZ = np.diff(Z_orig)
    # perimeter_orig = np.sum(np.sqrt(dR**2 + dZ**2))
    # dR_off = np.diff(R_offset)
    # dZ_off = np.diff(Z_offset)
    # perimeter_offset = np.sum(np.sqrt(dR_off**2 + dZ_off**2))
    # print(f"\nApproximate perimeter:")
    # print(f"  Original: {perimeter_orig:.3f} m")
    # print(f"  Offset: {perimeter_offset:.3f} m")
    # print(f"  Difference: {perimeter_offset - perimeter_orig:.3f} m")
    # 
    # plt.show()