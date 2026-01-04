# import statements 
import os
# print export OPENMC_CROSS_SECTIONS=/Users/gretali/Desktop/research2025/FissileDependence/endfb-viii.0-hdf5/cross_sections.xml
# print openmc to get statepoint files
import numpy as np
import openmc

from Python.parameters import *
from Python.utilities  import *

# for this to work, python -m Python.wedgie_initialize

# set parameters ###################################################
hom = False ### heterogeneous - homogenize eurofer and breeder, but pack in discrete biso particles
fertile_bulk_density_kgm3 = 0.1
n_inner = 23.78
n_outer = 29.75 # [cm]
fertile_element = 'U'

current_dir = os.getcwd()
print(current_dir)

if hom == True: 
    save_folder = os.path.join(current_dir,"SpatialSelfShielding/HCPB_hom",f"{fertile_element}_{fertile_bulk_density_kgm3}")
else: 
    save_folder = os.path.join(current_dir,"SpatialSelfShielding/HCPB_het",f"{fertile_element}_{fertile_bulk_density_kgm3}")

os.makedirs(save_folder, exist_ok=True)

os.chdir(save_folder)



# wedge materials ##################################################
print("MATERIALS TIME!")
firstwall = openmc.Material(name='firstwall', temperature=TEMP_K)
firstwall.set_density('g/cm3',19.0)
firstwall.depletable = False
firstwall.add_element('W',1)

eurofer = openmc.Material(name='euroferwall', temperature=TEMP_K)
eurofer.depletable = False
eurofer.set_density('g/cm3', 7.8)
eurofer.add_element('Fe', 89.36, percent_type='wo')
eurofer.add_element('C' ,  0.11, percent_type='wo')
eurofer.add_element('Cr',  9.00, percent_type='wo')
eurofer.add_element('W' ,  1.10, percent_type='wo')
eurofer.add_element('Mn',  0.40, percent_type='wo')
eurofer.add_element('N' ,  0.03, percent_type='wo')

li4sio4 = openmc.Material(name='Li4SiO4', temperature=TEMP_K) 
li4sio4.set_density('g/cm3', DENSITY_LI4SIO4)  # normalized for ceramic porosity and 900 K (pure, room temp g/cm3 = 2.42)
li4sio4.add_elements_from_formula('Li4SiO4', enrichment_target='Li6', enrichment_type='ao', enrichment=ENRICH_HCPB) 
be = openmc.Material(name='Beryllium') 
be.set_density('g/cm3', DENSITY_BE)  # normalized from 1.85 for 900 K
be.add_element('Be', 1, percent_type='wo') 
he = openmc.Material(name='Helium') 
he.set_density('g/cm3', 0.004279) # Helium density at 900 K ~80 bar 
he.add_element('He', 1)  # if it's helium cooled...shouldn't it be colder? 

if fertile_element == 'U':
    kernel = openmc.Material(name='UO2')
    kernel.add_elements_from_formula('UO2', enrichment=ENRICH_U)
    kernel.set_density('g/cm3', DENSITY_UO2)  

elif fertile_element == 'Th': 
    kernel = openmc.Material(name='ThO2') 
    kernel.add_elements_from_formula('ThO2') 
    kernel.set_density('g/cm3', DENSITY_ThO2) 

sic = openmc.Material(name='SiC')
sic.add_elements_from_formula('SiC')
sic.set_density('g/cm3', 3.2)

vf_li4sio4 = 0.1304 ; vf_be = 0.3790 ; vf_eurofer = 0.1176 ; vf_he = 1 - (vf_li4sio4 + vf_be + vf_eurofer) # lu et al, og
vf_biso, vf_li4sio4_be, biso_per_cc_bv = calc_biso_breeder_vol_fracs(fertile_bulk_density_kgm3, fertile_isotope='U238')
print(f"vf_biso % {vf_biso*100}")
print(f"biso density per breeding VOLUME {biso_per_cc_bv}")
print(f"biso density per breeding REGION {biso_per_cc_bv/2}")
# vf_biso is not the total volume fraction of the entire blanket, just of the breeding fraction
# biso displace the lithium orthosilicate and beryllium pebbles

# integrate biso case ###############################################
if hom == True: 
    biso_hom = openmc.Material.mix_materials([kernel, sic], [BISO_KERNEL_VOL_FRAC, BISO_COAT_VOL_FRAC], 'vo') 
    homogeneous_breeder = openmc.Material.mix_materials([li4sio4, be, biso_hom], [vf_li4sio4/(vf_li4sio4+vf_be)*vf_li4sio4_be, vf_be/(vf_li4sio4+vf_be)*vf_li4sio4_be, vf_biso], 'vo') 
    homogeneous_blanket = openmc.Material.mix_materials([homogeneous_breeder, eurofer, he], [(vf_li4sio4+vf_be), vf_eurofer, vf_he], 'vo') # breeder - still original vol frac
    homogeneous_materials = openmc.Materials([firstwall, eurofer, homogeneous_blanket]) 
    homogeneous_materials.export_to_xml()

else: # heterogeneous case - must renormalize the blanket material without biso
    packing_fraction = vf_biso * (vf_li4sio4 + vf_be) # biso: of total blanket volume
    vf_li4sio4_disp = vf_li4sio4 * vf_li4sio4_be # li4sio4: fraction of total blanket volume, after displacement
    vf_be_disp = vf_be * vf_li4sio4_be # be: fraction of total blanket volume, after displacement 
    print(packing_fraction + vf_li4sio4_disp + vf_be_disp)
    print(f"biso fraction {vf_biso} wrt total breeder volume")
    print(f"packing fraction {packing_fraction} wrt total blanket volume")

    het_blanket_vfs=np.array([vf_li4sio4_disp, vf_be_disp, vf_eurofer, vf_he]) # renormalization, accounting for displacement 
    print(het_blanket_vfs/np.sum(het_blanket_vfs))

    heterogeneous_blanket = openmc.Material.mix_materials([li4sio4, be, eurofer, he], het_blanket_vfs/np.sum(het_blanket_vfs), 'vo') # breeder - homogenize without pebbles
    heterogeneous_materials = openmc.Materials([firstwall, eurofer, kernel, sic, heterogeneous_blanket])
    heterogeneous_materials.export_to_xml()

# geometry: create surfaces ##############################################
# WORK IN CENTIMETERS !
print("GEOMETRY SURFACE TIME!")
xmin_inner, xmax_inner = openmc.XPlane(x0=-n_inner/2, boundary_type='reflective'), openmc.XPlane(x0=n_inner/2, boundary_type='reflective')
ymin_inner, ymax_inner = openmc.YPlane(y0=-n_inner/2, boundary_type='reflective'), openmc.YPlane(y0=n_inner/2, boundary_type='reflective')

xmin_outer, xmax_outer = openmc.XPlane(x0=-n_outer/2, boundary_type='reflective'), openmc.XPlane(x0=n_outer/2, boundary_type='reflective')
ymin_outer, ymax_outer = openmc.YPlane(y0=-n_outer/2, boundary_type='reflective'), openmc.YPlane(y0=n_outer/2, boundary_type='reflective')


# set up trapezoid surfaces for frustrum
m = (n_outer - n_inner)/800
x_pos = openmc.Plane(a=1,b=0,c=-m,d=n_inner/2 - 50.2*m, name='frustum at +x', boundary_type='reflective')
x_neg = openmc.Plane(a=1,b=0,c=m,d=-n_inner/2 + 50.2*m, name='frustum at -x', boundary_type='reflective')
y_pos = openmc.Plane(a=0,b=1,c=-m,d=n_inner/2 - 50.2*m, name='frustum at +y', boundary_type='reflective')
y_neg = openmc.Plane(a=0,b=1,c=m,d=-n_inner/2 + 50.2*m, name='frustum at -y', boundary_type='reflective')

z_str1_l = openmc.ZPlane(z0=0, boundary_type='vacuum')
z_br1_l = openmc.ZPlane(z0=2.5)
z_str2_l = openmc.ZPlane(z0=47.5)
z_fw1_l = openmc.ZPlane(z0=50)
z_p_l = openmc.ZPlane(z0=50.2)
z_fw2_l = openmc.ZPlane(z0=450.2)
z_str3_l = openmc.ZPlane(z0=450.4)
z_br2_l = openmc.ZPlane(z0=452.9)
z_str4_l = openmc.ZPlane(z0=534.9)
z_end = openmc.ZPlane(z0=537.4, boundary_type='vacuum')

# geometry: create cells and fill w/ material ############################
print("GEOMETRY CELL TIME!")
inner_cuboid = +xmin_inner & -xmax_inner & +ymin_inner & -ymax_inner
outer_cuboid = +xmin_outer & -xmax_outer & +ymin_outer & -ymax_outer
plasma_trap = +x_neg & -x_pos & +y_neg & -y_pos

cell_str1 = openmc.Cell(name='eurofer wall LL', region=inner_cuboid & +z_str1_l & -z_br1_l, fill=eurofer)
if hom == True: 
    cell_br1 = openmc.Cell(name='breeder blanket L', region=inner_cuboid & +z_br1_l & -z_str2_l, fill=homogeneous_blanket)
    cell_br1 = [cell_br1]
else: 
    print("PACKING INNER")
    start_time = time.time()
    br1_centers = openmc.model.pack_spheres(radius=BISO_RADIUS, region=inner_cuboid & +z_br1_l & -z_str2_l,
                                            pf=packing_fraction, seed=5, initial_pf=0.1)
    print(len(br1_centers))
    br1_kernel_cells = []
    br1_coating_cells, br1_coating_surfaces = [], []
    for i, (x, y, z) in enumerate(br1_centers):
        kernel_surface = openmc.Sphere(x0=x, y0=y, z0=z, r=BISO_KERNEL_RADIUS)
        coating_surface = openmc.Sphere(x0=x, y0=y, z0=z, r=BISO_RADIUS)
        br1_coating_surfaces.append(coating_surface)

        kernel_cell = openmc.Cell(name=f'kernel_inner_{i}', region=-kernel_surface, fill=kernel)
        coating_cell = openmc.Cell(name=f'coating_inner_{i}', region=+kernel_surface & -coating_surface, fill=sic)
        br1_kernel_cells.append(kernel_cell)
        br1_coating_cells.append(coating_cell)
    
    swiss_cheese1_region = inner_cuboid & +z_br1_l & -z_str2_l
    for biso in br1_coating_surfaces:
        swiss_cheese1_region &= +biso # define as the complement, make it faster? 
    cell_swiss_cheese1 = openmc.Cell(name='inner breeder w/ biso scooped out', region=swiss_cheese1_region, 
                                    fill=heterogeneous_blanket)
    cell_br1 = [cell_swiss_cheese1] + br1_kernel_cells + br1_coating_cells
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Time taken: {elapsed_time} seconds")

cell_str2 = openmc.Cell(name='eurofer wall LR', region=inner_cuboid & +z_str2_l & -z_fw1_l, fill=eurofer)
cell_fw1 = openmc.Cell(name='first wall L', region=inner_cuboid & +z_fw1_l & -z_p_l, fill=firstwall)
cell_p = openmc.Cell(name='plasma chamber', region=plasma_trap & +z_p_l & -z_fw2_l)
cell_fw2 = openmc.Cell(name='first wall R', region=outer_cuboid & +z_fw2_l & -z_str3_l, fill=firstwall)
cell_str3 = openmc.Cell(name='eurofer wall RL', region=outer_cuboid & +z_str3_l & -z_br2_l, fill=eurofer)
if hom==True:
    cell_br2 = openmc.Cell(name='breeder blanket R', region=outer_cuboid & +z_br2_l & -z_str4_l, fill=homogeneous_blanket)
    cell_br2 = [cell_br2] # to put it in as a list 
else: 
    print("PACKING OUTER")
    start_time = time.time()
    br2_centers = openmc.model.pack_spheres(radius=BISO_RADIUS, region=outer_cuboid & +z_br2_l & -z_str4_l,
                                            pf=packing_fraction, seed=1, initial_pf=0.1)
    print(len(br2_centers))
    br2_kernel_cells = []
    br2_coating_cells, br2_coating_surfaces = [], []
    for i, (x, y, z) in enumerate(br2_centers):
        kernel_surface = openmc.Sphere(x0=x, y0=y, z0=z, r=BISO_KERNEL_RADIUS)
        coating_surface = openmc.Sphere(x0=x, y0=y, z0=z, r=BISO_RADIUS)
        br2_coating_surfaces.append(coating_surface)

        kernel_cell = openmc.Cell(name=f'kernel_outer_{i}', region=-kernel_surface, fill=kernel)
        coating_cell = openmc.Cell(name=f'coating_outer_{i}', region=+kernel_surface & -coating_surface, fill=sic)
        br2_kernel_cells.append(kernel_cell)
        br2_coating_cells.append(coating_cell)
    
    swiss_cheese2_region = outer_cuboid & +z_br2_l & -z_str4_l
    for biso in br2_coating_surfaces:
        swiss_cheese2_region &= +biso 
    cell_swiss_cheese2 = openmc.Cell(name='inner breeder w/ biso scooped out', region=swiss_cheese2_region, 
                                    fill=heterogeneous_blanket)
    cell_br2 = [cell_swiss_cheese2] + br2_kernel_cells + br2_coating_cells
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Time taken: {elapsed_time} seconds")

cell_str4 = openmc.Cell(name='eurofer wall RR', region=outer_cuboid & +z_str4_l & -z_end, fill=eurofer)

cells_all = [cell_str1,] + cell_br1 + [cell_str2, cell_fw1, cell_p, cell_fw2, cell_str3,] + cell_br2 + [cell_str4]
geometry = openmc.Geometry(openmc.Universe(cells=cells_all))
geometry.export_to_xml()

# # PLOTTING ############################
# print("PLOTTING TIME!")
# plot = openmc.Plot()
# plot.filename = 'geometry_plot_xz'
# plot.origin = (0, 0, 268.7)
# plot.width = (60, 550)
# plot.pixels = (300, 600)
# plot.basis = 'xz'           
# plot.color_by = 'material'    

# plots = openmc.Plots([plot])
# plots.export_to_xml()
# openmc.plot_geometry()

# NEUTRON SOURCE ####################
print("SETTINGS TIME!")
source = openmc.IndependentSource()
source.space = openmc.stats.Point((0.0, 0.0, 250.2))
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([14.0e6], [1.0])  # 14 MeV                                                
source.particle = 'neutron'

settings = openmc.Settings()
settings.source = source
settings.run_mode = 'fixed source'
settings.particles = int(1e5)
settings.batches = 10

settings.export_to_xml()

# FIND URANIUM MATERIALS ####################
u_mat_ids = []

if hom == True: 
    materials = homogeneous_materials 

else: 
    materials = heterogeneous_materials

for mat in materials:
    for nuclide, _, _ in mat.nuclides:
        if nuclide == 'U238':
            u_mat_ids.append(mat.id)

u_cell_ids = [cell.id for cell in cells_all if 'kernel' in cell.name]
print(f"Number of U kernel cells: {len(u_cell_ids)}")
print(f"density of biso {len(u_cell_ids) / (n_inner**2 * 45 + n_outer**2 * 85)}")

# TALLIES #######################
print("TALLIES TIME!")
tallies = openmc.Tallies()

u238_ngamma = openmc.Tally(name='U238_n_gamma')
u238_ngamma.filters = [openmc.MaterialFilter(u_mat_ids)]
# u238_ngamma.filters = [openmc.CellFilter(u_cell_ids)]
u238_ngamma.nuclides = ['U238']
u238_ngamma.scores = ['(n,gamma)']

tallies.append(u238_ngamma)
tallies.export_to_xml()