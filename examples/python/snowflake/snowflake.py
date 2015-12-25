import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import openmc

def build_inputs():
    """Build OpenMC inputs and export them to XML"""
    # OpenMC simulation parameters
    batches = 15
    inactive = 5
    particles = 10000

    ####################
    ## Materials
    ####################

    # Instantiate some Nuclides
    h1 = openmc.Nuclide('H-1')
    o16 = openmc.Nuclide('O-16')
    u235 = openmc.Nuclide('U-235')

    # Instantiate some Materials and register the appropriate Nuclides
    moderator = openmc.Material(material_id=41, name='moderator')
    moderator.set_density('g/cc', 2.0)
    moderator.add_nuclide(h1, 2.)
    moderator.add_nuclide(o16, 1.)
    moderator.add_s_alpha_beta('HH2O', '71t')

    fuel = openmc.Material(material_id=40, name='fuel')
    fuel.set_density('g/cc', 4.5)
    fuel.add_nuclide(u235, 1.)

    # Instantiate a MaterialsFile, register all Materials, and export to XML
    materials_file = openmc.MaterialsFile()
    materials_file.default_xs = '71c'
    materials_file.add_materials([moderator, fuel])
    materials_file.export_to_xml()


    ####################
    ## Geometry
    ####################

    angle_30 = np.deg2rad(30)
    angle_60 = np.deg2rad(60)

    # Central trunk.
    s1 = openmc.XPlane(x0=-0.5)
    s2 = openmc.XPlane(x0=0.5)
    s3 = openmc.YPlane(y0=15.0)
    central_trunk = +s1 & -s2 & -s3

    # Outer branches.
    s11 = openmc.Plane(A=np.cos(-angle_60), B=np.sin(-angle_60), C=0., D=-10)
    s12 = openmc.Plane(A=np.cos(-angle_60), B=np.sin(-angle_60), C=0., D=-11)
    s13 = openmc.Plane(A=np.cos(angle_30), B=np.sin(angle_30), C=0., D=8.5)
    branch1 = -s11 & +s12 & -s13 & +s2

    s14 = openmc.Plane(A=np.cos(angle_60), B=np.sin(angle_60), C=0., D=10)
    s15 = openmc.Plane(A=np.cos(angle_60), B=np.sin(angle_60), C=0., D=11)
    s16 = openmc.Plane(A=np.cos(-angle_30), B=np.sin(-angle_30), C=0., D=-8.5)
    branch2 = +s14 & -s15 & +s16 & -s1

    # Middle branches
    s21 = openmc.Plane(A=np.cos(-angle_60), B=np.sin(-angle_60), C=0., D=-8)
    s22 = openmc.Plane(A=np.cos(-angle_60), B=np.sin(-angle_60), C=0., D=-9)
    s23 = openmc.Plane(A=np.cos(angle_30), B=np.sin(angle_30), C=0., D=8.5)
    branch3 = -s21 & +s22 & -s23 & +s2

    s24 = openmc.Plane(A=np.cos(angle_60), B=np.sin(angle_60), C=0., D=8)
    s25 = openmc.Plane(A=np.cos(angle_60), B=np.sin(angle_60), C=0., D=9)
    s26 = openmc.Plane(A=np.cos(-angle_30), B=np.sin(-angle_30), C=0., D=-8.5)
    branch4 = +s24 & -s25 & +s26 & -s1

    # Inner branches
    s31 = openmc.Plane(A=np.cos(-angle_60), B=np.sin(-angle_60), C=0., D=-6)
    s32 = openmc.Plane(A=np.cos(-angle_60), B=np.sin(-angle_60), C=0., D=-7)
    s33 = openmc.Plane(A=np.cos(angle_30), B=np.sin(angle_30), C=0., D=7)
    branch5 = -s31 & +s32 & -s33 & +s2

    s34 = openmc.Plane(A=np.cos(angle_60), B=np.sin(angle_60), C=0., D=6)
    s35 = openmc.Plane(A=np.cos(angle_60), B=np.sin(angle_60), C=0., D=7)
    s36 = openmc.Plane(A=np.cos(-angle_30), B=np.sin(-angle_30), C=0., D=-7)
    branch6 = +s34 & -s35 & +s36 & -s1

    # Center ring
    s41 = openmc.Plane(A=np.cos(angle_60), B=np.sin(angle_60), C=0., D=4)
    s42 = openmc.Plane(A=np.cos(angle_60), B=np.sin(angle_60), C=0., D=5)
    ring1 = +s41 & -s42 & +s2

    s43 = openmc.Plane(A=np.cos(-angle_60), B=np.sin(-angle_60), C=0., D=-4)
    s44 = openmc.Plane(A=np.cos(-angle_60), B=np.sin(-angle_60), C=0., D=-5)
    ring2 = -s43 & +s44 & -s1

    # Create the tree universe
    flake_c = openmc.Cell()
    flake_c.region = (central_trunk | branch1 | branch2 | branch3 | branch4
         | branch5 | branch6 | ring1 | ring2)
    flake_c.fill = fuel
    outside_c = openmc.Cell()
    outside_c.region = ~flake_c.region
    outside_c.fill = moderator
    tree_univ = openmc.Universe()
    tree_univ.add_cells([flake_c, outside_c])

    # Make azimuthal divisions
    s101 = openmc.XPlane(x0=0.)
    s102 = openmc.Plane(A=np.cos(-angle_60), B=np.sin(-angle_60), C=0., D=0.)
    s103 = openmc.Plane(A=np.cos(angle_60), B=np.sin(angle_60), C=0., D=0.)
    s104 = openmc.ZCylinder(R=16.0)
    s104.boundary_type = 'vacuum'

    # Place 6 trees around the flake
    tree1 = openmc.Cell()
    tree1.region = +s101 & -s102 & -s104
    tree1.fill = tree_univ
    tree1.rotation = [0.0, 0.0, -30.0]

    tree2 = openmc.Cell()
    tree2.region = +s102 & +s103 & -s104
    tree2.fill = tree_univ
    tree2.rotation = [0.0, 0.0, -90.0]

    tree3 = openmc.Cell()
    tree3.region = -s103 & +s101 & -s104
    tree3.fill = tree_univ
    tree3.rotation = [0.0, 0.0, -150.0]

    tree4 = openmc.Cell()
    tree4.region = -s101 & +s102 & -s104
    tree4.fill = tree_univ
    tree4.rotation = [0.0, 0.0, -210.0]

    tree5 = openmc.Cell()
    tree5.region = -s102 & -s103 & -s104
    tree5.fill = tree_univ
    tree5.rotation = [0.0, 0.0, -270.0]

    tree6 = openmc.Cell()
    tree6.region = +s103 & -s101 & -s104
    tree6.fill = tree_univ
    tree6.rotation = [0.0, 0.0, -330.0]

    root = openmc.Universe(universe_id=0, name='root universe')
    root.add_cells([tree1, tree2, tree3, tree4, tree5, tree6])

    # Create and export the geometry and geometry file.
    geometry = openmc.Geometry()
    geometry.root_universe = root

    geometry_file = openmc.GeometryFile()
    geometry_file.geometry = geometry
    geometry_file.export_to_xml()


    ####################
    ## Settings
    ####################

    # Instantiate a SettingsFile, set all runtime parameters, and export to XML
    settings_file = openmc.SettingsFile()
    settings_file.batches = batches
    settings_file.inactive = inactive
    settings_file.particles = particles
    settings_file.set_source_space('box', [-4, -4, -4, 4, 4, 4])
    settings_file.output = {'summary': True}
    settings_file.export_to_xml()


    ####################
    ## Plots
    ####################

    # Create a material plot
    plot = openmc.Plot()
    plot.filename = 'matplot'
    plot.width = (33, 33)
    plot.pixels = (400, 400)
    plot.color = 'mat'

    plotfile = openmc.PlotsFile()
    plotfile.add_plot(plot)
    plotfile.export_to_xml()


    ####################
    ## Tallies
    ####################

    # Instantiate tally filters and meshes
    e_bins = [0.0, 5.8e-8, 1.4e-7, 2.8e-7, 6.25e-7, 4e-6, 0.00553, 0.821, 20.0]
    energy_filter = openmc.Filter(type='energy', bins=e_bins)
    mesh = openmc.Mesh()
    mesh.dimension = (100, 100)
    mesh.lower_left = (-16, -16)
    mesh.upper_right = (16, 16)
    mesh_filter = openmc.Filter(type='mesh', bins=[mesh.id])

    # Instantiate the tally
    tally = openmc.Tally(tally_id=1)
    tally.add_filter(energy_filter)
    tally.add_filter(mesh_filter)
    tally.add_score('flux')
    tally.add_score('total')
    tally.add_score('fission')

    # Instantiate a tally file, register the tally, and export to XML
    tallies_file = openmc.TalliesFile()
    tallies_file.add_tally(tally)
    tallies_file.add_mesh(mesh)
    tallies_file.export_to_xml()


def convert_e_to_g(df):
    """Convert energy bins to group numbers in a dataframe"""
    # Find all of the possible energy bins.
    e_bins = []
    for x in df['energy [MeV]']:
        if x in e_bins: continue
        if type(x) == float:
            if math.isnan(x): continue
        e_bins.append(x)

    # Find the indices that would sort the energy bins in descending order
    # (i.e. group numbers)
    vals = [float(x.split(' - ')[0].strip('(')) for x in e_bins]
    sorted_vals = sorted(vals)[::-1]
    inds = [sorted_vals.index(x) for x in vals]

    # Make a map from energy bins to group numbers.
    e2g_map = {}
    for i in range(len(e_bins)): e2g_map[e_bins[i]] = inds[i]

    # Make a subfunction to replace energy bins with group numbers in a
    # dataframe.
    def energy_to_group(df, old_label, new_label):
        groups = []
        for ebin in df[old_label]:
            if type(ebin) == float:
                assert math.isnan(ebin)
                groups.append(ebin)
            else:
                groups.append(e2g_map[ebin])
        df[new_label] = groups
        ind = df.columns.tolist().index(old_label)
        del df[old_label]
        return df

    # Make the conversion and return.
    df = energy_to_group(df, 'energy [MeV]', 'group')
    if 'energyout [MeV]' in df:
        df = energy_to_group(df, 'energyout [MeV]', 'group out')
    return df


if __name__ == '__main__':
    # Build inputs and run OpenMC.
    build_inputs()
    executor = openmc.Executor()
    executor.run_simulation()

    # Read the tally data.
    sp = openmc.StatePoint('statepoint.15.h5')
    su = openmc.Summary('summary.h5')
    sp.link_with_summary(su)
    t = sp.tallies[1]
    df = t.get_pandas_dataframe()
    df = convert_e_to_g(df)

    # Extract the fast flux and shape it into y-x image format.
    df = df[df['group'] == 0]
    df = df[df['score'] == 'flux']
    flux = df['mean'].reshape((100, 100)).swapaxes(0, 1)

    # Make a blue-to-white colormap.
    cdict = {'red':   ((0.0, 0.0, 0.0),
                       (1.0, 1.0, 1.0)),
             'green': ((0.0, 0.0, 0.0),
                       (1.0, 1.0, 1.0)),
             'blue':  ((0.0, 1.0, 1.0),
                       (1.0, 1.0, 1.0))}
    cmap = matplotlib.colors.LinearSegmentedColormap('My CMap', cdict)

    # Plot the flux.
    plt.imshow(flux, cmap=cmap)
    plt.savefig('flux.png')
    plt.show()
