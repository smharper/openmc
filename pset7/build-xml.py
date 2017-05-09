from __future__ import print_function

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import openmc


def build_inputs(U235_content, particles, batches, sp_every_batch=False,
                 output_to_tmp=False):
    # Materials.
    fuel = openmc.Material()
    fuel.set_density('g/cc', 1.0)
    fuel.add_nuclide('H1', 1.0)
    fuel.add_nuclide('U235', U235_content)

    materials = openmc.Materials([fuel])
    materials.export_to_xml()

    # Geometry.
    x0 = openmc.XPlane(x0=-150, boundary_type='reflective')
    x1 = openmc.XPlane(x0=150, boundary_type='reflective')
    y0 = openmc.YPlane(y0=-150, boundary_type='reflective')
    y1 = openmc.YPlane(y0=150, boundary_type='reflective')
    z0 = openmc.ZPlane(z0=-150, boundary_type='reflective')
    z1 = openmc.ZPlane(z0=10, boundary_type='reflective')

    cell = openmc.Cell()
    cell.region = +x0 & -x1 & +y0 & -y1 & +z0 & -z1
    cell.fill = fuel

    root = openmc.Universe()
    root.add_cell(cell)

    geometry = openmc.Geometry(root)
    geometry.export_to_xml()

    # Settings.
    settings = openmc.Settings()
    settings.batches = batches
    settings.inactive = 0
    settings.particles = particles
    settings.source = openmc.Source(space=openmc.stats.Point())
    if sp_every_batch:
        settings.statepoint = {'batches':[b+1 for b in range(batches)]}
    if output_to_tmp:
        settings.output = {'path':'/tmp/sp_scratch'}
    settings.export_to_xml()

    # Tallies.
    tallies = openmc.Tallies()
    #t = openmc.Tally(tally_id=1)
    #t.filters = [openmc.EnergyFilter(np.logspace(-5, np.log10(20e6), 200))]
    #t.scores = ['flux']
    #tallies = openmc.Tallies([t])
    t = openmc.Tally(tally_id=1)
    t.scores = ['absorption-age', 'absorption']
    tallies.append(t)

    t = openmc.Tally(tally_id=2)
    t.scores = ['absorption']
    t.filters = [openmc.TimeFilter(np.linspace(0, 1e3, int(1e3)))]
    tallies.append(t)

    t = openmc.Tally(tally_id=3)
    t.scores = ['absorption']
    t.filters = [openmc.TimeFilter(np.linspace(0, 1500*1.31e-5, int(1e3)))]
    tallies.append(t)
    tallies.export_to_xml()


if __name__ == '__main__':
    # This quantity of U-235 makes the reactor within ~10 pcm of critical.  If
    # you want the batch-by-batch plots, set sp_every_batch to True to output a
    # statepoint at every batch.  I would recommend also setting output_to_tmp
    # to True which puts all the statepoints in the /tmp/sp_scratch directory
    # (if you get an HDF5 error, mkdir that directory) so that the ~1 GB of
    # files are deleted on restart.
    U235 = 4.620e-4
    particles = int(10e3)
    batches = 1000
    build_inputs(U235, particles, batches, sp_every_batch=False,
                 output_to_tmp=False)
