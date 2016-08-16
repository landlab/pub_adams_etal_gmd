"""
overland_flow_driver.py

OverlandFlow component example, initializing a 30 m square watershed from an
ESRI ASCII file, simulating a 5 mm/hr rainfall intensity over 2 hours, the
standard storm example from Adams et al., in prep for Geoscientific Model Development.

Written by Jordan Adams, August 2016


## Step 1
from landlab.components import OverlandFlow, SinkFiller # SinkFiller is optional, Step #4
from landlab.io import read_esri_ascii
import numpy as np

from landlab.plot import imshow_grid  # plotter functions are optional
from matplotlib import pyplot as plt

## Step 2: Reading in a watershed DEM
watershed_dem = 'Square_TestBasin.asc'
(rmg, z) = read_esri_ascii(watershed_dem, name='topographic__elevation')

## --OR-- Step 2: Setting up a generic RasterModelGrid
#rmg = RasterModelGrid((number_of_node_rows, number_of_node_columns), dx)
#z = user_defined_elevation_data        # length of number_of_nodes
#rmg[‘node’][‘topographic__elevation’] = z

## Step 3
rmg.set_watershed_boundary_condition(z)

## Step 4: Preprocessing the DEM (Optional)
# sf = SinkFiller(rmg, routing='D4', apply_slope=True, fill_slope=1.e-5)
# sf.fill_pits()

## Step 5
of = OverlandFlow(rmg, steep_slopes=True)

## Step 6
elapsed_time = 1.0
model_run_time = 7200.

storm_duration = 7200.0
rainfall_mmhr = 5.

## For plotting (Optional)
#hydrograph_time = []
#discharge_at_outlet = []
#outlet_link = 299

## Step 7
while elapsed_time < model_run_time:

    # Adaptive time step
    of.dt = of.calc_time_step()

    if elapsed_time < (storm_duration):
        of.rainfall_intensity =  rainfall_mmhr * (2.777778 * 10**-7)

    else:
        of.rainfall_intensity = 0.0

    of.overland_flow()

## For plotting a hydrograph (Optional)
#    hydrograph_time.append(elapsed_time / 3600.)
#    discharge_at_outlet.append(np.abs(of.q[outlet_link]) * rmg.dx)

    elapsed_time += of.dt

## Plotting a hydrograph
#plt.plot(hydrograph_time, discharge_at_outlet)
#plt.ylabel('Time (hr)')
#plt.xlabel('Discharge, (cms)')
#plt.title('Outlet Hydrograph, Rainfall: 5 mm/hr in 2 hr')

## Plotting water depth maps
imshow_grid(rmg, 'water__depth', plot_name='Water depth at time = 2 hr',
        var_name='Water Depth', var_units='m', grid_units=('m', 'm'),
        cmap='Blues')
