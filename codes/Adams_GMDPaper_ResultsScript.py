"""
overland_flow_driver.py

OverlandFlow component example, initializing a 30 m square watershed from an
ESRI ASCII file, simulating one of the rainfall events from Adams et al.,
in prep for Geoscientific Model Development.

Written by Jordan Adams, August 2016

"""

from landlab.components import OverlandFlow, DetachmentLtdErosion
from landlab.io import read_esri_ascii
from landlab.grid.raster_mappers import map_max_of_inlinks_to_node
import numpy as np
from matplotlib import pyplot as plt

### For each storm event described in Adams et al., in prep, set both the
### basin flag and the storm flag for the event.
basin_flag = 'Square' # 'Long'
storm_flag = 'Base' # 'HigherIntensity' # 'LongerDuration'

### If the basin flag matches one of the two selec basins, the filename and
### the outlet link to sample discharge values from for plotting.
if basin_flag == 'Square':
    watershed_dem = 'Square_TestBasin.asc'
    link_to_sample = 299

else:
    watershed_dem = 'Long_TestBasin.asc'
    link_to_sample = 149

### Reading in the ESRI ASCII DEM given the filename from above, and setting
### the topographic__elevation field.
(rmg, z) = read_esri_ascii(watershed_dem, name='topographic__elevation')

### Setting the boundary conditions: all NODATA nodes in the DEM are closed
### boundaries, and the outlet is set to an open boundary.
rmg.set_watershed_boundary_condition(z)

### Initializing the OverlandFlow and DetachmentLtdErosion components,
### including stability criteria and parameters that can be tuned.
of = OverlandFlow(rmg, steep_slopes=True)
dle = DetachmentLtdErosion(rmg, K_sp = 1.259162261 * (10**-7))

### Setting precipitation parameters based on the storm__flag named above.
### The three storms are "Base" (5 mm/hr for 2 hr), "HigherIntensity (10 mm/hr
### for 2 hr) and "LongerDuration" (5 mm/hr for 4 hr).
if storm_flag == 'Base':
    starting_precip_mmhr = 5.0
    starting_precip_ms = starting_precip_mmhr * (2.77778 * 10 ** -7)
    storm_duration = 7200.
elif storm_flag == 'HigherIntensity':
    starting_precip_mmhr = 10.0
    starting_precip_ms = starting_precip_mmhr * (2.77778 * 10 ** -7)
    storm_duration = 7200.
elif storm_flag == 'LongerDuration':
    starting_precip_mmhr = 5.0
    starting_precip_ms = starting_precip_mmhr * (2.77778 * 10 ** -7)
    storm_duration = 14400.

uplift_rate = 3.170979 * (10**-10) # m/s

## Model Run Time in seconds
elapsed_time = 1.0
model_run_time = 43200.0

## Lists for saving data
discharge_at_outlet = []
hydrograph_time = []

## Setting initial fields...
rmg['node']['surface_water__discharge'] = np.zeros(rmg.number_of_nodes)
rmg['node']['water_surface__slope'] = np.zeros(rmg.number_of_nodes)

## Running the overland flow component.
while elapsed_time < model_run_time:

    # Setting the adaptive time step
    of.dt = of.calc_time_step()

    ## The storm starts when the model starts. While the elapsed time is less
    ## than the storm duration, we add water to the system as rainfall.
    if elapsed_time < (storm_duration):

        of.rainfall_intensity =  starting_precip_ms

    ## Then the elapsed time exceeds the storm duration, rainfall ceases.
    else:

        of.rainfall_intensity = 0.0

    ## Generating overland flow based on the deAlmeida solution.
    of.overland_flow()

    ## Mapping water discharge from links (m^2/s) to nodes (m^3/s) for use
    ## in the DetachmentLtdErosion component.
    rmg['node']['surface_water__discharge'] = (map_max_of_inlinks_to_node(
                                                rmg, np.abs(of.q)) * rmg.dx)

    ## Calculating water surface slope from the OverlandFlow component.
    rmg['node']['water_surface__slope'] = ((of.slope[rmg.links_at_node] *
                                    rmg.active_link_dirs_at_node).max(axis=1))

    ## Eroding topographic__elevation using DetachmentLtdErosion component.
    dle.erode(of.dt, slope='water_surface__slope')

    ## Updating topographic__elevation after surface was eroded in
    ## DetachmentLtdErosion component.
    rmg['node']['topographic__elevation'] += uplift_rate * of.dt

    ## Append time and discharge to their lists to save data and for plotting.
    hydrograph_time.append(elapsed_time)
    discharge_at_outlet.append(np.abs(of.q[link_to_sample]) * rmg.dx)

    print(elapsed_time)
    ## Updating elapsed_time
    elapsed_time += of.dt

## Plotting the hydrograph at the outlet
plt.figure(1)
plt.plot(hydrograph_time, discharge_at_outlet, 'k-')
plt.ylabel('Discharge (cms)')
plt.xlabel('Time (hr)')
plt.title('Hydrograph at Basin Outlet')
