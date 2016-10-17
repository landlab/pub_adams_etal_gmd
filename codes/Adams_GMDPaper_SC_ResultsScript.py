"""
Adams_GMDPaper_Synthetic_ResultsScript.py

OverlandFlow component example, initializing a 30 m square watershed from an
ESRI ASCII file, simulating one of the rainfall events from Adams et al.,
in prep for Geoscientific Model Development.

Written by Jordan Adams, October 2016

"""

from landlab.components import OverlandFlow
from landlab.io import read_esri_ascii
import numpy as np
from matplotlib import pyplot as plt

### Setting the DEM filename and the outlet link to sample discharge values
### for plotting.
watershed_dem = 'SpringCreek_DEM.asc'
link_to_sample = 110756

### Reading in the ESRI ASCII DEM given the filename from above, and setting
### the topographic__elevation field.
(rmg, z) = read_esri_ascii(watershed_dem, name='topographic__elevation')

### Setting the boundary conditions: all NODATA nodes in the DEM are closed
### boundaries, and the outlet is set to an open boundary.
rmg.set_watershed_boundary_condition(z)

### Initializing the OverlandFlow component,
### including stability criteria and parameters that can be tuned.
of = OverlandFlow(rmg, steep_slopes=True, alpha=0.05, mannings_n=0.055)

### Setting precipitation parameters based on the storm__flag named above.
### The storm modeled here is the "Base" (5 mm/hr for 2 hr),
### from Adams et al., in prep.
starting_precip_mmhr = 5.0
starting_precip_ms = starting_precip_mmhr * (2.77778 * 10 ** -7)
storm_duration = 7200.

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
