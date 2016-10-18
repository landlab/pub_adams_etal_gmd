"""
Adams_GMDPaper_SampleScript.py

OverlandFlow component example, initializing a 30 m square watershed from an
ESRI ASCII file, simulating the base rainfall event from Adams et al.,
in prep for Geoscientific Model Development.

Written by Jordan Adams, October 2016

"""

## Import necessary libraries from Landlab
from landlab.components import OverlandFlow, DetachmentLtdErosion, SinkFiller
from landlab.io import read_esri_ascii

## Read in the DEM and create a RasterModelGrid instance
(rmg, z) = read_esri_ascii('SpringCreek_DEM.asc', name='topographic__elevation')

## Set the boundary conditions
rmg.set_watershed_boundary_condition(z)

## Set rainfall parameters
starting_precip_mmhr = 5.0
starting_precip_ms = starting_precip_mmhr * (2.77778 * 10 ** -7)

## Initialize all Landlab components
of = OverlandFlow(rmg, steep_slopes=True, rainfall_intensity = starting_precip_ms)
dle = DetachmentLtdErosion(rmg, K_sp = 1.259162261 * (10**-7))
sf = SinkFiller(rmg, routing='D4', apply_slope=True, fill_slope=1.e-5)

## Fill pits in D4 scheme. This is computationally intensive and is 
## commented out by default - the SpringCreek_DEM is already processed.
#sf.fill_pits()


## Time loop parameters
elapsed_time = 1.0
model_run_time = 36000.0

## Loop through time - the storm parameters are continuous (line 27)
while elapsed_time < model_run_time:

    # Calculate stable time step.
    of.dt = of.calc_time_step()

    # Generate flow
    of.overland_flow()

    # Map discharge from links to nodes.
    rmg['node']['surface_water__discharge'] = of.discharge_mapper(of.q, 
                                                        convert_to_volume=True)
    
    # Map slope at steepest descent from links to nodes.
    rmg['node']['water_surface__slope'] = (of.water_surface_slope[
                rmg.links_at_node] * rmg.active_link_dirs_at_node).max(axis=1)

    # Erode
    dle.erode(of.dt, discharge_cms = 'surface_water__discharge', 
                                                  slope='water_surface__slope')

    elapsed_time += of.dt