# Theoretical DDM generation with Zavorotny-Voronovich model for ocean reflected signals

Generation of delay-Doppler maps (DDM) or correlation waveforms (WF) for ocean reflected GPS C/A signals according to the Zavorotny and Voronovich model from V. U. Zavorotny and A. G. Voronovich, "Scattering of GPS signals from the ocean with wind remote sensing application," in IEEE Transactions on Geoscience and Remote Sensing, vol. 38, no. 2, pp. 951-964, March 2000, doi: 10.1109/36.841977.

## Description

These scripts solve the integral expresion of the mean DDM/WF obtained according to the Zavorotny-Voronovich model. This implementation has the bistatic geometry as input as well as the ocean state.

### Executing program

Define the simulation scenario in the config.m script. The input parameters are:

* OCEAN MODEL: wind speed, direction and ocean wave age. Optionally it can be set to simulate the reflection over an ideally flat surface.
* GEOMETRY: GPS satellite elevation angle, moving direction of GPS and LEO satellites, LEO altitude.
* RECEIVER: receiver front-end noise figure
* DDM: define the Delay and Doppler axis for the DDM, or only the delay axis for the WF case.

After configuration, run main.m
