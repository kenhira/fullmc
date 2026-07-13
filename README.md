# FullMC
This is a fortran-based test code for forward and backward Monte Carlo radiative transfer simulation.

## Requirements
 - `gfortran` Fortran compiler
 - Python 3 with NumPy and Matplotlib

### Questions?
Contact Ken Hirata (@kenhira) for any questions.

---
### Internal
TODO

a. Incorporate photon path statistics within the trajectory recording

b. Cloud fraction implementation using cloudy/clear distinction

c. Backward MC sensor initial point update

d. Flexible dz

e. 1D vs 3D atmosphere implementation

Example

a. Output histogram of photon-averaged path lengths (source to sensor)

b. Cloudiness (input of cloud and clear states (n-th atmosphere))

c. Assume point spread function and detach (x, y) from domain (x, y, z)

d. Coarser dz toward TOA

e. Imhomogeneous atmosphere only for cloud regions (x, y, z)