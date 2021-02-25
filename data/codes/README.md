This directory contains the codes that generates the data. The codes create data for *line segment barrier only*. Due to floating point errors, code does *not* work for some coordinate points as start and end point barrier. Some exploration with different points is required to find working cases.

`set_aux_SRD_hbox.f90` is the setup file where you prescribe the barrier (e.g. starting point, end point) that then calls the codes to compute the required data.

`aux_module_SRD2.f90` is the module that computes the geometry of small cells caused by the barrier.

`hbox_module.f90` is the module that computes the h-box geometry (e.g. required to compute h-box averages along the barrier).

`aux_hbox2.f90` is the module that provides auxiliary functions to do computation in h-box module.

`geom_calc.sh` is the compilation script that provides the executable to run `set_aux_SRD_hbox.f90`.
