# 2D_hbox
Solving 2D SWE in PyClaw with zero width barriers: flat and diagonal barriers

The setup files are named `barrier_flat.py` and `barrier_diagonal.py`.

The file `barrier_flat.py` runs an example of horizontal barrier just off grid edge by distance of 0.2 * dx with a dam break producing an oblique wave.

The file `barrier_diagonal.py` runs an example of diagonal barrier (45 degree to the grid) with a dam break producing a radial wave headed towards the barrier.
In `barrier_diagonal.py` there are two options (one commented out) for bathymetric setting, one being flat and the other being a sloping beach. Note that around the barrier in the sloping beach example the bathymetry is flat. 

The solver and object/class files required to run the examples are in "PyClaw_hbox2D" repo.


