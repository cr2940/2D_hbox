# 2D hbox codes for barrier simulations

   ![image](https://user-images.githubusercontent.com/36740525/115186258-5872b200-a11c-11eb-86ab-cba8a3acbcd8.png)
(Picture of 20 degree barrier blocking a wave)

Solving 2D SWE in PyClaw with zero width barriers: flat, diagonal, slanted (~20 degree) barriers

The setup files are named `barrier_flat.py` , `barrier_diagonal.py` , and `barrier_slanted.py`.

The file `barrier_flat.py` runs an example of horizontal barrier just off grid edge by distance of 0.2 * dx with a dam break producing an oblique wave.

The file `barrier_diagonal.py` runs an example of diagonal barrier (45 degree to the grid) with a dam break producing a radial wave headed towards the barrier.
In `barrier_diagonal.py` there are two options (one commented out) for bathymetric setting, one being flat and the other being a sloping beach. Note that around the barrier in the sloping beach example the bathymetry is flat. 

The file `barrier_slanted.py` runs an example of slanted barrier as an example of arbitrarily angled barrier (~20 degree to the grid) with a dam break producing a planar wave headed towards the barrier. The bathymetry is flat. To run `barrier_slanted.py` you need the data (in .txt) in the `data` directory, to be in the same directory as this setup file.

The solver and object/class files required to run the examples are in "PyClaw_hbox2D" [repo](https://github.com/cr2940/PyClaw_hbox2D).


