#!/usr/bin/env python
# encoding: utf-8
r"""
2D shallow water: flow over a sill
==================================

Solve the 2D shallow water equations with horizontal zero width barrier and
variable bathymetry:

.. :math:
    h_t + (hu)_x + (hv)_y & = 0 \\
    (hu)_t + (hu^2 + \frac{1}{2}gh^2)_x + (huv)_y & = -g h b_x \\
    (hv)_t + (huv)_x + (hv^2 + \frac{1}{2}gh^2)_y & = -g h b_y.

The bathymetry is flat. The commented setting is the rotated equivalent example cf. diagonal barrier example
"""
from __future__ import absolute_import
from clawpack import riemann
from clawpack import pyclaw
from clawpack.riemann.shallow_roe_with_efix_2D_constants import depth, x_momentum, y_momentum, num_eqn
import numpy as np
from clawpack.pyclaw.plot import plot
import matplotlib.pyplot as plt

## barrier information:
bar_height = 1.64
bar_index_i = 50 #70
bar_loc =0.502 #0.7

def bathymetry(x,y):
    # r2 = 1.4285714285714285714285714285714*y-2.714285714285714
    # r2[:,bar_index_i+1] = r2[:,bar_index_i]
    # r2[:,bar_index_i-1] = r2[:,bar_index_i]
    return -2 #r2


def gauge_height(q,aux):
    h = q[0]
    return h

def setup(kernel_language='Fortran', solver_type='classic', use_petsc=False,
          outdir='./_output'):

    solver = pyclaw.ClawSolver2D(riemann.sw_aug_2D)
    solver.dimensional_split = True # No transverse solver available
    solver.order = 1

    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap
    solver.bc_lower[1] = pyclaw.BC.wall
    solver.bc_upper[1] = pyclaw.BC.extrap

    solver.aux_bc_lower[0] = pyclaw.BC.extrap
    solver.aux_bc_upper[0] = pyclaw.BC.extrap
    solver.aux_bc_lower[1] = pyclaw.BC.wall
    solver.aux_bc_upper[1] = pyclaw.BC.extrap

    my = 100 # 140
    mx = my

    x = pyclaw.Dimension(0.,1.,mx,name='x') # pyclaw.Dimension(-0.2,1.2,mx,name='x')
    y = pyclaw.Dimension(0.,1.,my,name='y') # pyclaw.Dimension(-0.2,1.2,my,name='y')
    domain = pyclaw.Domain([x,y])
    state = pyclaw.State(domain,num_eqn,num_aux=1) # the second and third aux var are small cell values

    X, Y = state.p_centers
    state.aux[0,:,:] = bathymetry(X,Y)
    state.q[depth,:,:] = 1.0 #-0.55-state.aux[0,:,:]
    state.q[depth,:20,:20 ] += 1.5
    # for i in range(int(mx/4),int(mx/2)):
    #     for j in range(int(mx/4)):
    #         if j+int(mx/4) <= i:
    #             state.q[depth,i,j] += 1.2
    # for i in range(int(mx/2),int(3*mx/4)):
    #     for j in range(int(mx/2)):
    #         if j <= -i+int(3*mx/4):
    #             state.q[depth,i,j] += 1.2
    state.q[x_momentum,:,:] = 0.
    state.q[y_momentum,:,:] = 0.

    state.problem_data['grav'] = 1.0
    state.problem_data['dry_tolerance'] = 1.e-14
    state.problem_data['sea_level'] = 0.
    state.problem_data['bar_index_i'] = bar_index_i
    state.problem_data['bar_ht'] = bar_height
    state.problem_data['bar_loc'] = bar_loc # be sure to reflect this location in the bar_index_i e.g. 0.602 would be i=30 since i=30 ~ xe(i)=0.6
    state.problem_data['orientation'] ='horz'
    state.problem_data['method'] = 'hbox'
    state.problem_data['alpha'] = 0.2


    claw = pyclaw.Controller()
    claw.tfinal = 1
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.num_output_times = 10
    claw.output_style = 1

    state.grid.add_gauges([(0.5,-0.025),(0.5,0.15),(0.5,0.325),(0.5,0.5),(0.5,0.675),(0.5,0.85),(0.5,1.025)])
    solver.compute_gauge_values = gauge_height
    state.keep_gauges = True
    claw.setplot = setplot
    claw.keep_copy = True

    claw.write_aux_always = True

    claw.run()
    plot(setplot=setplot,outdir='./_output',plotdir='./_plots',iplot=False,htmlplot=True)
    #return claw

def barrier_draw(current_data):
    x_1 = 1 #1.2
    x_0 = 0 #-0.2
    y_1 = 1 #1.2
    y_0 = 0 #-0.2
    bar_loc = 0.502 # 0.5
    axis = plt.gca()
    axis.plot([x_0,x_1],[bar_loc,bar_loc],'g',linewidth=1.5)
    return
def barrier_draw_1d(current_data):
    x_1 = 1.2
    x_0 = -0.2
    bar_loc = 0.502#0.5#
    bar_index = 50
    b = bathymetry(current_data.x,current_data.y)
    aux_wall = -2
    axis = plt.gca()
    axis.plot([bar_loc,bar_loc],[aux_wall,aux_wall+bar_height],'g',linewidth=1.5)
    # axis.plot(np.linspace(x_0,x_1,140),b[70,:],'k-')
    return

def surface_height(current_data):
    h = current_data.q[0,:,:]
    b = bathymetry(current_data.x,current_data.y)
    return h+b

def gauge_spots(current_data):
    gauge_points = [(0.5,-0.025),(0.5,0.15),(0.5,0.325),(0.5,0.5),(0.5,0.675),(0.5,0.85),(0.5,1.025)]
    axis = plt.gca()
    x_0 = -0.2 ; x_1 = 1.2
    bar_loc = 0.5
    axis.plot([x_0,x_1],[bar_loc,bar_loc],'g',linewidth=1.5)
    for i in range(len(gauge_points)):
        axis.plot(gauge_points[i][0],gauge_points[i][1],'k*')
        axis.annotate(str(i+1),(gauge_points[i][0],gauge_points[i][1]))
    return

def momentum_x(current_data):
    hu = current_data.q[1,:,5]
    return hu

def height_x(current_data):
    h = current_data.q[0,:,50]
    b = bathymetry(current_data.x,current_data.y)
    return h+b[:,50]

def height_y(current_data):
    h = current_data.q[0,70,:]
    b = bathymetry(current_data.x,current_data.y)
    return h + b[70,:]# -1 #

def setplot(plotdata):
    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for q[0]
    plotfigure = plotdata.new_plotfigure(name='Water height', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Water height with barrier'
    plotaxes.scaled = False
    # plotaxes.xlimits = [0,1]
    # plotaxes.ylimits = [0,1]

# cellby celll approach for afteraxes
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contourf')
    plotitem.plot_var =  surface_height
    plotitem.add_colorbar = True
    plotitem.contour_min = -0.6
    plotitem.contour_max = 0.6
    plotitem.contour_colors = 'blue'
    plotaxes.afteraxes = barrier_draw
    # plotaxes.xlimits = [-0.2,1.2]
    # plotaxes.ylimits = [-0.2,1.2]
    # plotaxes.afteraxes = gauge_spots

    plotfigure = plotdata.new_plotfigure(name="Momentum",figno=1)
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Momentum in x"
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = momentum_x


    plotfigure = plotdata.new_plotfigure(name="Height",figno=2)
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Height in y"
    plotaxes.xlimits = [0,1] #[-0.2,1.2]
    plotaxes.ylimits = [-2.1, 1.2] #[-3.2,1.1]
    plotaxes.afteraxes = barrier_draw_1d

    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = height_y

    return plotdata

if __name__=="__main__":
    setup()
