# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 15:34:44 2018

@author: L1817
"""
import numpy as np
import fipy as fp
import matplotlib.pyplot as plt
import copy
import os
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable # for plots


import hydro_utils, utilities

"""
   SET FIPY SOLVER
"""
fp.solvers.DefaultSolver = fp.solvers.LinearLUSolver

# Define coefficients for WTD-CO2 emissions relationship (Equation 8 and 9 in the manuscript)
ALPHA = 74.11  # Mg/ha/yr/m
BETA = 29.34   # Mg/ha/yr


def big_4_raster_plot(title, raster1, raster2, raster3, raster4, output_folder): 
    """
    Plots four raster arrays in a 2x2 grid and saves the plot as a PNG file.

    Args:
        title (str): Title for the entire plot.
        raster1 (numpy.ndarray): Raster array for the top left plot (D).
        raster2 (numpy.ndarray): Raster array for the top right plot (canal water level).
        raster3 (numpy.ndarray): Raster array for the bottom left plot (DEM).
        raster4 (numpy.ndarray): Raster array for the bottom right plot (elevation - phi).
        output_folder (str): Path to the folder where the plot should be saved.
    """
    fig1, axes1 = plt.subplots(nrows=2, ncols=2, figsize=(16, 12), dpi=80)
    fig1.suptitle(title)

    a = axes1[1, 0].imshow(raster1, cmap='pink', interpolation='nearest')
    ax1_divider = make_axes_locatable(axes1[1, 0])
    cax1 = ax1_divider.append_axes('right', size='7%', pad='2%')
    plt.colorbar(a, cax=cax1)
    axes1[1, 0].set(title="D")

    b = axes1[0, 1].imshow(raster2, cmap='viridis')
    ax2_divider = make_axes_locatable(axes1[0, 1])
    cax2 = ax2_divider.append_axes('right', size='7%', pad='2%')
    plt.colorbar(b, cax=cax2)
    axes1[0, 1].set(title="canal water level")

    c = axes1[0, 0].imshow(raster3, cmap='viridis')
    ax3_divider = make_axes_locatable(axes1[0, 0])
    cax3 = ax3_divider.append_axes('right', size='7%', pad='2%')
    plt.colorbar(c, cax=cax3)
    axes1[0, 0].set(title="DEM")

    d = axes1[1, 1].imshow(raster4, cmap='pink')
    ax4_divider = make_axes_locatable(axes1[1, 1])
    cax4 = ax4_divider.append_axes('right', size='7%', pad='2%')
    plt.colorbar(d, cax=cax4)
    axes1[1, 1].set(title="elevation - phi")
    
    # Save the plot as PNG
    plot_filename = f"{title.replace(' ', '_')}.png" # Replace spaces with underscores in filename
    plot_filepath = os.path.join(output_folder, plot_filename) 
    plt.savefig(plot_filepath)
    plt.close(fig1) # Close the figure to release resources


def plot_line_of_peat(raster, y_value, title, nx, ny, label, color, 
    linewidth=1., output_folder=None):
    """
    Plots a cross-section of a raster array at a specified y-value.

    Args:
        raster (numpy.ndarray): The raster array to plot.
        y_value (int): The y-value (row) at which to take the cross-section.
        title (str): Title for the plot.
        nx (int): Number of columns in the raster.
        ny (int): Number of rows in the raster.
        label (str): Label for the plot line.
        color (str): Color of the plot line.
        linewidth (float, optional): Width of the plot line. Defaults to 1.0.
        output_folder (str, optional): Path to the folder where the plot should be saved. Defaults to None (no saving).
    """
    plt.figure(10)
    plt.title(title)
    plt.plot(raster[y_value, :], label=label, color=color, linewidth=linewidth)
    plt.legend(fontsize='x-large')
    
    if output_folder is not None:
        plot_filename = f"{title.replace(' ', '_')}.png"
        plot_filepath = os.path.join(output_folder, plot_filename)
        plt.savefig(plot_filepath)
        plt.close() # Close the figure 

def plot_hydro_params(httd, avg_wt, output_folder):
    """
    Plots hydraulic parameters and average water table depth over time, saving the plot to a PNG file.

    Args:
        httd (dict): Dictionary containing hydraulic properties.
        avg_wt (list): Average water table depth over time.
        output_folder (str): Path to the folder where the plot should be saved.
    """
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 9), dpi=80)
    x = np.arange(-21, 1, 0.1)

    axes[0, 0].plot(httd[1]['hToTra'](x), x)
    axes[0, 0].set(title='hToTra', ylabel='depth')

    axes[0, 1].plot(httd[1]['C'](x), x)
    axes[0, 1].set(title='C')

    axes[1, 0].plot()  # This subplot seems to be intentionally left empty
    axes[1, 0].set(title="Nothing") 

    axes[1, 1].plot(avg_wt)
    axes[1, 1].set(title="avg_wt_over_time")

    # Save the plot as a PNG file
    plot_filename = "hydrological_parameters.png"
    plot_filepath = os.path.join(output_folder, plot_filename)
    plt.savefig(plot_filepath)
    plt.close(fig) # Close the figure
    
def get_rasters(ny, nx, peat_type_mask, gwt, httd, tra_to_cut, cmask, drmask_not, ele, wt_canal_arr, dr, catchment_mask):
    """
    Compute and return the four raster arrays used in the big_4_raster_plot function.

    Args:
        ny (int): Number of rows in the raster arrays.
        nx (int): Number of columns in the raster arrays.
        peat_type_mask (numpy.ndarray): Mask array indicating the peat type.
        gwt (numpy.ndarray): Groundwater table elevation array.
        httd (dict): Dictionary containing peat hydraulic properties.
        tra_to_cut (numpy.ndarray): Transmissivity values to cut.
        cmask (fipy.CellVariable): Cell variable mask array.
        drmask_not (fipy.CellVariable): Cell variable mask array for non-drain cells.
        ele (numpy.ndarray): Elevation array.
        wt_canal_arr (numpy.ndarray): Canal water table elevation array.
        dr (numpy.ndarray): Drain mask array.
        catchment_mask (numpy.ndarray): Catchment mask array.

    Returns:
        tuple: A tuple containing the four raster arrays (rast_D, rast_cwl, rast_dem, rast_elev_phi).
    """
    rast_D = ((hydro_utils.peat_map_h_to_tra(soil_type_mask=peat_type_mask, gwt=gwt, h_to_tra_and_C_dict=httd) - tra_to_cut) * cmask.value * drmask_not.value).reshape(ny, nx)
    rast_cwl = (ele.reshape(ny, nx) - wt_canal_arr) * dr * catchment_mask
    rast_dem = ele.reshape(ny, nx)
    rast_elev_phi = gwt.reshape(ny, nx)
    return rast_D, rast_cwl, rast_dem, rast_elev_phi



def hydrology(solve_mode, nx, ny, dx, dy, days, ele, phi_initial, catchment_mask, wt_canal_arr, boundary_arr,
              peat_type_mask, httd, tra_to_cut, sto_to_cut, track_WT_drained_area, track_WT_notdrained_area,
              diri_bc=0.0, neumann_bc=None, plotOpt=False, remove_ponding_water=True, P=0.0, ET=0.0, dt=1.0, 
              n_day_avg=3, output_folder=None, y_value=270):
    """
    Simulates peatland hydrology using the Boussinesq equation.

    Args:
        solve_mode (str): 'steadystate' or 'transient'. 'transient' was used in the study.
        nx (int): Number of grid cells in the x-direction.
        ny (int): Number of grid cells in the y-direction.
        dx (float): Grid cell size in the x-direction (meters).
        dy (float): Grid cell size in the y-direction (meters).
        days (int): Number of days to simulate.
        ele (numpy.ndarray): Elevation array (meters above sea level).
        phi_initial (numpy.ndarray): Initial hydraulic head array (meters above sea level).
        catchment_mask (numpy.ndarray): Boolean mask indicating the catchment area.
        wt_canal_arr (numpy.ndarray): Array of canal water table elevations (meters above sea level).
        boundary_arr (numpy.ndarray): Array of Dirichlet boundary conditions.
        peat_type_mask (numpy.ndarray): Array indicating peat type for each grid cell.
        httd (dict): Dictionary containing hydraulic properties for each peat type.
        tra_to_cut (numpy.ndarray): Transmissivity to be subtracted based on impermeable bottom depth.
        sto_to_cut (numpy.ndarray): Storage coefficient to be subtracted based on impermeable bottom depth.
        track_WT_drained_area (tuple): Coordinates of a point near a drained area for WTD tracking.
        track_WT_notdrained_area (tuple): Coordinates of a point away from drained areas for WTD tracking.
        diri_bc (float, optional): Dirichlet boundary condition value. Defaults to 0.0.
        neumann_bc (float, optional): Neumann boundary condition value. Defaults to None (not used in the study).
        plotOpt (bool, optional): Whether to generate plots. Defaults to False.
        remove_ponding_water (bool, optional): Whether to remove water above the surface. Defaults to True.
        P (float or numpy.ndarray, optional): Precipitation (mm/day). Can be a single value or an array for each day. 
            Defaults to 0.0.
        ET (float or numpy.ndarray, optional): Evapotranspiration (mm/day). Can be a single value or an array for each 
            day. Defaults to 0.0.
        dt (float, optional): Time step (days). Defaults to 1.0.
        n_day_avg (int, optional): Number of days to average WTD over. Defaults to 3.

    Returns:
        tuple: A tuple containing:
            - cumulative_Vdp (float): Cumulative dry peat volume over the simulation (m^3).
            - dry_peat_volume (float): Dry peat volume at the end of the simulation (m^3).
            - avg_wt_over_time (float): Average WTD over the entire simulation (m).
            - wt_track_drained (list): WTD values at the drained location over time (m).
            - wt_track_notdrained (list): WTD values at the undrained location over time (m).
            - timestep_data (list): List of lists containing timestep data.
            - rast_... (tuples): Raster arrays before and after the simulation.
            - avg_wtd_nday (list): N-day average WTD values (m).
    """

    ele[~catchment_mask] = 0.
    ele = ele.flatten()
    phi_initial = (phi_initial + 0.0 * np.zeros((ny, nx))) * catchment_mask
    #    phi_initial = phi_initial * catchment_mask
    phi_initial = phi_initial.flatten()

    if len(ele) != nx * ny or len(phi_initial) != nx * ny:
        raise ValueError("ele or Hinitial are not of dim nx*ny")

    mesh = fp.Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)
    phi = fp.CellVariable(name='computed H', mesh=mesh, value=phi_initial, hasOld=True)  # response variable H in meters above reference level

    if diri_bc is not None and neumann_bc is None:
        phi.constrain(diri_bc, mesh.exteriorFaces)

    elif diri_bc is None and neumann_bc is not None:
        phi.faceGrad.constrain(neumann_bc * mesh.faceNormals, where=mesh.exteriorFaces)

    else:
        raise ValueError(
            "Cannot apply Dirichlet and Neumann boundary values at the same time. Contradictory values.")

    # *******omit areas outside the catchment. c is input
    cmask = fp.CellVariable(mesh=mesh, value=np.ravel(catchment_mask))
    cmask_not = fp.CellVariable(mesh=mesh, value=np.array(~ cmask.value, dtype=int))

    # *** drain mask or canal mask
    dr = np.array(wt_canal_arr, dtype=bool)
    dr[np.array(wt_canal_arr, dtype=bool) * np.array(boundary_arr, dtype=bool)] = False  # Pixels cannot be canals and boundaries at the same time. Everytime a conflict appears, boundaries win. This overwrites any canal water level info if the canal is in the boundary.
    drmask = fp.CellVariable(mesh=mesh, value=np.ravel(dr))
    drmask_not = fp.CellVariable(mesh=mesh, value=np.array(
        ~ drmask.value, dtype=int))  # Complementary of the drains mask, but with ints {0,1}

    # mask away unnecesary stuff
    #    phi.setValue(np.ravel(H)*cmask.value)
    #    ele = ele * cmask.value

    source = fp.CellVariable(mesh=mesh, value=0.)  # cell variable for source/sink
    #    CC=fp.CellVariable(mesh=mesh, value=C(phi.value-ele))                   # differential water capacity

    def D_value(phi, ele, tra_to_cut, cmask, drmask_not):
        # Some inputs are in fipy CellVariable type
        gwt = phi.value * cmask.value - ele

        d = hydro_utils.peat_map_h_to_tra(
            soil_type_mask=peat_type_mask, gwt=gwt, h_to_tra_and_C_dict=httd) - tra_to_cut

        # d <0 means tra_to_cut is greater than the other transmissivity, which in turn means that
        # phi is below the impermeable bottom. We allow phi to have those values, but
        # the transmissivity is in those points is equal to zero (as if phi was exactly at the impermeable bottom).
        d[d < 0] = 1e-3  # Small non-zero value not to wreck the computation

        dcell = fp.CellVariable(
            mesh=mesh, value=d)  # diffusion coefficient, transmissivity. As a cell variable.
        dface = fp.FaceVariable(mesh=mesh, value=dcell.arithmeticFaceValue.value)  # The correct Face variable.

        return dface.value

    def C_value(phi, ele, sto_to_cut, cmask, drmask_not):
        # Some inputs are in fipy CellVariable type
        gwt = phi.value * cmask.value - ele

        c = hydro_utils.peat_map_h_to_sto(
            soil_type_mask=peat_type_mask, gwt=gwt, h_to_tra_and_C_dict=httd) - sto_to_cut
        c[c < 0] = 1e-3  # Same reasons as for D

        ccell = fp.CellVariable(
            mesh=mesh, value=c)  # diffusion coefficient, transmissivity. As a cell variable.
        return ccell.value

    D = fp.FaceVariable(mesh=mesh, value=D_value(
        phi, ele, tra_to_cut, cmask, drmask_not))  # The correct Face variable.
    C = fp.CellVariable(mesh=mesh, value=C_value(
        phi, ele, sto_to_cut, cmask, drmask_not))  # differential water capacity

    largeValue = 1e20  # value needed in implicit source term to apply internal boundaries

    # Get raster arrays before the computation
    rast_D_before, rast_cwl_before, rast_dem_before, rast_elev_phi_before = get_rasters(
        ny, nx, peat_type_mask, phi.value - ele, httd, tra_to_cut, cmask, drmask_not, ele, wt_canal_arr, dr, catchment_mask)

    if plotOpt:

        # Plot raster arrays before the computation
        big_4_raster_plot(title='Before the computation',
                          raster1=rast_D_before,
                          raster2=rast_cwl_before,
                          raster3=rast_dem_before,
                          raster4=rast_elev_phi_before,
                          output_folder=output_folder)
        plt.show()

        #        print "first cross-section plot"
        #        ele_with_can = copy.copy(ele).reshape(ny,nx)
        #        ele_with_can = ele_with_can * catchment_mask
        #        ele_with_can[wt_canal_arr > 0] = wt_canal_arr[wt_canal_arr > 0]
        #        plot_line_of_peat(ele_with_can, y_value=y_value, title="cross-section", color='green', nx=nx, ny=ny, label="ele")

    # ********************************** PDE, STEADY STATE **********************************
    if solve_mode == 'steadystate':
        if diri_bc is not None:
            #        diri_boundary = fp.CellVariable(mesh=mesh, value= np.ravel(diri_boundary_value(boundary_mask, ele2d, diri_bc)))

            eq = 0. == (fp.DiffusionTerm(coeff=D)
                        + source * cmask * drmask_not
                        - fp.ImplicitSourceTerm(
                            cmask_not * largeValue) + cmask_not * largeValue * np.ravel(boundary_arr)
                        - fp.ImplicitSourceTerm(
                            drmask * largeValue) + drmask * largeValue * (np.ravel(wt_canal_arr))
                        #                    - fp.ImplicitSourceTerm(bmask_not*largeValue) + bmask_not*largeValue*(boundary_arr)
                        )

        elif neumann_bc is not None:
            raise NotImplementedError("Neumann BC not implemented yet!")
            cmask_face = fp.FaceVariable(
                mesh=mesh, value=np.array(cmask.arithmeticFaceValue.value, dtype=bool))
            D[cmask_face.value] = 0.
            eq = 0. == (fp.DiffusionTerm(coeff=D) + source * cmask * drmask_not
                        - fp.ImplicitSourceTerm(
                            cmask_not * largeValue) + cmask_not * largeValue * (diri_bc)
                        - fp.ImplicitSourceTerm(
                            drmask * largeValue) + drmask * largeValue * (np.ravel(wt_canal_arr))
                        #                + fp.DiffusionTerm(coeff=largeValue * bmask_face)
                        #                - fp.ImplicitSourceTerm((bmask_face * largeValue *neumann_bc * mesh.faceNormals).divergence)
                        )

    elif solve_mode == 'transient':
        if diri_bc is not None:
            #        diri_boundary = fp.CellVariable(mesh=mesh, value= np.ravel(diri_boundary_value(boundary_mask, ele2d, diri_bc)))

            eq = fp.TransientTerm(coeff=C) == (fp.DiffusionTerm(coeff=D)
                                              + source * cmask * drmask_not
                                              - fp.ImplicitSourceTerm(
                                                  cmask_not * largeValue) + cmask_not * largeValue * np.ravel(boundary_arr)
                                              - fp.ImplicitSourceTerm(
                                                  drmask * largeValue) + drmask * largeValue * (np.ravel(wt_canal_arr))
                                              #                        - fp.ImplicitSourceTerm(bmask_not*largeValue) + bmask_not*largeValue*(boundary_arr)
                                              )
        elif neumann_bc is not None:
            raise NotImplementedError("Neumann BC not implemented yet!")

    # ********************************************************

    max_sweeps = 10  # inner loop.

    avg_wt = []
    wt_track_drained = []
    wt_track_notdrained = []
    # Create a list to store timestep data
    timestep_data = []
    # List to store n-day average WTD
    avg_wtd_nday = [] 

    cumulative_Vdp = 0.
    dry_peat_volume = 0.0

    # ********Finite volume computation******************
    for d in range(days):

        if isinstance(P, np.ndarray):  # assume it is a numpy array
            source.setValue(
                (P[d] - ET[d]) * .001 * np.ones(ny * nx))  # source/sink, in mm/day. The factor of 10^-3 takes into account that there are 100 x 100 m^2 in one pixel
            print("(d,P) = ", (d, (P[d] - ET[d]) * 1.))
        else:
            source.setValue(
                (P - ET) * 10. * np.ones(ny * nx))
            print("(d,P) = ", (d, (P - ET) * 10.))

        if plotOpt and d != 0:
            # print "one more cross-section plot"
            plot_line_of_peat(phi.value.reshape(
                ny, nx), y_value=y_value, title="cross-section", color='cornflowerblue', nx=nx, ny=ny, label=d, 
                             output_folder=output_folder)

        res = 0.0

        phi.updateOld()

        D.setValue(D_value(phi, ele, tra_to_cut, cmask, drmask_not))
        C.setValue(C_value(phi, ele, tra_to_cut, cmask, drmask_not))

        for r in range(max_sweeps):
            resOld = res

            res = eq.sweep(var=phi, dt=dt)  # solve linearization of PDE

            # print "sum of Ds: ", np.sum(D.value)/1e8
            # print "average wt: ", np.average(phi.value-ele)

            # print 'residue diference:    ', res - resOld

            if abs(res - resOld) < 1e-7:
                break  # it has reached to the solution of the linear system

        if solve_mode == 'transient':  # solving in steadystate will remove water only at the very end
            if remove_ponding_water:
                s = np.where(phi.value > ele, ele,
                             phi.value)  # remove the surface water. This also removes phi in those masked values (for the plot only)
                phi.setValue(s)  # set new values for water table

        if (D.value < 0.).any():
            print("Some value in D is negative!")

        # For some plots
        avg_wt.append(np.average(phi.value - ele))
        wt_track_drained.append(
            (phi.value - ele).reshape(ny, nx)[track_WT_drained_area])
        wt_track_notdrained.append(
            (phi.value - ele).reshape(ny, nx)[track_WT_notdrained_area])

        """ Volume of dry peat calc."""
        not_peat = np.ones(shape=peat_type_mask.shape)  # Not all soil is peat!
        not_peat[peat_type_mask == 4] = 0  # NotPeat
        not_peat[peat_type_mask == 5] = 0  # OpenWater
        peat_vol_weights = utilities.PeatV_weight_calc(
            np.array(~dr * catchment_mask * not_peat, dtype=int))
        dry_peat_volume = utilities.PeatVolume(
            peat_vol_weights, (ele - phi.value).reshape(ny, nx))
        cumulative_Vdp = cumulative_Vdp + dry_peat_volume
        # Append timestep data to the list
        timestep_data.append([
            d,
            P[d] if isinstance(P, np.ndarray) else P,
            ET[d] if isinstance(ET, np.ndarray) else ET,
            dry_peat_volume,
            wt_track_drained[-1],
            wt_track_notdrained[-1],
            avg_wt[-1]
        ])
        print("avg_wt  = ", np.average(phi.value - ele))
        print("Cumulative vdp = ", cumulative_Vdp)

        # Calculate n-day average WTD
        if d >= n_day_avg - 1:  # Start calculating after n_day_avg - 1 days
            avg_wtd_nday.append(np.mean(avg_wt[-n_day_avg:])) # Average of the last n_day_avg values

    if solve_mode == 'steadystate':  # solving in steadystate we remove water only at the very end
        if remove_ponding_water:
            s = np.where(phi.value > ele, ele,
                             phi.value)  # remove the surface water. This also removes phi in those masked values (for the plot only)
            phi.setValue(s)  # set new values for water table

        # Areas with WT <-1.0; areas with WT >0
        # plot_raster_by_value((ele-phi.value).reshape(ny,nx), title="ele-phi in the end, colour keys", bottom_value=0.5, top_value=0.01)

    # Get raster arrays after the computation
    rast_D_after, rast_cwl_after, rast_dem_after, rast_elev_phi_after = get_rasters(
        ny, nx, peat_type_mask, phi.value - ele, httd, tra_to_cut, cmask, drmask_not, ele, wt_canal_arr, dr, catchment_mask)

    if plotOpt:

        # Plot raster arrays after the computation
        big_4_raster_plot(title='After the computation',
                          raster1=rast_D_after,
                          raster2=rast_cwl_after,
                          raster3=rast_dem_after,
                          raster4=rast_elev_phi_after,
                          output_folder=output_folder)

        plt.show()

        plot_hydro_params(httd, avg_wt, output_folder=output_folder)
        
        plt.show()

        # plot surface in cross-section
        ele_with_can = copy.copy(ele).reshape(ny, nx)
        ele_with_can = ele_with_can * catchment_mask
        ele_with_can[wt_canal_arr > 0] = wt_canal_arr[wt_canal_arr > 0]
        plot_line_of_peat(ele_with_can, y_value=y_value, title="cross-section", nx=nx,
                         ny=ny, label="surface", color='peru', linewidth=2.0,
                         output_folder=output_folder)

        plt.show()

    #    change_in_canals = (ele-phi.value).reshape(ny,nx)*(drmask.value.reshape(ny,nx)) - ((ele-H)*drmask.value).reshape(ny,nx)
    #    resulting_phi = phi.value.reshape(ny,nx)
    avg_wt_over_time = np.mean(np.array(avg_wt))
    print("Length of timestep_data:", len(timestep_data))
    return (
        cumulative_Vdp,
        dry_peat_volume,
        avg_wt_over_time,
        wt_track_drained,
        wt_track_notdrained,
        timestep_data,
        rast_D_before,
        rast_cwl_before,
        rast_dem_before,
        rast_elev_phi_before,
        rast_D_after,
        rast_cwl_after,
        rast_elev_phi_after,
        avg_wtd_nday  # Return n-day average WTD
    )
