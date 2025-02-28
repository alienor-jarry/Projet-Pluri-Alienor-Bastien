import math
import os
import sys
import numpy as np
from scipy.integrate import *
from matplotlib import pylab as plt
import bempp.api                    ## Need to add first, second, third and fourth derivate of BEMPP
import gmsh
import pyopencl 
import RK45
from einzel import einzel_lens

el = einzel_lens(show_mesh=True)
el.set_parameters(r1=)

def simu(V, Fichier, I, el):
    el.setup()    ## Setup the basic configuration for BEMPP

    ###  Setup the basic inital condition.
    
    Fichier_dat=Fichier    ## Switch to True if you already have a .dat file in the directory or False if you need one
    Write_file_name = '\\Data\\JBLENS.dat'
    read_file_name = '\\Data\\Lentille_Simion_10000points_-100_1700.dat'
    file_path_msh='\\msh\\einzel_jb.msh'  ## Use to define path of the .msh file (where it will be created and read from)
    MeshSizeFactor=1.7      ## Decide the finess of the msh, the lower the finest
    n_z_points = 10000       ## The number of points on which BEMPP will calculate data
    zmin=-100               ## Starting point of data calculated (potential and its derivates)
    zmax=1700                ## Ending point of data calculated 


    E0=30000           # Energy of ion beam
    Qref=1             # Charge of reference: 1 for +charge (proton), -1 for -charge(electron)
    deltaE0=7          # Dispersion of energy of ion beam
    Voltage_E1 = 0     # Voltage of electrode 1
    Voltage_E2 = V     # Voltage of electrode 2
    Voltage_E3 = 0     # Voltage of electrode 3
    initial=I          # Inital condition of ray departure : first its height, second its angle in degree

    ###

    ##
    ### If a .dat is not in the directory, it will calculate data using BEM method with BEMPP
    ### It will first create the geometry (.msh), calculate function space and boundary operators for BEMPP use.
    ### Then it calculate for each electrode the dirichlet and neumann values that is used to calculate potential and its negatives derivates.
    ### It is then written in a .dat file that can be used for later time without need to re-calculate all data with BEM.
    ##

    if not(Fichier_dat):
        # Geometrie(True, MeshSizeFactor,file_path_msh) #If you want a cylinder around, switch Masse to True 
        el.GeometrieTest(MeshSizeFactor,file_path_msh) #If you want a cylinder around, switch Masse to True 
        (dp0_space,p1_space,identity,slp,dlp)=el.pre_dirichlet(file_path_msh)   # Calcul function space and boundary operators from msh created at file_path_msh
        
        ##  Create variable in which will be stored data
        neumann=[]
        Potential_evaluated_z=[]
        E_evaluated_z=[]
        D2_evaluated_z=[]
        D3_evaluated_z=[]
        D4_evaluated_z=[]
        dz=[]
        z_axis=[]
        for i in range(3):  # For each electrode from 1 to 3
            neumann.append(el.Dirichlet_Neumann(p1_space,identity,slp,dlp,i+1))    # Calcul Dirichelt and neumanns values for electrode i+1 (from 1 to 3)
            values=el.Potential_Derivate(dp0_space,neumann[i],zmin,zmax,n=n_z_points)  # Calcul potential and its negatives derivates from zmin to zmax
            dz.append(values[0])                        # Add dz to variable from Potential_Derivate (shoud be the same for each i)
            z_axis.append(values[1])                    # Add z_axis values to variable from Potential_Derivate (should be the same for each i)
            Potential_evaluated_z.append(values[2])     # Add potential on z-axis to variable from Potential_Derivate
            E_evaluated_z.append(values[3])             # Add first negative derivate on z-axis to variable from Potential_Derivate
            D2_evaluated_z.append(values[4])            # Add second negative derivate on z-axis to variable from Potential_Derivate
            D3_evaluated_z.append(values[5])            # Add third negative derivate on z-axis to variable from Potential_Derivate
            D4_evaluated_z.append(values[6])            # Add fourth negative derivate on z-axis to variable from Potential_Derivate

        z_axis=z_axis[0]    #Only need one array of z_axis values. [0] is the same as [1] or [2]
        el.write_data(Write_file_name,z_axis,Potential_evaluated_z,E_evaluated_z,D2_evaluated_z,D3_evaluated_z,D4_evaluated_z) #Write data in the Write_file_name

    ##
    ### If a .dat is in the directory, it will simply read the data in the .dat and store in variable that will be used.
    ##

    if (Fichier_dat): #If a .dat file is already in the directory, it will stock the data in correct variable 
        Potential_evaluated_z=[0,0,0]
        E_evaluated_z=[0,0,0]
        D2_evaluated_z=[0,0,0]
        D3_evaluated_z=[0,0,0]
        D4_evaluated_z=[0,0,0]
        (r,z_axis,Potential_evaluated_z[0],E_evaluated_z[0],D2_evaluated_z[0],D3_evaluated_z[0],D4_evaluated_z[0],Potential_evaluated_z[1],E_evaluated_z[1],D2_evaluated_z[1],D3_evaluated_z[1],D4_evaluated_z[1],Potential_evaluated_z[2],E_evaluated_z[2],D2_evaluated_z[2],D3_evaluated_z[2],D4_evaluated_z[2])=el.read_data(read_file_name) # Read the data written in the .dat data file
    
    ##
    ### Starting here, the script goes through, even if data was not written in .dat before (it was created)
    ##

    for i in range(1): ## if needed to check for multiple angle values at once, uncomment block below and change 'in range' value.
        # Voltage_E2=-55000-i*1000    
        ##
        ### The derivate here are, all, negative derivative of potential
        ##
        Potential_evaluated_z_TOT_0,E_evaluated_z_TOT,D2_evaluated_z_TOT,D3_evaluated_z_TOT,D4_evaluated_z_TOT=el.SumAndVoltage(Potential_evaluated_z,E_evaluated_z,D2_evaluated_z,D3_evaluated_z,D4_evaluated_z,[Voltage_E1,Voltage_E2,Voltage_E3]) # Sum of voltage created by differents electrode and multiplied from unitary voltage to chosen one (Voltage_E*)
        Potential_evaluated_z_TOT=E0/Qref-Potential_evaluated_z_TOT_0 ## We use the kinetic energy divided by e. Acceleration voltage divided by charge reference minus potential (calculated by BEMPP). It is necessary to use it for formulae such as Spheric aberration.
        
        ### You can uncomment below if you want to plot potential and its negative derivates
        # Plot(z_axis,Potential_evaluated_z_TOT,5)
        # Plot(z_axis,E_evaluated_z_TOT,6)
        # Plot(z_axis,D2_evaluated_z_TOT,7)
        # Plot(z_axis,D3_evaluated_z_TOT,8)
        # Plot(z_axis,D4_evaluated_z_TOT,9)

        ##
        ### Calculation of principal ray and marginal ray solving paraxial equation with a Runge Kunta 4. It return y and y' for Principal and Marginal. (y : its height along the z axis, y': its derivate (slope))
        ##
        [Principal,Marginal]=el.PrincipalAndMarginalRays(z_axis,Potential_evaluated_z_TOT,E_evaluated_z_TOT,D2_evaluated_z_TOT)
        # Plot(z_axis,Principal[0],10)
        # Plot(z_axis,Marginal[0],10)

        ## 
        ### Calculation of trajectory using principal and marginal ray with initial condition
        ##
        r2=el.Trajectory(Principal[0],Marginal[0],initial) ## Trajectory of a ray with initial condition
        # Plot(z_axis,r2,11)

        for j in range(1): ## if needed to check for multiple angle values at once, uncomment block below and change 'in range' valule. 
        # initial=[0,0.5*j+0.5]

            print("")
            print("Voltage Electrode 2 at", Voltage_E2,"V: ")
            print("Initial state : h =",initial[0]," and alpha =",initial[1],"degree")
            ##
            ### Calculation of image point (zi)
            ##
            halfz= 10   ## A number > 0 but small so it will search the root of marginal after index 10 in the marginal table
            zi= np.argmin(np.abs(Marginal[0][halfz:]))+halfz ## Calculate the point where the marginal cross the axis
            print("Image Plan : ",z_axis[zi])
            M=el.Magnification(Principal[0],Marginal[0])   ## Give Magnification: when marginal ray cross z=0
            print("Magnification : ", M)

            ##
            ### Calculation of spherical aberration and chromatic aberration. 
            ### Cs is the spherical aberration coefficient, DXS its delta X. Cc is the chromatic aberration coefficient, DXC its delta X.
            ##
            Cs,DXS=el.Spherical_Aberration(Marginal[0][0:zi],Marginal[1][0:zi],Potential_evaluated_z_TOT[0:zi],D2_evaluated_z_TOT[0:zi],D4_evaluated_z_TOT[0:zi],z_axis[0:zi],initial) ## Calculatation of Spherical Aberration coefficient and its delta X (transverse)
            print("Spherical Abberation coefficient :", Cs)
            print("Coef CS JB :",get_CoeffSpheric(Potential_evaluated_z_TOT[0:zi],E_evaluated_z_TOT[0:zi],D2_evaluated_z_TOT[0:zi],z_axis[0:zi],Marginal[0][0:zi],Marginal[1][0:zi])) ## Spherical Coefficient aberration (JB)
            print("Spherical Delta X : ", DXS)

            Cc,DXC=el.Chromatic_Aberration(Marginal[0][0:zi],Marginal[1][0:zi],Potential_evaluated_z_TOT[0:zi],z_axis[0:zi],initial,deltaE0)  ## Calculatation of Chromatic Aberration coefficient and its delta X (transverse)
            print("Chromatic Abberation coefficient :", Cc)
            print("Chromatic Delta X : ", DXC)
            
            ##
            ### Calculation of spot Diameter
            ##
            print("Spot Diameter :",el.Spot_Diameter(DXS,DXC,M,2))
