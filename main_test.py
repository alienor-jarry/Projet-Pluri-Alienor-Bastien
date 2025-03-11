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
el.set_parameters(mesh_step=1.5, r1=1, d1=2,d2=2,t1=28,t2=26,t3=32)


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

    ##for i in range(1): ## if needed to check for multiple angle values at once, uncomment block below and change 'in range' value.
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

        ##for j in range(1): ## if needed to check for multiple angle values at once, uncomment block below and change 'in range' valule. 
        # initial=[0,0.5*j+0.5]

    
            ##
            ### Calculation of image point (zi)
            ##
    halfz= 10   ## A number > 0 but small so it will search the root of marginal after index 10 in the marginal table
    zi= np.argmin(np.abs(Marginal[0][halfz:]))+halfz ## Calculate the point where the marginal cross the axis
    
    M=el.Magnification(Principal[0],Marginal[0])   ## Give Magnification: when marginal ray cross z=0
    

            ##
            ### Calculation of spherical aberration and chromatic aberration. 
            ### Cs is the spherical aberration coefficient, DXS its delta X. Cc is the chromatic aberration coefficient, DXC its delta X.
            ##
    Cs,DXS=el.Spherical_Aberration(Marginal[0][0:zi],Marginal[1][0:zi],Potential_evaluated_z_TOT[0:zi],D2_evaluated_z_TOT[0:zi],D4_evaluated_z_TOT[0:zi],z_axis[0:zi],initial) ## Calculatation of Spherical Aberration coefficient and its delta X (transverse)
    

    Cc,DXC=el.Chromatic_Aberration(Marginal[0][0:zi],Marginal[1][0:zi],Potential_evaluated_z_TOT[0:zi],z_axis[0:zi],initial,deltaE0)  ## Calculatation of Chromatic Aberration coefficient and its delta X (transverse)
    
            
            ##
            ### Calculation of spot Diameter
            ##
    
    spot=el.Spot_Diameter(DXS,DXC,M,2)
    
    return z_axis[zi],DXS,DXC,spot



focale,DXS,DXC,spot=simu(27000,False,[0,.5],el)##Values of the original lens before optimizing

def autofocus(efl):
    return (efl-focale)

print(autofocus(focale))

def jacob_focus(F,var,h):
 elt=einzel_lens(show_mesh=True)
 n=len(var)  
 J=np.zeros(n)
 for j in range(n):
    var_step=var.copy()
    var_step[j]+= h
    elt.set_parameters(mesh_step=1.5, r1=1, d1=var_step[1],d2=var_step[2],t1=28,t2=26,t3=32)
    foc,S,C,spot_size=simu(var_step[0],False,[0,.5],elt)
    focus=(autofocus(foc))
    ##print(focus)
    J[j]=(focus-F)/h
 return J

def DLS(F,p,Var,J):
    
    a=np.dot(np.transpose(J),J)+p*np.eye(len(Var))
    b=np.dot(-np.transpose(J),F)
    Delta_Var=np.linalg.solve(a,b)

    return Delta_Var


def optim(p,F,var):
   f_start=F
   v0=var
   elt=einzel_lens(show_mesh=True)
   k=0
   Vec_f=Vec_f.append(f_start)
   for k in range(10):
      v=v0
      J=jacob_focus(f_start,v,h=1)
      D_Var=DLS(f_start,p,v,J)
      v=v+D_Var
      elt.set_parameters(mesh_step=1.5, r1=1, d1=v[1],d2=v[2],t1=28,t2=26,t3=32)
      ft,xst,xct,spt=simu(v[0],False,[0,1],elt)
      f=autofocus(ft)
      while (f**2)>(f_start**2) and p<10000:
         k+=1
         p=10*p
         D_Var=DLS(f_start,p,v,J)
         v=v+D_Var
         elt.set_parameters(mesh_step=1.5, r1=1, d1=v[1]+D_Var[1],d2=v[2]+D_Var[2],t1=28,t2=26,t3=32)
         ft,xst,xct,spt=simu(v[0]+D_Var[0],False,[0,1],elt)
         f=autofocus(ft)
         Vec_f=Vec_f.append(f)
      
      if (f**2)<=(f_start**2):
         while True:
            k+=1
            p=0.1*p
            D_Var=DLS(f_start,p,v,J)
            v=v+D_Var
            elt.set_parameters(mesh_step=1.5, r1=1, d1=v[1]+D_Var[1],d2=v[2]+D_Var[2],t1=28,t2=26,t3=32)
            ft,xst,xct,spt=simu(v[0]+D_Var[0],False,[0,1],elt)
            f=autofocus(ft)
            Vec_f=Vec_f.append(f)
            if not ((Vec_f[k]**2)<(Vec_f[k-1]**2)):
               break
         
         v0=v-D_Var
         p=p*10
         f_start=Vec_f[k-1]
      print(f)
   print(ft)


el0 = einzel_lens(show_mesh=True)
el0.set_parameters(mesh_step=1.5, r1=1, d1=20,d2=20,t1=28,t2=26,t3=32)
fz,xs,xc,sp=simu(27000,False,[0,.5],el0)
##print(fz)
##print(jacob_focus(autofocus(focale),[27000,20,20],h=1))
##D_var=DLS(autofocus(fz),1,[27000,20,20],jacob_focus(autofocus(fz),[27000,20,20],1))
##elf=einzel_lens(show_mesh=True)
##elf.set_parameters(mesh_step=1.5, r1=1, d1=20+D_var[1],d2=20+D_var[2],t1=28,t2=26,t3=32)
##ff,xsf,xcs,spf=simu(27000+D_var[0],False,[0,.5],elf)
##print(ff)
optim(0.1,fz,[27000,20,20])
