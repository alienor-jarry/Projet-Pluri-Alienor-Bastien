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




class einzel_lens:
    def __init__(self, show_mesh=False):
        self._show_mesh = show_mesh
        self.Focale_Voulue=None
    
    def set_parameters(self, mesh_step, t1, d1, d2, r1, r2=None, r3=None, t2=None, t3=None, rc=None, tc=None):
        self._d1=d1                   ## d => espacement électrode
        self._d2=d2                   ## t => longueur électrode
        self._r1=r1                   ## r => rayon électrode
        self._t1 = t1
        self._t2=t2
        self._t3=t3
        if r2 == None :
            self._r2 = r1
        if r3 == None :
            self._r3 = r1
        if t2 == None :
            self._t2 = t1
        if t3 == None :
            self._t3 = t1
        if tc == None :
            self._tc = t1*1.5 # a definir
        if rc == None :
            self._rc = t1 # a definir 

    def setup(self):
     """
     Setup the basic need. 
     Might need to change the bempp.core.opencl_kernels.set_default_cpu_device(1, 0) to your computer configuration.
     To do that, install openCL for your processor and use bempp.core.opencl_kernels.show_available_platforms_and_devices() to see your own configuration.

     """

     bempp.core.opencl_kernels.set_default_cpu_device(1, 0) 
     self.virtual_env_path = sys.prefix  # Use to define gmsh path
     self.relative_path = "\Scripts\gmsh.bat" # Use to define gmsh path
     full_path = self.virtual_env_path+self.relative_path # Use to define gmsh path
     bempp.api.GMSH_PATH = full_path #Give bempp gmsh path
     #gmsh.option.setString('Geometry.OCCTargetUnit', 'M')
     scal = 1e-3

    def GeometrieTest(self, MeshSizeX=1.5, file_path='\\msh\\test.msh'):
    
     ## Lens Geometry for simple Einzel lens: Simion lens in particular.

    
     offset=18
     lc = 0.15
     rayon=1
     longueur_electrode_centrale=self._t2
     l1=self._t1
     l3=self._t3   
     espacement1=self._d1
     espacement2=self._d2

     gmsh.initialize()
     gmsh.clear()
     gmsh.model.add("three electrode")

     #1ST ELECTRODE
     gmsh.model.occ.addPoint(0, 0+offset, 0, lc, 1)
     gmsh.model.occ.addPoint(0, rayon+offset, 0, lc, 2)
     gmsh.model.occ.addPoint(0, rayon+offset, l1, lc, 3)
     gmsh.model.occ.addPoint(0, 0+offset, l1, lc, 4)
     gmsh.model.occ.addLine(1, 2, 1)
     gmsh.model.occ.addLine(2, 3, 2)
     gmsh.model.occ.addLine(3, 4, 3)
     gmsh.model.occ.addLine(4, 1, 4)
     gmsh.model.occ.addCurveLoop([1, 2, 3, 4], 1)
     gmsh.model.occ.addPlaneSurface([1], 1)

     #2ND ELECTRODE
     gmsh.model.occ.addPoint(0, 0+offset, espacement1+l1, lc, 5)
     gmsh.model.occ.addPoint(0, rayon+offset, espacement1+l1, lc, 6)
     gmsh.model.occ.addPoint(0, rayon+offset, longueur_electrode_centrale+espacement1+l1, lc, 7)
     gmsh.model.occ.addPoint(0, 0+offset, longueur_electrode_centrale+espacement1+l1, lc, 8)
     gmsh.model.occ.addLine(5, 6, 5)
     gmsh.model.occ.addLine(6, 7, 6)
     gmsh.model.occ.addLine(7, 8, 7)
     gmsh.model.occ.addLine(8, 5, 8)
     gmsh.model.occ.addCurveLoop([5, 6, 7, 8], 2)
     gmsh.model.occ.addPlaneSurface([2], 2)

     #3ND ELECTRODE
     gmsh.model.occ.addPoint(0, 0+offset, espacement1+espacement2+l1+longueur_electrode_centrale, lc, 9)
     gmsh.model.occ.addPoint(0, rayon+offset, espacement1+espacement2+l1+longueur_electrode_centrale, lc, 10)
     gmsh.model.occ.addPoint(0, rayon+offset, espacement1+espacement2+l1+longueur_electrode_centrale+l3, lc, 11)
     gmsh.model.occ.addPoint(0, 0+offset, espacement1+espacement2+l1+longueur_electrode_centrale+l3, lc, 12)
     gmsh.model.occ.addLine(9, 10, 9)
     gmsh.model.occ.addLine(10, 11, 10)
     gmsh.model.occ.addLine(11, 12, 11)
     gmsh.model.occ.addLine(12, 9, 12)
     gmsh.model.occ.addCurveLoop([9, 10, 11, 12], 3)
     gmsh.model.occ.addPlaneSurface([3], 3)

     #Revolution
     e1=gmsh.model.occ.revolve([(2, 1)], 0, 0, 0, 0, 0, 1, 2*math.pi)
     e2=gmsh.model.occ.revolve([(2, 2)], 0, 0, 0, 0, 0, 1, 2*math.pi)
     e3=gmsh.model.occ.revolve([(2, 3)], 0, 0, 0, 0, 0, 1, 2*math.pi)

     #Add of cylinder around
     a=gmsh.model.occ.addCylinder(0, 0, -longueur_electrode_centrale, 0, 0, 1.5*(l1+l3+longueur_electrode_centrale), 1.2*(rayon+offset))
     b=gmsh.model.occ.addCylinder(0, 0, -longueur_electrode_centrale, 0, 0, 1.5*(l1+l3+longueur_electrode_centrale), 1.1*(rayon+offset))
     res=gmsh.model.occ.cut([(3,a)],[(3,b)])

     #Clean-up and synchronization
     gmsh.model.occ.remove(gmsh.model.occ.getEntities(2), True)
     gmsh.model.occ.synchronize()


     #Add of physical group
     Post=gmsh.model.occ.getEntities(2)
     set_post=set(Post)
     set_e1=set(e1)
     set_e2=set(e2)
     set_e3=set(e3)
     intersection_e1=list(set_post.intersection(set_e1))
     intersection_e2=list(set_post.intersection(set_e2))
     intersection_e3=list(set_post.intersection(set_e3))

     list0=[]
     for i in gmsh.model.occ.getEntities(0):
         list0.append(i[1])
     list1=[]
     for i in gmsh.model.occ.getEntities(1):
         list1.append(i[1])
     liste1=[]
     for i in intersection_e1:
         liste1.append(i[1])
     liste2=[]
     for i in intersection_e2:
         liste2.append(i[1])
     liste3=[]
     for i in intersection_e3:
         liste3.append(i[1])

     #print(liste1,"   " , liste2, "   ",liste3)
     # gmsh.model.addPhysicalGroup(0,list0,1,name="Point")
     # gmsh.model.addPhysicalGroup(1,list1,1,name="Line")
     # gmsh.model.addPhysicalGroup(2, liste3, 3, name="Electrode 3")
     # gmsh.model.addPhysicalGroup(2, liste2, 2,name="Electrode 2")  
     # gmsh.model.addPhysicalGroup(2, liste1, 1,name="Electrode 1")
     # gmsh.model.addPhysicalGroup(2, [19,20,21,22], 4)
     gmsh.model.addPhysicalGroup(0,list0,1)
     gmsh.model.addPhysicalGroup(1,list1,1)
     gmsh.model.addPhysicalGroup(2, [4,5,6,7], 1)
     gmsh.model.addPhysicalGroup(2, [8,9,10,11], 2)  
     gmsh.model.addPhysicalGroup(2, [12,13,14,15], 3)
     gmsh.model.addPhysicalGroup(2, [19,20,21,22], 4)

     # gmsh.model.addPhysicalGroup(3, [1], tag=3,name="Electrode 3")
     # gmsh.model.addPhysicalGroup(3, [2], tag=2,name="Electrode 2")
     # gmsh.model.addPhysicalGroup(3, [3], tag=1,name="Electrode 1")
     # gmsh.model.addPhysicalGroup(3, [4], tag=4,name="Masse")
     # gmsh.model.addPhysicalGroup(3, [box], tag=box,name="Sample")


     #Set of Mesh size 
     gmsh.option.set_number("Mesh.CharacteristicLengthMin",1.7)
     gmsh.option.set_number("Mesh.CharacteristicLengthMax",1.7)
     gmsh.option.setNumber("Mesh.SaveAll", 1)
     gmsh.option.setNumber("Mesh.MeshSizeFactor", MeshSizeX)

     #Mesh of 2D 
     gmsh.model.mesh.generate(2)

     #Writing as jb.msh in a folder named 'msh' in the directory of this file .py
     file_path_all=os.path.dirname(os.path.realpath(__file__))+file_path
     gmsh.write(file_path_all)
     ##gmsh.fltk.run()
     gmsh.finalize()
    
    def pre_dirichlet(self,file_path):
     """
     Importation of grid made with gmsh.
     Importation of boundary of said grid.
     Calculation of function space and boundary operators needed for calculation of Dirichlet and Neumanns values respectively.

     """
     file_path_all=os.path.dirname(os.path.realpath(__file__))+file_path
     grid = bempp.api.import_grid(file_path_all)
     ((xmin, xmax),(ymin, ymax),(zmin,zmax)) =grid.bounding_box
     #xmin, ymin, zmin, xmax, ymax, zmax = -100. ,-100., -100., 100., 100., 100.
     grid.plot()

     # function spaces for BEMPP
     dp0_space = bempp.api.function_space(grid, "DP", 0)
     p1_space = bempp.api.function_space(grid, "P", 1)

     # boundary operators
     identity = bempp.api.operators.boundary.sparse.identity(
        p1_space, p1_space, dp0_space)
     dlp = bempp.api.operators.boundary.laplace.double_layer(
        p1_space, p1_space, dp0_space, device_interface="opencl")
     slp = bempp.api.operators.boundary.laplace.single_layer(
        dp0_space, p1_space, dp0_space, device_interface="opencl")
    
     return(dp0_space,p1_space,identity,slp,dlp)
    
    def Dirichlet_Neumann(self,p1_space,identity,slp,dlp,numero_electrode):
     """
     Setting the Dirichlet values of 'numero_electrode' to 1 and others electrode to 0. 
     Electrode 1 is the last electrode the ion will see . 
     Calculation of Neumann values for these dirichlet Data
    
     """
     # boundary conditions
     @bempp.api.real_callable
     def dirichlet_data(x, n, domain_index, result):

        # Assign voltage to electrodes
        if domain_index==numero_electrode:
            result[0]=1
        else:
            result[0]=0

     dirichlet_fun = bempp.api.GridFunction(p1_space, fun=dirichlet_data)
     dirichlet_fun.plot()

     rhs = (-0.5 * identity + dlp) * dirichlet_fun
     neumann_fun, info = bempp.api.linalg.cg(slp, rhs, tol=1E-3)

     return neumann_fun  

    def Potential_Derivate(self,dp0_space,neumann_fun,zmin,zmax,n):
     """
     Calculation of potential and its derivates up to 4nd order using BEM method, BEMPP.
     Careful, the derivate obtained is the negative derivate of potential.
     E=-Potential' ; D2=-Potential'' ; D3=-Potential''' ; D4=-Potential''''

     """
     # Derivatives on optical axis
     # numerical diff. distance
     dz = (zmax-zmin)/n

     grid_z =  np.mgrid[0:1:1j,0:1:1j,zmin:zmax:n*1j]
     points_z = np.vstack((grid_z[0].ravel(), grid_z[1].ravel(), grid_z[2].ravel()))
     zzz=grid_z[2].ravel()#[0:n])

     #Get the potential
     slp_pot = bempp.api.operators.potential.laplace.single_layer(dp0_space, points_z,device_interface="opencl") #Calculation of potential operator.
     u_zb = -slp_pot * neumann_fun
     Potential_eval_z=u_zb.reshape(n)

     #Get the electrid field
     E = bempp.api.operators.potential.laplace.single_layer_gradient(dp0_space, points_z, device_interface="opencl")
     E_eval = -E * neumann_fun
     E_eval_z=E_eval[-1,:] 

     #Get the minus second derivate
     D2 = bempp.api.operators.potential.laplace.single_layer_2nd_deriv(dp0_space, points_z, device_interface="opencl")
     D2_eval = -D2 * neumann_fun
     D2_eval_z=D2_eval[-1,:]

     #Get the minus third derivate
     D3 = bempp.api.operators.potential.laplace.single_layer_3rd_deriv(dp0_space, points_z, device_interface="opencl")
     D3_eval = -D3 * neumann_fun
     D3_eval_z=D3_eval[-1,:]

     #Get the minus fourth derivate
     D4 = bempp.api.operators.potential.laplace.single_layer_4th_deriv(dp0_space, points_z, device_interface="opencl")
     D4_eval = -D4 * neumann_fun
     D4_eval_z=D4_eval[-1,:]

     return (dz,zzz,Potential_eval_z,E_eval_z,D2_eval_z,D3_eval_z,D4_eval_z)

    def write_data(self,file_name,zzz,Potential_eval_z,E_eval_z,D2_eval_z,D3_eval_z,D4_eval_z):
     """"
     Write the data as .dat with a column of r=0,z position, potential when Electrode1=1, D1,D2,D3,D4, potential when E2=1,D1,D2,D3,D4 potential when E3=1,D1,D2,D3,D4.
     Careful, the derivate written is the negative derivate of potential.
     It write data at file_name.
    
     """
     r=np.zeros(len(zzz))
     tab=np.array([r,zzz,Potential_eval_z[0],E_eval_z[0],D2_eval_z[0],D3_eval_z[0],D4_eval_z[0],Potential_eval_z[1],E_eval_z[1],D2_eval_z[1],D3_eval_z[1],D4_eval_z[1],Potential_eval_z[2],E_eval_z[2],D2_eval_z[2],D3_eval_z[2],D4_eval_z[2]]).transpose()
     np.savetxt(os.path.dirname(os.path.realpath(__file__))+ file_name,tab,delimiter='\t',fmt='%.12f',comments="")

    def read_data(self,file_name):
     """
     Read a data set .dat with a column of r=0,z position, potential when Electrode1=1 then D1,D2,D3,D4, potential when E2=1 then D1,D2,D3,D4, potential when E3=1 then D1,D2,D3,D4
     Careful, the derivate written is the negative derivate of potential.
     It read data from file_name.

     """
     # Charger le fichier .dat en utilisant pandas en spécifiant le délimiteur comme une tabulation
     data = np.loadtxt(os.path.dirname(os.path.realpath(__file__))+ file_name, delimiter="\t")
     #data = pd.read_csv((os.path.dirname(os.path.realpath(__file__))+ file_name), delimiter="\t", header=None)
     # Accéder à chaque colonne
     r = data[:, 0]  # r=0 : Colonne 0
     z = data[:, 1]  # z position : Colonne 1
     ue1 = data[:, 2]  # potential when Electrode1=1 : Colonne 2
     d1e1 = data[:, 3]  # D1 when Electrode1=1 : Colonne 3
     d2e1 = data[:, 4]  # D2 when Electrode1=1 : Colonne 4
     d3e1 = data[:, 5]  # D3 when Electrode1=1 : Colonne 5
     d4e1 = data[:, 6]  # D4 when Electrode1=1 : Colonne 6
     ue2 = data[:, 7]  # potential when Electrode2=1 : Colonne 7
     d1e2 = data[:, 8]  # D1 when Electrode2=1 : Colonne 8
     d2e2 = data[:, 9]  # D2 when Electrode2=1 : Colonne 9
     d3e2 = data[:, 10]  # D3 when Electrode2=1 : Colonne 10
     d4e2 = data[:, 11]  # D4 when Electrode2=1 : Colonne 11
     ue3 = data[:, 12]  # potential when Electrode3=1 : Colonne 12
     d1e3 = data[:, 13]  # D1 when Electrode3=1 : Colonne 13
     d2e3 = data[:, 14]  # D2 when Electrode3=1 : Colonne 14
     d3e3 = data[:, 15]  # D3 when Electrode3=1 : Colonne 15
     d4e3 = data[:, 16]  # D4 when Electrode3=1 : Colonne 16
     return (r, z, ue1, d1e1, d2e1, d3e1, d4e1, ue2, d1e2, d2e2, d3e2, d4e2, ue3, d1e3, d2e3, d3e3, d4e3)

    def SumAndVoltage(self,u_z_line,E_eval_z,D2_eval_z,D3_eval_z,D4_eval_z,voltage):
     """
     Multiplication of unitary voltage to choosen voltage and make the sum. 
     Voltage is a table like [0,30000,0] where Electrode 1&3 are set to 0 and electrode 2 to 30000 for example.
     u_z_line is the potential. 
     E_eval_z is the electrid field. 
     D2_eval_z is the second derivative.
     D3_eval_z is the third derivative. 
     D4_eval_z is the fourth derivative. 
     Voltage is the table which contained the values to which you want to multiply the voltage and its derivatives.

     """
    
     #Sum of voltage and derivative
     u_z_linetot=u_z_line[0]*voltage[0]+u_z_line[1]*voltage[1]+u_z_line[2]*voltage[2]#uu_z_linefinal[0]+u_z_linefinal[1]+u_z_linefinal[2]
     E_eval_ztot=E_eval_z[0]*voltage[0]+E_eval_z[1]*voltage[1]+E_eval_z[2]*voltage[2]#E_eval_zfinal[0]+E_eval_zfinal[1]+E_eval_zfinal[2]
     D2_eval_ztot=D2_eval_z[0]*voltage[0]+D2_eval_z[1]*voltage[1]+D2_eval_z[2]*voltage[2]#D2_eval_zfinal[0]+D2_eval_zfinal[1]+D2_eval_zfinal[2]
     D3_eval_ztot=D3_eval_z[0]*voltage[0]+D3_eval_z[1]*voltage[1]+D3_eval_z[2]*voltage[2]#D2_eval_zfinal[0]+D2_eval_zfinal[1]+D2_eval_zfinal[2]
     D4_eval_ztot=D4_eval_z[0]*voltage[0]+D4_eval_z[1]*voltage[1]+D4_eval_z[2]*voltage[2]#D2_eval_zfinal[0]+D2_eval_zfinal[1]+D2_eval_zfinal[2]

     return (u_z_linetot,E_eval_ztot,D2_eval_ztot,D3_eval_ztot,D4_eval_ztot)

    def PrincipalAndMarginalRays(self,zzz,Ec,E_eval_ztot,D2_eval_ztot):
     """
     Calculation of Trajectory and plot it. Call a Runge Kunta method to solve paraxial equation.
     zzz is the z-axis values. 
     Ec correspond to the potential of acceleration (E0) minus potential on z-axis due to electrode (calculate with BEMPP). 
     E_eval_ztot is the total electrid field. 
     D2_eval_ztot is the total negative second derivative of potential. 

     RK45 return the trajectory (its height along the z axis) along with its derivate (slope)
     """
   
     dz=abs(zzz[0]-zzz[1])
     n=len(zzz)
     PrincipalRay=RK45.RK45([1,0],n,dz,Ec,E_eval_ztot,D2_eval_ztot) #We have Ec=E0-U so we take negative derivate of U 
     MarginalRay=RK45.RK45([0,1],n,dz,Ec,E_eval_ztot,D2_eval_ztot)

     ##
     ### If you want to write data in a file
     ##
     # tab=np.array([zzz,PrincipalRay[0],PrincipalRay[1],MarginalRay[0],MarginalRay[1]]).transpose()
     # file_name = '\\trajectory.dat'
     # np.savetxt(os.path.dirname(os.path.realpath(__file__))+ file_name,tab,delimiter='\t',fmt='%.12f',comments="")

     return (PrincipalRay,MarginalRay)

    def Trajectory(self,PrincipalRay,MarginalRay,initial_state):
     """
     Calculate the trajectory of a ray with particular inital condition.
     PrincipalRay is the trajectory of the principal ray.
     MarginalRay is the trajectory of the marginal ray.
     initial_state correspond to initial condition. It is a table like [0,3] where 0 correspond to the height of the ray when it starts. 3 correspond to the demi-angle of departure in degree.

     """
     anglerad=np.tan(np.deg2rad(initial_state[1]))
     f=initial_state[0]*PrincipalRay+MarginalRay*anglerad
     return (f)
    
    def Magnification(self,PrincipalRay,MarginalRay):
     """
     Close approximation of Magnification.
     PrincipalRay is the trajectory of the principal ray.
     MarginalRay is the trajectory of the marginal ray.

     """
     #zero_crossings = np.where(np.diff(np.sign(MarginalRay)))[0]
     halfz=10 ## A number > 0 but small so it will search the root of marginal after index 10 in the marginal table
     index_val= np.argmin(np.abs(MarginalRay[halfz:]))+halfz
     M=PrincipalRay[index_val]/PrincipalRay[0]
     return M

    def Spherical_Aberration(self,r,s,U,D2,D4,z,initial_state):
     """
     Calculate the spherical coefficient aberration and its delta X associated.
     It use the formulae of 'Unified theory for electrostatic round lenses multipoles lenses and deflectors' by Munro and Smith (1986). Formula (37) and (40) and table 1
     r correspond to the trajectory of the marginal ray.
     s correspond to the derivative of the marginal ray (its slope)
     U correspond to the potential. 
        Careful, it is not only the potenial calculated with BEMPP but it should be the potential of acceleration (E0) minus potential calulated with BEMPP. (for + charge like protons)
     D2 correspond to the negative second derivative of the potential.
     D4 correspond to the negative fourth derivative of the potential.
     z correspond to the z-axis values
     initial_state correspond to initial condition. It is a table like [0,3] where 0 correspond to the height of the ray when it starts. 3 correspond to the demi-angle of departure in degree with which will be calculated the deltaX.
    
        All these parameters except inital_state should be given to the function from object point to image point so that U[-1] correspond to potential at image point or that r[-1] should be null for example.

     """

     angleradinit=np.tan(np.deg2rad(initial_state[1])) # We use the tan of radian value instead of degree

     ### Switch name to follow Unified theory for electrostatic round lenses multipoles lenses and deflectors' by Munro and Smith (1986)
     xa = r 
     xa_prime = s
     xr_prime=angleradinit*xa_prime[-1] #Formula (37a)

     ## Method with D4 (formula (40))

     f=U**(1/2)*(-D4*xa**4/(32*U)+D2**2*xa**4/(32*U**2)+D2*xa**2*xa_prime**2/(4*U)+xa_prime**4/2) 
     Cs=simpson(f,x=z)/((U[-1]**(1/2)*xa_prime[-1]**4)) #Formula (40) of 'Unified theory for electrostatic round lenses multipoles lenses and deflectors' by Munro and Smith (1986).
     I=xr_prime**3*Cs #Table (1) of 'Unified theory for electrostatic round lenses multipoles lenses and deflectors' by Munro and Smith (1986).

     # Method with IBP (D2 max) (formula (43))
     # xa_2nd=np.gradient(xa_prime,z)
     # fbis=U**(1/2)*(D2**2*xa**4/(32*U**2)+D2*xa**2*xa_prime**2/(4*U)+xa_prime**4/2-D2/(32*U)*(3*D1**2*xa**4/(4*U**2)-D2*xa**4/(2*U)-4*D1*xa**3*xa_prime/U+4*xa**3*xa_2nd+12*xa_prime**2*xa**2))
     # Csbis=1/(xa_prime[-1]**4*U[-1]**(1/2))*(simpson(fbis,x=z)+(U[-1]**(1/2)*D2[-1]*xa[-1]**3*xa_prime[-1]/(32*U[-1]))-(U[0]**(1/2)*D2[0]*xa[0]**3*xa_prime[0]/(32*U[0])))
     # Ibis=xr_prime**3*Csbis

     return Cs,I

    def Chromatic_Aberration(self,r,s,U,z,initial_state,DELTAE0):
     """
     Calculate the chromatic coefficient aberration and its delta X associated.
     It use the formulae of 'Unified theory for electrostatic round lenses multipoles lenses and deflectors' by Munro and Smith (1986). Formula (48),(49),(50),(51)
     r correspond to the trajectory of the marginal ray.
     s correspond to the derivative of the marginal ray (its slope)
     U correspond to the potential. 
        Careful, it is not only the potenial calculated with BEMPP but it should be the potential of acceleration (E0) minus potential calulated with BEMPP. (for + charge like protons)
     z correspond to the z-axis values
     initial_state correspond to initial condition. It is a table like [0,3] where 0 correspond to the height of the ray when it starts. 3 correspond to the demi-angle of departure in degree with which will be calculated the deltaX.
     DELTAE0 correspond to the dispersion of energy of the ion beam.
        
        All these parameters except inital_state should be given to the function from object point to image point so that U[-1] correspond to potential at image point or that r[-1] should be null.

     """

     angleradinit=np.tan(np.deg2rad(initial_state[1])) # We use the tan of radian value instead of degree

     ### Switch name to follow Unified theory for electrostatic round lenses multipoles lenses and deflectors' by Munro and Smith (1986)
     xa = r
     xa_prime = s
     xa_2nd = np.gradient(xa_prime,z) #Formula 51
     xr_prime =angleradinit*xa_prime[-1] # Formula (37a)

     ### Coefficient Cc (50),(51) 
     f=xa*xa_2nd/U**(1/2)
     C1=U[-1]**(1/2)/(xa_prime[-1]**2)*simpson(f,x=z)

     ### Axial chromatic aberration delta X. Formula (49)
     deltaXc=C1*DELTAE0*xr_prime/U[-1]
     return -C1,deltaXc

    def Spot_Diameter(self,DXS,DXC,M,alpha):
     """
     Calcul of spot diameter: quadrature of magnification, spherical delta x and chromatic delta x
    
     """
     #alpha_rad=alpha*math.pi/180
     ds=M
     # dsa=Cs*(alpha_rad)**3
     # dc=deltaE0*Cc*alpha_rad/E0
     return np.sqrt(ds**2+DXS**2+DXC**2)
    


    def simu(self,V, Fichier, I):
     self.setup()    ## Setup the basic configuration for BEMPP

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
         self.GeometrieTest(MeshSizeFactor,file_path_msh) #If you want a cylinder around, switch Masse to True 
         (dp0_space,p1_space,identity,slp,dlp)=self.pre_dirichlet(file_path_msh)   # Calcul function space and boundary operators from msh created at file_path_msh
        
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
             neumann.append(self.Dirichlet_Neumann(p1_space,identity,slp,dlp,i+1))    # Calcul Dirichelt and neumanns values for electrode i+1 (from 1 to 3)
             values=self.Potential_Derivate(dp0_space,neumann[i],zmin,zmax,n=n_z_points)  # Calcul potential and its negatives derivates from zmin to zmax
             dz.append(values[0])                        # Add dz to variable from Potential_Derivate (shoud be the same for each i)
             z_axis.append(values[1])                    # Add z_axis values to variable from Potential_Derivate (should be the same for each i)
             Potential_evaluated_z.append(values[2])     # Add potential on z-axis to variable from Potential_Derivate
             E_evaluated_z.append(values[3])             # Add first negative derivate on z-axis to variable from Potential_Derivate
             D2_evaluated_z.append(values[4])            # Add second negative derivate on z-axis to variable from Potential_Derivate
             D3_evaluated_z.append(values[5])            # Add third negative derivate on z-axis to variable from Potential_Derivate
             D4_evaluated_z.append(values[6])            # Add fourth negative derivate on z-axis to variable from Potential_Derivate

         z_axis=z_axis[0]    #Only need one array of z_axis values. [0] is the same as [1] or [2]
         self.write_data(Write_file_name,z_axis,Potential_evaluated_z,E_evaluated_z,D2_evaluated_z,D3_evaluated_z,D4_evaluated_z) #Write data in the Write_file_name

     ##
     ### If a .dat is in the directory, it will simply read the data in the .dat and store in variable that will be used.
     ##

     if (Fichier_dat): #If a .dat file is already in the directory, it will stock the data in correct variable 
         Potential_evaluated_z=[0,0,0]
         E_evaluated_z=[0,0,0]
         D2_evaluated_z=[0,0,0]
         D3_evaluated_z=[0,0,0]
         D4_evaluated_z=[0,0,0]
         (r,z_axis,Potential_evaluated_z[0],E_evaluated_z[0],D2_evaluated_z[0],D3_evaluated_z[0],D4_evaluated_z[0],Potential_evaluated_z[1],E_evaluated_z[1],D2_evaluated_z[1],D3_evaluated_z[1],D4_evaluated_z[1],Potential_evaluated_z[2],E_evaluated_z[2],D2_evaluated_z[2],D3_evaluated_z[2],D4_evaluated_z[2])=self.read_data(read_file_name) # Read the data written in the .dat data file
    
     ##
     ### Starting here, the script goes through, even if data was not written in .dat before (it was created)
     ##

     ##for i in range(1): ## if needed to check for multiple angle values at once, uncomment block below and change 'in range' value.
         # Voltage_E2=-55000-i*1000    
         ##
         ### The derivate here are, all, negative derivative of potential
         ##
     Potential_evaluated_z_TOT_0,E_evaluated_z_TOT,D2_evaluated_z_TOT,D3_evaluated_z_TOT,D4_evaluated_z_TOT=self.SumAndVoltage(Potential_evaluated_z,E_evaluated_z,D2_evaluated_z,D3_evaluated_z,D4_evaluated_z,[Voltage_E1,Voltage_E2,Voltage_E3]) # Sum of voltage created by differents electrode and multiplied from unitary voltage to chosen one (Voltage_E*)
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
     [Principal,Marginal]=self.PrincipalAndMarginalRays(z_axis,Potential_evaluated_z_TOT,E_evaluated_z_TOT,D2_evaluated_z_TOT)
         # Plot(z_axis,Principal[0],10)
         # Plot(z_axis,Marginal[0],10)

         ## 
         ### Calculation of trajectory using principal and marginal ray with initial condition
         ##
     r2=self.Trajectory(Principal[0],Marginal[0],initial) ## Trajectory of a ray with initial condition
         # Plot(z_axis,r2,11)

         ##for j in range(1): ## if needed to check for multiple angle values at once, uncomment block below and change 'in range' valule. 
         # initial=[0,0.5*j+0.5]

    
             ##
             ### Calculation of image point (zi)
             ##
     halfz= 10   ## A number > 0 but small so it will search the root of marginal after index 10 in the marginal table
     zi= np.argmin(np.abs(Marginal[0][halfz:]))+halfz ## Calculate the point where the marginal cross the axis
     
     M=self.Magnification(Principal[0],Marginal[0])   ## Give Magnification: when marginal ray cross z=0
    

             ##
             ### Calculation of spherical aberration and chromatic aberration. 
             ### Cs is the spherical aberration coefficient, DXS its delta X. Cc is the chromatic aberration coefficient, DXC its delta X.
             ##
     Cs,DXS=self.Spherical_Aberration(Marginal[0][0:zi],Marginal[1][0:zi],Potential_evaluated_z_TOT[0:zi],D2_evaluated_z_TOT[0:zi],D4_evaluated_z_TOT[0:zi],z_axis[0:zi],initial) ## Calculatation of Spherical Aberration coefficient and its delta X (transverse)
    

     Cc,DXC=self.Chromatic_Aberration(Marginal[0][0:zi],Marginal[1][0:zi],Potential_evaluated_z_TOT[0:zi],z_axis[0:zi],initial,deltaE0)  ## Calculatation of Chromatic Aberration coefficient and its delta X (transverse)
    
            
             ##
             ### Calculation of spot Diameter
             ##
    
     spot=self.Spot_Diameter(DXS,DXC,M,2)
    
     return z_axis[zi],DXS,DXC,spot ## On récupère les variables qui nous sont utiles pour nos opérandes
                                    ## z_axis[zi] sera pour nous la focale
    

    def autofocus(self,focale_actuelle): ## Opérande autofocus
     return (focale_actuelle-self.Focale_Voulue)##self.Focale_Voulue est définie dans la fonction __init__ 
                                                ##pour que sa valeur soit initiée dans optim mais aussi utilisée dans d'autres fonctions



    def jacob_focus(self,F,var,h): ##Jacobienne spécifique pour autofocus => F est le vecteur d'opérandes (là juste l'autofocus) et var est le vecteur de variable
     elt=einzel_lens(show_mesh=True)##Lentille intermédiare
     n=len(var)  
     J=np.zeros(n)
     for j in range(n): ##Calculs des dérivées partielles selon chaque variable
         var_step=var.copy() ## A chaque itération on reprend les valeurs de départ de var
         var_step[j]+= h ## On fait varier la variable d'indice j
         elt.set_parameters(mesh_step=1.5, r1=1, d1=var_step[1],d2=var_step[2],t1=28,t2=26,t3=32)
         ft,S,C,spot_size=elt.simu(var_step[0],False,[0,.5],elt)                                 ## On recalcule la focale (ft) après variation de la variable
         focus=(self.autofocus(ft)) ## On obtient f(xj+h)
         ##print(focus)
         J[j]=(focus-F)/h ## Formule dérivée => [f(x+h)-f(x)]/h
     return J
    

    def DLS(self,F,p,Var,J): ## p => facteur d'ammortissement
    
     a=np.dot(np.transpose(J),J)+p*np.eye(len(Var))
     b=np.dot(-np.transpose(J),F)
     Delta_Var=np.linalg.solve(a,b)

     return Delta_Var
    

    def Ite(self,p,var):
     


     "Version pour autofocus => On ne peut pas partir de la lentille de JB => On doit utiliser une lentille avec focale différente"
     
     
     elt=einzel_lens(show_mesh=True)##Lentille intermédiare
     elt.set_parameters(mesh_step=1.5, r1=1, d1=var[1],d2=var[2],t1=28,t2=26,t3=32)##Caréctristiques de la lentilles de départ
     self.Focale_Voulue,DXS,DXC,Spot_Size=self.simu(27000,False,[0,.5])##On récupère la focale de la lentille de JB pour avoir la valeur de l'opérande autofocus
     ##elt.Focale_Voulue=self.Focale_Voulue
     ft,xst,xct,spt=elt.simu(v[0],False,[0,.5])##On récupère la focale de la lentille départ
     F_start=self.autofocus(ft)
     v0=var    ## v0 stocke nos variables de départ => On modifie v0 au lieu de var
     k=0
     Vec_f=[] ##Création tableau des vecteurs d'opérandes (pour l'instant juste autofocus)
     Vec_f.append(F_start)
     for i in range(1):## On fait tourner la fonction un certains nombres d'itérations pour obtenir plusieurs minima 
                       ## A chaque itération, on charche un nouveau minimum à partir de celui de l'itération précedente
         v=v0
         J=self.jacob_focus(F_start,v,h=1)                                             ## Cette première partie permet d'avoir un facteur d'ammortissement assez grand
         D_Var=self.DLS(F_start,p,v,J)                                                 ## pour se rapprocher d'un extremum
         v=v+D_Var
         elt.set_parameters(mesh_step=1.5, r1=1, d1=v[1],d2=v[2],t1=28,t2=26,t3=32)    ## Pour cela on effectue en boucle un DLS avec F_start en gardant le même J
         ft,xst,xct,spt=elt.simu(v[0],False,[0,.5],elt)                                ## jusqu'à que l'on est un p pour lequel (F**2)<(F_start**2) ou que p est trop grand
         F=self.autofocus(ft)
         while (F**2)>(F_start**2) and p<10000:
             k+=1
             p=10*p
             D_Var=self.DLS(F_start,p,v,J)
             v=v+D_Var
             elt.set_parameters(mesh_step=1.5, r1=1, d1=v[1],d2=v[2],t1=28,t2=26,t3=32)
             ft,xst,xct,spt=elt.simu(v[0],False,[0,.5],elt)
             F=self.autofocus(ft)
             Vec_f.append(F)
      
         if (F**2)<=(F_start**2):                                                ## Lorsque l'une des deux conditions au dessus est validée on regarde si (F**2)<=(F_start**2)
             while True:
                 p=0.1*p
                 J=self.jacob_focus(Vec_f[k],v,h=1)                              ## Si (F**2)>(F_start**2) alors cela veut dire que F_start est l'extremum que l'on cherche
                 D_Var=self.DLS(Vec_f[k],p,v,J)
                 v=v+D_Var
                 elt.set_parameters(mesh_step=1.5, r1=1, d1=v[1],d2=v[2],t1=28,t2=26,t3=32)  ## Dans le cas contraire, on réaffine p tout en effectuant des DLS entre les deux derniers vecteurs d'opérandes
                 ft,xst,xct,spt=elt.simu(v[0],False,[0,.5],elt)                              ## que l'on a stockés depuis le début dans le tableau Vec_f
                 F=self.autofocus(ft)                                                        ## On continue ainsi jusqu'à que l'on a dépassé le minimum que l'on veut (Vec_f[k-1])
                 Vec_f.append(F)
                 k+=1
                 if not ((Vec_f[k]**2)<(Vec_f[k-1]**2)):
                  break
         
             v0=v-D_Var
             p=p*10
             F_start=Vec_f[k-1]                                 ## Le nouveau minimum servira de F_start pour la prochaine itération
         print(F_start)
     print(ft)
     print(Vec_f)



if __name__ == "__main__":
    el = einzel_lens(show_mesh=True)
