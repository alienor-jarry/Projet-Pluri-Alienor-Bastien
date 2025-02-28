# Projet-Pluri-Alienor-Bastien


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
    
    def set_parameters(self, d1, d2, r1, r2=None, r3=None,
                       t1, t2=None, t3=None, rc=None, tc=None, mesh_step):
        self._d1=d1 
        self._d2=d2
        self._r1=r1
        self._t1 = t1
        if r2 == None :
            self._r2 = r1
        if r3 == None :
            self._r3 = r1
        if t2 == None :
            self._t2 = t1
        if r3 == None :
            self._t3 = t1
        if tc == None :
            self._tc = t1*1.5 # a definir
        if rc == None :
            self._rc = t1 # a definir    

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
     gmsh.fltk.run()
     gmsh.finalize()


if __name__ == "__main__":
    el = einzel_lens(show_mesh=True)
