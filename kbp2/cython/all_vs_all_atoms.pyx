# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 12:03:07 2012

@author: Tim Meyer
"""
#import pyximport; pyximport.install()
#import all_vs_all_atoms

import numpy as np
# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np

ctypedef np.float_t DTYPE_t
ctypedef np.float32_t FTYPE_t

#%cython
cdef extern from "math.h":
    double sqrt(double)


cimport cython

@cython.boundscheck(False) # turn of bounds-checking for entire function
cpdef list find_collisions(np.ndarray[double, ndim=2] atom_coord, double min_dist):
    
    cdef int i, j, d
    
    cdef int noa = len(atom_coord)
    
    cdef np.ndarray[DTYPE_t, ndim=1] coord       = np.zeros(3)
    cdef np.ndarray[DTYPE_t, ndim=1] coord_comp  = np.zeros(3)
    
    cdef float dist
    cdef float min_dist_2 = min_dist**2
    cdef np.ndarray[DTYPE_t, ndim=1] dist_v = np.zeros(3)
    
    atoms_delete_list = []
    
#    cdef float atom_coord[2,3]
#    for coor in atom_coord_np:
#        for d in range(3):
#            atom_coord[i,d] = atom_coord_np[i,d]

    for i in range( noa - 1):
        for d in range(3):
            coord[d] = atom_coord[i,d]

        #print i
        
        for j in range(i+1, noa):
            for d in range(3):
                coord_comp[d] = atom_coord[j,d]

            for d in range(3):
                dist_v[d] = coord[d] - coord_comp[d]


#            dist = np.linalg.norm(a-b)
            dist = 0
            for d in range(3):
                dist = dist + dist_v[d] * dist_v[d]
            #dist = np.sqrt(np.dot( dist_v, dist_v ))
#            dist = sqrt(dist)
            if dist < min_dist_2:
                atoms_delete_list.append((i,j))

    return atoms_delete_list



@cython.boundscheck(False) # turn of bounds-checking for entire function
@cython.cdivision(True) # turn of bounds-checking for entire function
cpdef float coulomb(np.ndarray[np.float32_t, ndim=2] atom_coord, np.ndarray[np.float32_t] charges):

    cdef int i, j, d

    cdef int noa = len(atom_coord)

    cdef np.ndarray[np.float32_t, ndim=1] coord       = np.zeros(3, dtype=np.float32)
    cdef np.ndarray[np.float32_t, ndim=1] coord_comp  = np.zeros(3, dtype=np.float32)

    cdef np.ndarray[np.float32_t, ndim=1] dist_v = np.zeros(3, dtype=np.float32)

    cdef float dist

    cdef double total_energy
    total_energy = 0.0


    cdef float charge1, coord1_x, coord1_y, coord1_z
    cdef float dx, dy, dz

#    cdef float atom_coord[2,3]
#    for coor in atom_coord_np:
#        for d in range(3):
#            atom_coord[i,d] = atom_coord_np[i,d]

    for i in range( noa - 1):
        # for d in range(3):
        #     coord[d] = atom_coord[i,d]
        coord1_x = atom_coord[i,0]
        coord1_y = atom_coord[i,1]
        coord1_z = atom_coord[i,2]
        charge1 = charges[i]

        for j in range(i+1, noa):
            # for d in range(3):
            #     coord_comp[d] = atom_coord[j,d]

            # for d in range(3):
            #     dist_v[d] = coord[d] - coord_comp[d]
                # dist_v[d] = atom_coord[i,d] - atom_coord[j,d]



            # dist = np.linalg.norm(coord-coord_comp)
            # dist = 0.0
            # for d in range(3):
                ## dist = dist + dist_v[d] * dist_v[d]
                # dist += (atom_coord[i,d] - atom_coord[j,d]) * (atom_coord[i,d] - atom_coord[j,d])
            # dist = sqrt((atom_coord[i,0] - atom_coord[j,0]) * (atom_coord[i,0] - atom_coord[j,0]) \
            #             + (atom_coord[i,1] - atom_coord[j,1]) * (atom_coord[i,1] - atom_coord[j,1]) \
            #             + (atom_coord[i,2] - atom_coord[j,2]) * (atom_coord[i,2] - atom_coord[j,2]))

            dx = coord1_x - atom_coord[j,0]
            dy = coord1_y - atom_coord[j,1]
            dz = coord1_z - atom_coord[j,2]
            dist = sqrt((dx * dx) + (dy * dy) + (dz * dz))

            # dist = np.sqrt(np.dot( dist_v, dist_v ))
            # dist = sqrt(dist)

            total_energy += charge1 * charges[j] / dist


    Na = 6.022141e23           # Avogadro constant [1/mol]
    e = 1.602176565e-19        # Elementary charge [C]
    coul = 8.987551787e9       # Coulomb constant [Mm^2/C^2]
    to_kJ_per_mol = (e**2/1e-10 * coul * Na) / 1000.0

    return total_energy * to_kJ_per_mol


@cython.boundscheck(False) # turn of bounds-checking for entire function
@cython.cdivision(True) # turn of bounds-checking for entire function
cpdef float coulomb_interaction(np.ndarray[np.float32_t, ndim=2] atom_coord1, np.ndarray[np.float32_t] charges1, \
                                 np.ndarray[np.float32_t, ndim=2] atom_coord2, np.ndarray[np.float32_t] charges2):

    cdef int i, j

    cdef int noa1 = len(atom_coord1)
    cdef int noa2 = len(atom_coord2)

    cdef float dist

    cdef double total_energy
    total_energy = 0.0


    cdef float charge1, coord1_x, coord1_y, coord1_z
    cdef float dx, dy, dz

    for i in range(noa1):
        coord1_x = atom_coord1[i,0]
        coord1_y = atom_coord1[i,1]
        coord1_z = atom_coord1[i,2]
        charge1 = charges1[i]

        for j in range(noa2):
            dx = coord1_x - atom_coord2[j,0]
            dy = coord1_y - atom_coord2[j,1]
            dz = coord1_z - atom_coord2[j,2]
            dist = sqrt((dx * dx) + (dy * dy) + (dz * dz))

            if abs(dist) > 0.01:
                total_energy += charge1 * charges2[j] / dist


    Na = 6.022141e23           # Avogadro constant [1/mol]
    e = 1.602176565e-19        # Elementary charge [C]
    coul = 8.987551787e9       # Coulomb constant [Mm^2/C^2]
    to_kJ_per_mol = (e**2/1e-10 * coul * Na) / 1000.0

    return total_energy * to_kJ_per_mol / 2.0
