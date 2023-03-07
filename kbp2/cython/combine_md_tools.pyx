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

from collections import defaultdict

import cPickle as pickle
import os, sys, time

from kbp2.kbp_tools import parse_titratable_residues

# sys.path.append('/user/tmeyer/workspace/script/protein_toolbox')
# from combine_md import parse_titratable_residues, interaction_network


from kbp2.kbp_results import *



class interaction_network(object):
    def __init__(self, nr_of_residues, ph_values=[None], nr_of_structures=1, residue_names=None):
        self.nr_of_residues = nr_of_residues
        self.ph_values = ph_values
        self.nr_of_structures = nr_of_structures
        self.residue_names = residue_names

        # Create condesed matrices.
        from scipy.misc import comb
        nr_of_interactions = int(comb(nr_of_residues, 2, exact=True))

        # Contains all interactions [structure, ph, condensed interaction matrix]
        self.data = np.zeros((nr_of_structures, len(ph_values), nr_of_interactions), dtype=np.float)


        self.last_threshold = None
        self.last_clusters = None
        self.last_interactions = None

    def __get_index(self, first_residue, second_residue):
        # res1 = first_residue
        # res2 = second_residue
        res1 = min(first_residue, second_residue)
        res2 = max(first_residue, second_residue)
        #
        # assert(first_residue != second_residue)

        # index = 0
        # for i in range(res1):
        #     index += self.nr_of_residues - i - 1
        #index += res2 - res1 - 1

        # index = (self.nr_of_residues-1) * res1
        # # The sum of all natural numbers up to (res1 - 1)
        # index -= ((res1-1)*(res1)) / 2
        # index += res2 - res1 - 1

        # index = self.nr_of_residues*res1 - res1 - (res1-1)*res1/2 + res2 - res1 - 1
        index = res1 * (self.nr_of_residues - (res1-1)/2.0 - 2) + res2 - 1

        return index

    def mod_interaction(self, first_residue, second_residue, value, ph_number=0, struct=0, set=False):
        index = self.__get_index(first_residue, second_residue)
        # index = 0

        # if type(ph_number) is not int:
        #     ph_number = self.ph_values.index(ph_number)

        # if not set:
        #     self.data[struct, ph_number, index] += value
        # else:
        #     self.data[struct, ph_number, index] = value

        self.data[struct, ph_number, index] += value

    def get_interaction(self, first_residue, second_residue, ph_number=0, struct=0):
        index = self.__get_index(first_residue, second_residue)
        return self.data[struct, ph_number, index]

    def cluster(self, threshold, ph_numbers=-1, structs=-1):
        if ph_numbers == -1:
            ph_numbers = range(len(self.ph_values))
        if type(ph_numbers) is not list:
            ph_numbers = [ph_numbers]
        if structs == -1:
            structs = range(self.nr_of_structures)
        if type(structs) is not list:
            structs = [structs]

        interactions = {}
        clusters = []

        seed_residues = range(self.nr_of_residues)
        while len(seed_residues) > 0:
            start_residues = [seed_residues[0]]
            # Create new cluster.
            cluster = []
            while len(start_residues) > 0:
                start_res = start_residues.pop(0)
                seed_residues.pop(seed_residues.index(start_res))
                cluster.append(start_res)

                # Start to search for interaction partners.
                for res in seed_residues:
                    # Residue is only added once, but the  interactions at all ph and in all structures are checked to
                    # find the strongest one.
                    for struct in structs:
                        for ph_number in ph_numbers:
                            strength = self.get_interaction(start_res, res, ph_number, struct)
                            if strength >= threshold:
                                if res not in start_residues:
                                    start_residues.append(res)
                                # Store the interaction.
                                res1 = min(start_res, res)
                                res2 = max(start_res, res)
                                if interactions.has_key(res1):
                                    if interactions[res1].has_key(res2):
                                        current_strength = interactions[res1][res2]
                                        if strength > current_strength:
                                            interactions[res1][res2] = strength
                                    else:
                                        interactions[res1][res2] = strength
                                else:
                                    interactions[res1] = {}
                                    interactions[res1][res2] = strength

            clusters.append(sorted(cluster))
        self.last_interactions = interactions

        #
        # clusters = []
        # for res1 in range(self.nr_of_residues):
        #     # Was the residue be added to an existing cluster?
        #     hit = False
        #     for struct in structs:
        #         for ph_number in ph_numbers:
        #             if type(ph_number) is not int:
        #                 ph_number = self.ph_values.index(ph_number)
        #
        #             for cluster in clusters:
        #                 for cluster_res in cluster:
        #                     if res1 == cluster_res:
        #                         continue
        #                     strength = self.get_interaction(res1, cluster_res, ph_number, struct)
        #                     if strength >= threshold:
        #                         # Add residue to an existing cluster.
        #                         cluster.append(res1)
        #                         hit = True
        #                         break
        #                 if hit:
        #                     break
        #             if hit:
        #                 break
        #         if hit:
        #             break
        #     else:
        #         # Create new cluster for the residue.
        #         clusters.append([res1])

        self.last_threshold = threshold
        self.last_clusters  = clusters

        return clusters


# cdef class c_container(object):
#     # cdef public np.ndarray[np.float_t, ndim=3] data
#     # cdef public np.ndarray[np.int_t, ndim=2] index_lookup
#     cdef public np.ndarray data
#     cdef public np.ndarray index_lookup
#     cdef public int nr_of_residues

# @cython.boundscheck(False) # turn of bounds-checking for entire function
# @cython.wraparound(False)
# cpdef mod_interaction(self, int first_residue, int second_residue, float value, int ph_number=0, int struct=0, set=False):
# #     cdef int res1 = min(first_residue, second_residue)
#     cdef int res2 = max(first_residue, second_residue)
#     # cdef nr_of_residues = self.nr_of_residues
#     cdef int index
#     # index = nr_of_residues*res1 - res1 - (res1-1)*res1/2 + res2 - res1 - 1
#     # index = res1 * (nr_of_residues - (res1-1.0)/2.0 - 2.0) + res2 - 1
#     index = self.nr_of_residues*res1 - (res1+1)*res1/2 + res2 - res1 - 1
#     # print "%i | %i -> %s" % (res1, res2, index)
#
#     # cdef np.ndarray[np.int_t, ndim=2] index_lookup
#     # index_lookup = self.index_lookup
#     # index = self.index_lookup[res1, res2]
#
#     # if type(ph_number) is not int:
#     #     ph_number = self.ph_values.index(ph_number)
#
#     # if not set:
#     #     self.data[struct, ph_number, index] += value
#     # else:
#     #     self.data[struct, ph_number, index] = value
#
#     cdef np.ndarray[np.float_t, ndim=3] data
#     data = self.data
#     data[struct, ph_number, index] += value
#
# class c_interaction_network(object):
#     def __init__(self, nr_of_residues, ph_values=[None], nr_of_structures=1, residue_names=None):
#         self.nr_of_residues = nr_of_residues
#         self.ph_values = ph_values
#         self.nr_of_structures = nr_of_structures
#         self.residue_names = residue_names
#
#         # Create condesed matrices.
#         from scipy.misc import comb
#         nr_of_interactions = int(comb(nr_of_residues, 2, exact=True))
#
#         # Contains all interactions [structure, ph, condensed interaction matrix]
#         self.data = np.zeros((nr_of_structures, len(ph_values), nr_of_interactions), dtype=np.float)
#
#         # self.index_lookup = np.zeros((nr_of_residues, nr_of_residues), dtype=np.int)
#         # for i in range(nr_of_residues):
#         #     for j in range(i+1, nr_of_residues):
#         #         self.index_lookup[i, j] = nr_of_residues*i - (i+1)*i/2 + j - i - 1
#
#         # self.c_container = c_container()
#         # self.c_container.data = self.data
#         # self.c_container.index_lookup = self.index_lookup
#         # self.c_container.nr_of_residues = nr_of_residues
#
#         # self.mod_interaction = self.c_container.mod_interaction
#         # self.mod_interaction = mod_interaction
#
#
#
#     def get_interaction(self, first_residue, second_residue, ph_number=0, struct=0):
#         index = self.__get_index(first_residue, second_residue)
#         return self.data[struct, ph_number, index]
#
#     def cluster(self, threshold, ph_numbers=[0], structs=[0]):
#         if type(structs) is not list:
#             structs = [structs]
#         if type(ph_numbers) is not list:
#             ph_numbers = [ph_numbers]
#
#         clusters = []
#         for res1 in range(self.nr_of_residues):
#             for struct in structs:
#                 for ph_number in ph_numbers:
#                     if type(ph_number) is not int:
#                         ph_number = self.ph_values.index(ph_number)
#
#                     # Should the residue be added to an existing cluster?
#                     hit = False
#                     for cluster in clusters:
#                         for cluster_res in cluster:
#                             if res1 == cluster_res:
#                                 continue
#                             strength = self.get_interaction(res1, cluster_res, ph_number, struct)
#                             if strength >= threshold:
#                                 # Add residue to an existing cluster.
#                                 cluster.append(res1)
#                                 hit = True
#                                 break
#                         if hit:
#                             break
#                     if hit:
#                         break
#                 if hit:
#                     break
#             else:
#                 # Create new cluster for the residue.
#                 clusters.append([res1])
#
#         return clusters


# class KbpJobDescr(object):
#     """Object that contains information about a karlsberg+ job. Does not contain results.
#     """
#     def __init__(self):
#         # List of titrated residues, ordered by segment name and resid.
#         self.sorted_residue_list = None
#
#         # List of titrated residues, based on 'sorted_residue_list' but with additional information to allow efficient
#         # iteration over the resides. Each entry contains a tuple:
#         # (residue_num, residue_kbp, resname_kbp, nr_of_states)
#         # with: residue_num : Index of the residue in 'sorted_residue_list'
#         #       residue_kbp : Residue description within the framework of Karlsberg+ (e.g. 'ARG-12_A')
#         #       resname_kbp : Residue name within the framework of Karlsberg+ (e.g. ARG or EPP)
#         self.residue_list_ext = None
#
#         # Sorted list of strings containing pH values used for the titration.
#         self.ph_values = None
#
#         # Contains information about the type of resides that where titratable in the Karlsberg+ calculation and their
#         # properties. Format:
#         # {<Karlsberg+ residue name> : [ state : {'name'  : str <state name e.g. 'R', 'P', 'D', 'PP'>
#         #                                         'pka'   : float <energy of the residue at pH 0 in units of pKa.>
#         #                                         'patch' : str <name of the CHARMM patch that creates the state.>
#         #                                         'atoms' : dict {<atom name> : <charge>}
#         #                                       }
#         #                             ]
#         #         }
#         self.titratable_residues = None


# class KbpResult(object):
    # """Container to store the results of a Karlsberg+ calculation.
    # """
    #
    # def __init__(self, descr=None):
    #     """Constructor
    #     Parameters: descr : Object of type KbpJobDescr. The object is used directly, not copied!
    #     """
    #
    #     # Contains general properties of the Karlsberg+ job.
    #     if descr is None:
    #         self.descr = KbpJobDescr()
    #     else:
    #         assert(isinstance(descr, KbpJobDescr))
    #         self.descr = descr
    #
    #
    #     ### Karlsberg+ results ###
    #     # Intrinsic pKa values. (These are the pKa values for each residue, if all other residues would be in the
    #     # reference ('R') state. Format:
    #     # {str <Karlsberg+ residue description e.g. ARG-12_A> : numpy array [state : float <intrinsic pKa>]
    #     # }
    #     self.pkaint = None
    #
    #     # Numpy matrix that contains the interaction energy between all residues for all possible states.
    #     # Format:
    #     # 4-dimensional numpy array: [residue2_num, residue_num, state2, state]
    #     # with: residue_num  : Number of the first residue. This is the index of this residue in 'sorted_residue_list'.
    #     #       residue2_num : Number of the second residue. This is the index of this residue in 'sorted_residue_list'.
    #     #       state  : State of first residue. This is the index of the state in 'titratable_residues'.
    #     #       state2 : State of second residue. This is the index of the state in 'titratable_residues'.
    #     self.g = None
    #
    #     # Numpy matrix that contains the occupancy of all states of all residues for each pH.
    #     # Format:
    #     # 3-dimensional numpy array: [residue_num, state, ph_number]
    #     # with: residue_num  : Number of the first residue. This is the index of this residue in 'sorted_residue_list'.
    #     #       state        : State of first residue. This is the index of the state in 'titratable_residues'.
    #     #       ph_number    : Number of the pH value. This is the index of the pH in 'ph_values'.
    #     self.occs = None
    #
    #     # Contains the occupancy of each 'pH adopted conformation' (PAC).
    #     # Format: [int <ph_number> : float <occupancy of the PAC>]
    #     # with: pH number : Number of the pH value. This is the index of the pH in 'ph_values'.
    #     self.conformer_occs = None


# class FrameKbpResults(object):
#     """Object to store the results of Karlsberg+ calculations that are based on MD snapshots.
#     """
#
#     def __init__(self, descr=None):
#         """ Constructor
#         Parameters: descr : Object of type KbpJobDescr. The object is used directly, not copied!
#         """
#
#         # Contains general properites of the Kalrsberg+ job.
#         if descr is None:
#             self.descr = KbpJobDescr()
#         else:
#             assert(isinstance(descr, KbpJobDescr))
#             self.descr = descr
#
#         # Contains 'kbp_results' objects, each representing the results of the Karlsberg+ calculation for a MD snapshot.
#         # Format: {<frame number> : kbp_results <Karlsberg+ results>}
#         self.kbp_results = {}
#
#         # Float the translates frame number into nanoseconds. (optional)
#         self.timestep = None
#
#         # Specifies the frame range that is supposed to be analysed. Some functions use this variable, if set. (optional)
#         # Format: [<first frame number>, <last frame number>]
#         # -1 can be used for <last frame number> to select the last frame.
#         self.frame_range = [0, -1]
#
#     def add_job(self, kbp_result, frame_nr=None):
#         """ Adds a KbpResults object. If the property kbp_results.desc is None is is set to the self.descr.
#         """
#         if frame_nr is None:
#             frame_nr = len(self.kbp_results)
#
#         assert(not self.kbp_results.has_key(frame_nr))
#         assert(isinstance(kbp_result, KbpResult))
#
#         if kbp_result.descr is None:
#             kbp_result.descr = self.descr
#
#         self.kbp_results[frame_nr] = kbp_result
#
#     def is_in_range(self, frame_nr, frame_range=None):
#         """Checks if a given frame_nr is within a interval defined by frame_range. it returns the result of:
#         frame_range[0] < frame_nr < frame_range[1]
#         If frame_range[1] has the value -1, no upper boundary is considered.
#         @param frame_nr: int
#         @param frame_range: [<first frame number>, <last first frame number or -1>] if not specified, it is taken from
#                             self.frame_range
#         @return: bool
#         """
#         if frame_range is None:
#             frame_range = self.frame_range
#
#         if frame_range is not None:
#             x_min = frame_range[0]
#             x_max = frame_range[1]
#             if frame_nr < x_min:
#                 return False
#             if x_max != -1:
#                 if frame_nr > x_max:
#                     return False
#         return True
#
#     # Todo: average self.conformer_occs
#     def averge_over_frames(self, frame_range=None):
#         """Averages the following properties of the stored frames_kbp_results object over the frames: pkaint, g and occs
#         The frame selection of self.frame_range is used if given as parameter or set previously. The results are stored
#         as: self.pkaint, self.g and self.occs
#
#         @return: None
#         """
#         if frame_range is not None:
#             self.frame_range = list(frame_range)
#
#         nr_of_residues = len(self.descr.sorted_residue_list)
#         nr_of_ph_values = len(self.descr.ph_values)
#         max_states = 0
#         for residue_kbp in self.descr.sorted_residue_list:
#             resname_kbp = residue_kbp[0:residue_kbp.find('-')]
#             # Skip the R state.
#             nr_of_states = len(self.descr.titratable_residues[resname_kbp])
#             if max_states < nr_of_states:
#                 max_states = nr_of_states
#         self.g    = np.zeros((nr_of_residues, nr_of_residues, max_states, max_states), dtype=np.float)
#         self.occs = np.zeros((nr_of_residues, max_states, nr_of_ph_values), dtype=np.float)
#         self.pkaint = {}
#
#         nr_of_frames = 0
#         for frame_nr, kbp_result in self.kbp_results.iteritems():
#             if not self.is_in_range(frame_nr):
#                 continue
#
#             nr_of_frames += 1
#
#             self.g    += kbp_result.g
#             self.occs += kbp_result.occs
#             if not self.pkaint:
#                 for residue in kbp_result.pkaint.keys():
#                     self.pkaint[residue] = kbp_result.pkaint[residue]
#             else:
#                 for residue in kbp_result.pkaint.keys():
#                     self.pkaint[residue] += kbp_result.pkaint[residue]
#
#         self.g    /= nr_of_frames
#         self.occs /= nr_of_frames
#         for residue in self.pkaint.keys():
#             self.pkaint[residue] /= nr_of_frames

# class KbpTools(object):
def interactive_plot_vs_ph(descr, residue_property, title='', show=True):
    """Plots a residue specific property against the ph values.

    @param descr: KbpJobDescr object
    @param residue_property: The residues properties: [residue : [ph : <value to plot>]]
                             Each entry in the residue list corresponds to the entry with the same index in
                             descr.sorted_residue_list, each entry in the ph list corresponds to the entry with the
                             same index in descr.ph_values.
    @param title: str: Title shown on top of the graph.
    @param show: bool: plt.show(True) is called at the end of the function, if set to True.
    @return:
    """
    residue_list_ext = descr.residue_list_ext
    nr_of_residues = len(descr.sorted_residue_list)
    ph_values = descr.ph_values

    from matplotlib import pyplot as plt

    colors = plt.cm.Paired(np.linspace(0, 1, nr_of_residues))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.ion()
    line_names = {}
    lines = []

    residue_list_filtered = []
    for (residue_num, residue_kbp, resname_kbp, nr_of_states) in residue_list_ext:
        if resname_kbp not in ['EPP', 'DPP', 'HSP']:
            continue
        line, = ax.plot(ph_values, residue_property[residue_num], 'x-', color=colors[residue_num], picker=5)
        line_names[line] = residue_kbp
        lines.append(line)
        residue_list_filtered.append(residue_kbp)
    leg = ax.legend(residue_list_filtered, loc='best', prop={'size': 7}, ncol=4)
    leg.draw_frame(False)
    plt.title(title)
    ymin, ymax = plt.ylim()
    plt.ylim([ymin - 0.1, 1.0])

    lined = dict()
    for legline, origline in zip(leg.get_lines(), lines):
        legline.set_picker(5)  # 5 pts tolerance
        lined[legline] = origline

    onpick = OnPick(title, line_names, lined, lines, ax, plt, residue_list_filtered)
    # from matplotlib.artist import setp
    # def onpick(event):
    #     for line in line_names.keys():
    #         setp(line, linewidth=1, marker='x')
    #
    #     # lined[legline]
    #     if line_names.has_key(event.artist):
    #         thisLine = event.artist
    #         residue = line_names[thisLine]
    #     elif lined.has_key(event.artist):
    #         thisLine = lined[event.artist]
    #         residue = line_names[thisLine]
    #
    #     setp(thisLine, linewidth=2, marker='o')
    #     # event.artist.set_label('---')
    #
    #     # leg = ax.legend(residues, loc='best', prop={'size': 7}, ncol=4)
    #     leg = ax.legend(residue_list_filtered, loc='best', prop={'size': 7}, ncol=4)
    #     leg.draw_frame(False)
    #     # lined = dict()
    #     for legline, origline in zip(leg.get_lines(), lines):
    #         legline.set_picker(5)  # 5 pts tolerance
    #         lined[legline] = origline
    #
    #     if title is not None:
    #         plt.title(title + ' - ' +  residue)
    #     else:
    #         plt.title(residue)
    fig.canvas.mpl_connect('pick_event', onpick)

    if show:
        plt.show(True)

class OnPick(object):
    def __init__(self, title, line_names, lined, lines, ax, plt, residue_list_filtered):
        self.line_names = line_names
        self.lined = lined
        self.ax = ax
        self.plt = plt
        self.residue_list_filtered = residue_list_filtered
        self.title = title
        self.lines = lines
    def __call__(self, event):
        from matplotlib.artist import setp

        for line in self.line_names.keys():
            setp(line, linewidth=1, marker='x')

        # lined[legline]
        if self.line_names.has_key(event.artist):
            thisLine = event.artist
            residue = self.line_names[thisLine]
        elif self.lined.has_key(event.artist):
            thisLine = self.lined[event.artist]
            residue = self.line_names[thisLine]

        setp(thisLine, linewidth=2, marker='o')
        # event.artist.set_label('---')

        # leg = ax.legend(residues, loc='best', prop={'size': 7}, ncol=4)
        leg = self.ax.legend(self.residue_list_filtered, loc='best', prop={'size': 7}, ncol=4)
        leg.draw_frame(False)
        # self.lined = dict()
        for legline, origline in zip(leg.get_lines(), self.lines):
            legline.set_picker(5)  # 5 pts tolerance
            self.lined[legline] = origline

        if self.title is not None:
            self.plt.title(self.title + ' - ' +  residue)
        else:
            self.plt.title(residue)



@cython.boundscheck(False) # turn of bounds-checking for entire function
cpdef get_prot_energy(descr, pkaint, np.ndarray[np.float_t, ndim=4] g, np.ndarray[np.float_t, ndim=3] occs,
                      ref_residues=None, ref_ph=None):
    # Todo: write docstring
    """

    @param descr:
    @param pkaint:
    @param g:
    @param occs:
    @return:
    """
    cdef float occ, occ2
    cdef float RT_ln10_ph
    cdef float protE_sum, protE_res_sum
    cdef int state, state2
    cdef int nr_of_states, nr_of_states2
    cdef int residue_num, residue2_num
    cdef int ph_number
    cdef float dG_interact
    # cdef float dG_min, dG_max

    R = 8.3144621                             # J/(K*mol)
    # kb = 1.3806488e-23                      # J/K
    # RT = R * 300.0 / 1000.0                 # kJ/mol
    RT_ln10 = R * 300.0 / 1000.0 * np.log(10) # kJ/mol

    ph_values = descr.ph_values
    residue_list_ext = descr.residue_list_ext
    titratable_residues = descr.titratable_residues
    nr_of_residues = len(descr.sorted_residue_list)


    if ref_residues is None:
        residue_list_ext_ref = residue_list_ext
    else:
        residue_list_ext_ref = []
        for (residue_num, residue_kbp, resname_kbp, nr_of_states) in residue_list_ext:
            if residue_kbp in ref_residues:
                residue_list_ext_ref.append( (residue_num, residue_kbp, resname_kbp, nr_of_states) )

    protE     = {}
    protE_res = {}
    for ph_number, ph in enumerate(ph_values):
        if ref_ph is not None:
            if ph != ref_ph:
                continue
        protE_res[ph] = {}
        protE_sum = 0.0

        RT_ln10_ph = RT_ln10 * ph

        # for (residue_num, residue_kbp, resname_kbp, nr_of_states) in residue_list_ext:
        #     for (residue2_num, residue2_kbp, resname2_kbp, nr_of_states2) in residue_list_ext[residue_num+1:]:
        #         dG_min =  np.inf
        #         dG_max = -np.inf
        #         for state in range(1, nr_of_states):
        #             for state2 in  range(1, nr_of_states2):
        #                 dG_interact = g[residue_num, residue2_num, state, state2]
        #                 if dG_interact > dG_max:
        #                     dG_max = dG_interact
        #                 if dG_interact < dG_min:
        #                     dG_min = dG_interact
        #         index = nr_of_residues*residue_num - \
        #                 ((residue_num+1)*residue_num)/2 + residue2_num - residue_num - 1
        #         all_states_interactions[md, ph_number, index] = dG_max - dG_min


        for (residue_num, residue_kbp, resname_kbp, nr_of_states) in residue_list_ext_ref:
            protE_res_sum = 0.0

            #####################
            ### INTRINSIC PKA ###
            #####################
            for state in range(1, nr_of_states):
                occ = occs[residue_num, state, ph_number]

                state_name = titratable_residues[resname_kbp][state]['name']
                state_pka  = titratable_residues[resname_kbp][state]['pka']
                state_pka *= RT_ln10
                # slope_sign ->  0 for R/0 state
                #            ->  1 for P   state
                #            -> -1 for D   state
                pka_int = pkaint[residue_kbp][state]
                delta_pka = pka_int - state_pka

                # DEBUG #
                # delta_pka = 0.0

                if state_name == 'P':
                    # dG = (pka_int + 1.0 * RT_ln10_ph) * occ
                    dG = (delta_pka + (state_pka + 1.0 * RT_ln10_ph)) * occ
                elif state_name == 'D':
                    # dG = (pka_int - 1.0 * RT_ln10_ph) * occ
                    dG = (delta_pka + (state_pka - 1.0 * RT_ln10_ph)) * occ
                elif state_name == '0':
                    # dG = pka_int * occ
                    dG = (delta_pka + state_pka) * occ
                else:
                    assert("Sorry, state type %s is not yet implemented." % state_name)

                protE_sum     += dG
                protE_res_sum += dG

                ###########################
                ### Interaction energy. ###
                ###########################
                # Iterate over all other residues.
                # if state_type != 'R' and occ > 0.001:
                if  occ > 0.001:
                    for (residue2_num, residue2_kbp, resname2_kbp, nr_of_states2) in residue_list_ext:

                        if residue_num == residue2_num:
                            continue

                        for state2 in  range(1, nr_of_states2):
                            occ2 = occs[residue2_num, state2, ph_number]

                            ###################
                            ### INTERACTION ###
                            ###################
                            if residue2_num > residue_num:
                                dG_interact = g[residue_num, residue2_num, state, state2] * occ * occ2
                                protE_sum += dG_interact
                            else:
                                dG_interact = g[residue2_num, residue_num, state2, state] * occ * occ2
                            protE_res_sum += dG_interact

            protE_res[ph][residue_kbp] = protE_res_sum
        protE[ph] = protE_sum

    return protE, protE_res


cpdef analyse_tension3(kbp_result, sasa_factor):
    assert(isinstance(kbp_result, FrameKbpResults))

    ph_values = kbp_result.descr.ph_values
    residue_list_ext = kbp_result.descr.residue_list_ext
    md_weights = kbp_result.weights
    titratable_residues = kbp_result.descr.titratable_residues

    for ph_number, ph in enumerate(ph_values):
        for (residue_num, residue_kbp, resname_kbp, nr_of_states) in residue_list_ext:
            # Get current energy of the residue and its protonation state.

            occs = kbp_result.combined_results.occs[residue_num][:,ph_number]
            current_state = np.argmax(occs)
            current_state_name = titratable_residues[resname_kbp][current_state]['name']
            current_state_name.replace('R', '0')
            current_md = np.argmax(md_weights[ph_number])

            current_res_prot_energy = 0.0
            current_res_sasa_corr = 0.0
            current_res_confE = 0.0
            for md_number, md_kbp_result in kbp_result:
                current_res_prot_energy += md_weights[ph_number][md_number] * kbp_result.kbp_results[md_number].residue_prot_energy[residue_kbp][ph_number]
                current_res_sasa_corr += md_weights[ph_number][md_number] * sasa_factor * kbp_result.kbp_results[md_number].residue_sasa[residue_num]
                current_res_confE += md_weights[ph_number][md_number] * np.sum(kbp_result.kbp_results[md_number].conf_energies['residue_energy'][residue_num])

            # Compare with the energy of the individual MDs
            # If protonation state differs, suggest this protonation state for the next MD
            best_alternative_state = None
            best_alternative_state_name = None
            best_alternative_delta_e = None
            for md_number, md_kbp_result in kbp_result:
                occs = md_kbp_result.occs[residue_num][:,ph_number]
                md_state = np.argmax(occs)
                md_state_name = titratable_residues[resname_kbp][md_state]['name']
                md_state_name.replace('R', '0')

                res_prot_energy = md_kbp_result.residue_prot_energy[residue_kbp][ph_number]
                res_sasa_corr = sasa_factor * md_kbp_result.residue_sasa[residue_num]
                res_confE = np.sum(md_kbp_result.conf_energies['residue_energy'][residue_num])

                if current_state != md_state:
                    if (current_state_name != md_state_name) or (residue_kbp in ['HSP']):
                        delta_e = res_prot_energy - current_res_prot_energy
                        # delta_e += res_sasa_corr - current_res_sasa_corr
                        delta_e += res_confE - current_res_confE
                        if delta_e < -20.0:
                            if (best_alternative_delta_e is None) or (delta_e < best_alternative_delta_e):
                                best_alternative_state = md_state
                                best_alternative_state_name = md_state_name
                                best_alternative_delta_e = delta_e

            if best_alternative_state is not None:
                print "%.2f - %s  %i (%s) -> %i (%s): %.2f" % (ph, residue_kbp, current_state, current_state_name, best_alternative_state, best_alternative_state_name, best_alternative_delta_e)








cpdef analyse_tension2(md, md_kbp_results, confE, md_protE_res, md_weight, md_res_confE, md_res_near_sites, md_res_wat_energy, frame_range):
    """ Finds tensions in the results of a MD-pKa run.

    @param md:
    @param md_kbp_results:
    @param md_protE_res:
    @param md_weight:
    @param md_res_confE:
    @param md_res_near_sites:
    @param frame_range:
    @return:
    """

    threshold_for_state_occ = 0.2
    threshold_for_md_occ = 0.2

    R = 8.3144621                             # J/(K*mol)
    RT_ln10 = R * 300.0 / 1000.0 * np.log(10) # kJ/mol


    # Calculate tension for each pH value.
    ref_kbp_results = md_kbp_results[0]
    ph_values        = ref_kbp_results.descr.ph_values
    residue_list_ext = ref_kbp_results.descr.residue_list_ext
    nr_of_mds = len(md_kbp_results)
    titratable_residues = ref_kbp_results.descr.titratable_residues

    # Average all required values over the frames.
    for frame_kbp_results in md_kbp_results:
        frame_kbp_results.averge_over_frames()

    # Average res_confE for each MD over all frames.
    avg_res_confE_md = []
    for res_confE in md_res_confE:
        # print res_confE
        avg_res_confE = np.average(res_confE.values(), 0)
        # print avg_res_confE[5]
        # print np.std(res_confE.values(), 0)
        avg_res_confE_md.append(avg_res_confE)


    # Average confE over all frames.
    avg_confE = []
    for confE_md in confE:
        avg_confE.append(np.average(confE_md.values()))

    # DEBUG #
    # for avg_res_confE in avg_res_confE_md:
    #     avg_res_confE *= 0.0

    # DEBUG #
    # for kbp_results in md_kbp_results:
    #     kbp_results.g *= 0.0
    # for ph_number, ph in enumerate(ph_values):
    #     if ph > 1.0:
    #         md_weight[ph_number] = [1.0, 0.0]
    #     else:
    #         md_weight[ph_number] = [0.0, 1.0]

    res_ph_tension = defaultdict(dict)
    for ph_number, ph in enumerate(ph_values):
        # -> current_res_enery
        # -> Für jede andere MD mit abweichender protonierung -> md_res_energy

        for (residue_num, residue_kbp, resname_kbp, nr_of_states) in residue_list_ext:
            # if residue_kbp != 'EPP-38_A':
            #     continue
            # set fixed:     ph / residue
            # averaged over: frames
            # variable:      md / state
            # 1) Determine MDs with deviating protonation state of the current residue.
            best_md = 0
            best_md_occ = 0.0
            for md_nr in range(nr_of_mds):
                if md_weight[ph_number][md_nr] > best_md_occ:
                    best_md = md_nr
                    best_md_occ = md_weight[ph_number][md_nr]
            best_state = 0
            best_state_occ = 0.0
            for state in  range(nr_of_states):
                occ = md_kbp_results[best_md].occs[residue_num, state, ph_number]
                if occ > best_state_occ:
                    best_state = state
                    best_state_occ = occ
            current_state = best_state
            current_md = best_md
            current_state_type = titratable_residues[resname_kbp][current_state]['name']

            # 2) Get current_res_energy for the current situation
            current_protE = 0.0
            current_protE_res = 0.0
            current_confE = 0.0
            current_occs = md_kbp_results[md_nr].occs.copy() * 0.0
            # print "---"
            # print ph
            for md_nr in range(nr_of_mds):
                descr = md_kbp_results[md_nr].descr
                pkaint = md_kbp_results[md_nr].pkaint
                g = md_kbp_results[md_nr].g
                occs = md_kbp_results[md_nr].occs
                protE, protE_res = get_prot_energy(descr, pkaint, g, occs, ref_residues=[residue_kbp], ref_ph=ph)
                # current_protE_res += protE[ph] * md_weight[ph_number][md_nr]
                current_protE_res += protE_res[ph][residue_kbp] * md_weight[ph_number][md_nr]
                # current_confE += avg_confE[md_nr] * md_weight[ph_number][md_nr]
                current_confE += np.sum(avg_res_confE_md[md_nr][residue_num]) * md_weight[ph_number][md_nr]
                current_occs += occs * md_weight[ph_number][md_nr]

            current_res_enery = current_protE_res + current_confE + md_res_wat_energy[md_nr][residue_kbp]

            # 3) Find other MDs that other states are occupied
            #    Apply current protonation states to MDs from 3) except for the current residue.
            #    -> get md_res_energy
            res_ph_tension[residue_kbp][ph] = 0.0
            for md_nr in range(nr_of_mds):
                if md_nr == current_md:
                    continue

                best_state = 0
                best_state_occ = 0.0
                for state in range(nr_of_states):
                    occ = md_kbp_results[md_nr].occs[residue_num, state, ph_number]
                    if occ > best_state_occ:
                        best_state = state
                        best_state_occ = occ
                # if best_state == current_state:
                #     continue


                descr = md_kbp_results[md_nr].descr
                pkaint = md_kbp_results[md_nr].pkaint
                g = md_kbp_results[md_nr].g
                # g = np.minimum(md_kbp_results[md_nr].g, md_kbp_results[best_md].g)
                # all_g = [x.g for x in md_kbp_results]
                # g = reduce(np.minimum, all_g)

                occs = md_kbp_results[md_nr].occs

                # Try all states.
                # for state in range(1, nr_of_states):
                state = best_state
                if True:
                    state_type = titratable_residues[resname_kbp][state]['name']
                    # The following lines are required if the R state is relevant.
                    # if state_type == 'R':
                    #     state_type = '0'
                    if state_type == current_state_type:
                        continue
                    # if state == current_state:
                    #     continue

                    # current_occs[residue_num, :, ph_number] = 0.0
                    # current_occs[residue_num, state, ph_number] = 1.0
                    # protE, protE_res = get_prot_energy(descr, pkaint, g, current_occs, ref_residues=[residue_kbp],
                    #                                    ref_ph=ph)
                    protE, protE_res = get_prot_energy(descr, pkaint, g, occs, ref_residues=[residue_kbp],
                                                       ref_ph=ph)
                    # Todo: store all promising results

                    # res_energy = protE[ph] + avg_confE[md_nr]
                    # res_energy = protE_res[ph][residue_kbp] + avg_confE[md_nr]
                    res_energy = protE_res[ph][residue_kbp] + np.sum(avg_res_confE_md[md_nr][residue_num])  + md_res_wat_energy[md_nr][residue_kbp]
                    dg_res_energy = res_energy - current_res_enery
                    if dg_res_energy < res_ph_tension[residue_kbp][ph]:
                        res_ph_tension[residue_kbp][ph] = dg_res_energy



    descr = md_kbp_results[0].descr
    res_ph_tension_list = []
    for (residue_num, residue_kbp, resname_kbp, nr_of_states) in residue_list_ext:
        res_ph_tension_list.append([])
        for ph in ph_values:
            res_ph_tension_list[-1].append(res_ph_tension[residue_kbp][ph])
    interactive_plot_vs_ph(descr, res_ph_tension_list, 'tension')

    # from matplotlib import pyplot as plt
    # plt.figure()
    # legend = []
    # for residue in res_ph_tension.keys():
    #     x = sorted(res_ph_tension[residue].keys())
    #     y = [res_ph_tension[residue][ph] for ph in x]
    #     plt.plot(x, y, '-x')
    #     legend.append(residue)
    # plt.legend(legend)
    # plt.show()
    #
    # sys.exit()


@cython.boundscheck(False) # turn of bounds-checking for entire function
cpdef calc_prot_energy(self, frame_range, restart=False):

    cdef float occ, occ2
    cdef float dG_interact
    cdef float RT_ln10_ph
    cdef float protE_sum, protE_res_sum
    cdef float value
    cdef int state, state2
    cdef int nr_of_states, nr_of_states2
    cdef int residue_num, residue2_num
    cdef int frame_nr
    cdef int ph_number
    cdef int md
    cdef np.ndarray[np.float_t, ndim=4] g_frame_np
    cdef np.ndarray[np.float_t, ndim=3] occs_np
    cdef float dG_min, dG_max

    cdef int index, nr_of_residues

    if self.kbp_results is None:
        assert("self.kbp_results has to be set with 'get_titration_curves' before self.calc_prot_energy can be called.")

    R = 8.3144621                             # J/(K*mol)
    # kb = 1.3806488e-23                      # J/K
    # RT = R * 300.0 / 1000.0                 # kJ/mol
    RT_ln10 = R * 300.0 / 1000.0 * np.log(10) # kJ/mol

    # The titratable.yaml file is read only once for all MDs.
    titratable_residues = None


    # residue_list = self.md_g_pkint[0].values()[0]['residue_list']
    residue_list = self.sorted_residue_list
    nr_of_residues = len(residue_list)
    ph_values = self.kbp_results[0].values()[0]['ph_values']
    nr_of_mds = len(self.jobfolders)
    inet = interaction_network(nr_of_residues, ph_values, nr_of_mds, residue_list)

    from scipy.misc import comb
    nr_of_interactions = int(comb(nr_of_residues, 2, exact=True))
    cdef np.ndarray[np.float_t, ndim=3] all_states_interactions
    all_states_interactions = np.zeros((nr_of_mds, len(ph_values), nr_of_interactions), dtype=np.float)
    inet.data = all_states_interactions


    md_protE = []
    md_protE_res = []

    md_kbp_results = []

    for md, jobfolder in enumerate(self.jobfolders):
        if jobfolder[-1] != '/':
            jobfolder += '/'

        pickle_filename = jobfolder + 'prot_energies_detailed%s.pickle' % self.kbp_suffix
        if os.path.exists(pickle_filename) and not restart:
            print("Reading prot_energies_detailed cPickle file for: %s" % jobfolder)

            f = open(pickle_filename, 'rb')
            (protE, protE_res, frame_kbp_results) = pickle.load(f)
            f.close

        else:
            print "Calculating protonation energies for: " + jobfolder

            kbp_folder = jobfolder + 'kbp' + self.kbp_suffix + '/done/'

            if titratable_residues is None:
                titratable_residues = parse_titratable_residues(kbp_folder)

            frame_kbp_results = FrameKbpResults()
            frame_kbp_results.descr.ph_values = list(ph_values)
            frame_kbp_results.descr.sorted_residue_list = list(residue_list)
            frame_kbp_results.descr.titratable_residues = titratable_residues
            # The next entry is set later.
            frame_kbp_results.descr.residue_list_ext = None
            frame_kbp_results.frame_range = list(frame_range)


            protE     = {}
            protE_res = {}
            for frame_nr in self.md_g_pkint[md].keys():
                # start = time.time()
                if not self.is_in_range(frame_nr, frame_range):
                    continue

                # residue_list_ref = ['DPP-21_A']
                # residue_list_ref = []
                # residue_list_ref = ['EPP-48_A']
                # residue_list_ref = ['HSP-83_A', 'LYS-86_A', 'DPP-108_A']
                # residue_list_ref = ['DPP-10_A', 'EPP-48_A', 'DPP-70_A', 'DPP-134_A']
                # residue_list_ref = ['DPP-66_A', 'EPP-38_A']

                g_frame           = self.md_g_pkint[md][frame_nr]['g']
                pkaint_frame      = self.md_g_pkint[md][frame_nr]['pkint']
                kbp_results_frame = self.kbp_results[md][frame_nr]['curves']

                ph_values = self.kbp_results[0].values()[0]['ph_values']

                residue_list_ext = []
                residue_list_ext_ref = []
                max_states = 0
                for residue_num, residue_kbp in enumerate(residue_list):
                    resname_kbp = residue_kbp[0:residue_kbp.find('-')]
                    # Skip the R state.
                    nr_of_states = len(titratable_residues[resname_kbp])
                    if max_states < nr_of_states:
                        max_states = nr_of_states
                    # states       = range(1, len(titratable_residues[resname_kbp]))
                    # state_names  = [state_entry['name'] for state_entry in titratable_residues[resname_kbp][1:]]
                    residue_list_ext.append( (residue_num, residue_kbp, resname_kbp, nr_of_states) )
                    # if residue_kbp in residue_list_ref:
                    #     residue_list_ext_ref.append( (residue_num, residue_kbp, resname_kbp, nr_of_states) )
                residue_list_ext_ref = residue_list_ext

                frame_kbp_results.descr.residue_list_ext = list(residue_list_ext)

                nr_of_residues = len(residue_list)
                g_frame_np = np.zeros((nr_of_residues, nr_of_residues, max_states, max_states), dtype=np.float)
                nr_of_ph_values = len(ph_values)
                occs_np = np.zeros((nr_of_residues, max_states, nr_of_ph_values), dtype=np.float)
                for (residue_num, residue_kbp, resname_kbp, nr_of_states) in residue_list_ext:
                    for state in range(1, nr_of_states):
                        for (residue2_num, residue2_kbp, resname2_kbp, nr_of_states2) in residue_list_ext[residue_num+1:]:
                            for state2 in range(1, nr_of_states2):
                                value = g_frame[residue_kbp, residue2_kbp, state, state2]
                                g_frame_np[residue_num, residue2_num, state, state2] = value

                        for ph_number in range(len(ph_values)):
                            value = kbp_results_frame[residue_kbp]['states'][state][ph_number]
                            occs_np[residue_num, state, ph_number] = value


                protE[frame_nr]     = {}
                protE_res[frame_nr] = {}
                for ph_number, ph in enumerate(ph_values):
                    protE[frame_nr][ph] = 0.0
                    protE_sum = 0.0
                    protE_res[frame_nr][ph] = {}

                    RT_ln10_ph = RT_ln10 * ph

                    for (residue_num, residue_kbp, resname_kbp, nr_of_states) in residue_list_ext:
                        for (residue2_num, residue2_kbp, resname2_kbp, nr_of_states2) in residue_list_ext[residue_num+1:]:
                            dG_min =  np.inf
                            dG_max = -np.inf
                            for state in range(1, nr_of_states):
                                for state2 in  range(1, nr_of_states2):
                                    dG_interact = g_frame_np[residue_num, residue2_num, state, state2]
                                    if dG_interact > dG_max:
                                        dG_max = dG_interact
                                    if dG_interact < dG_min:
                                        dG_min = dG_interact
                            index = nr_of_residues*residue_num - \
                                    ((residue_num+1)*residue_num)/2 + residue2_num - residue_num - 1
                            # print "%i | %i -> %i" % (residue_num, residue2_num, index)
                            all_states_interactions[md, ph_number, index] = dG_max - dG_min


                    for (residue_num, residue_kbp, resname_kbp, nr_of_states) in residue_list_ext_ref:

                        # if residue_kbp in ['DPP-19_A', 'DPP-21_A', 'DPP-40_A']:
                        #     # print residue_kbp + "skipped"
                        #     continue

                        protE_res_sum = 0.0

                        #####################
                        ### INTRINSIC PKA ###
                        #####################
                        # letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O']
                        # for state in range(len(titratable_residues[resname_kbp])):
                        #     state_letter = letters[state]
                        #     occ = self.kbp_results[md][frame_nr]['curves'][residue_kbp][state_letter][ph_number]
                        # for state in range(1, len(titratable_residues[resname_kbp])):
                        for state in range(1, nr_of_states):
                            # occ = kbp_results_frame[residue_kbp]['states'][state][ph_number]
                            occ = occs_np[residue_num, state, ph_number]

                            state_name = titratable_residues[resname_kbp][state]['name']
                            state_pka  = titratable_residues[resname_kbp][state]['pka']
                            # slope_sign ->  0 for R/0 state
                            #            ->  1 for P   state
                            #            -> -1 for D   state
                            # pka_int = pkaint_frame[residue_kbp][state] - pkaint_frame[residue_kbp][1]
                            pka_int = pkaint_frame[residue_kbp][state]
                            # print "%i -> %.2f" % (state, pka_int)
                            if state_name == 'P':
                                dG = (pka_int + 1.0 * RT_ln10_ph) * occ
                            elif state_name == 'D':
                                dG = (pka_int - 1.0 * RT_ln10_ph) * occ
                            elif state_name == '0':
                                dG = pka_int * occ
                            else:
                                assert("Sorry, state type %s is not yet implemented." % state_name)

                            protE_sum     += dG
                            protE_res_sum += dG

                            # Interaction energy.
                            # Iterate over all other residues.
                            # if state_type != 'R' and occ > 0.001:
                            if  occ > 0.001:
                                for (residue2_num, residue2_kbp, resname2_kbp, nr_of_states2) in residue_list_ext:
                                # for residue2_num, residue2_kbp in enumerate(residue_list):

                                    # residue2_num += residue_num+1  # enumerate starts with 0.
                                    if residue_num == residue2_num:
                                        continue

                                    # resname2_kbp = residue2_kbp[0:residue2_kbp.find('-')]

                                    # if resname2_kbp in ['LYS', 'ARG', 'HSP']:
                                    # if resname2_kbp in ['DPP-98_E']:
                                    #     continue

                                    # for state2 in range(1, len(titratable_residues[resname2_kbp])):
                                    for state2 in  range(1, nr_of_states2):
                                        # occ2 = kbp_results_frame[residue2_kbp]['states'][state2][ph_number]
                                        occ2 = occs_np[residue2_num, state2, ph_number]

                                        # state2_type = titratable_residues[resname2_kbp][state2]['name']
                                        # if state2_type != 'R':
                                        #     # R <-> x and R <-> 0 are zero by definition.
                                        ###################
                                        ### INTERACTION ###
                                        ###################
                                        if residue2_num > residue_num:
                                            # dG_interact = g_frame[residue_kbp, residue2_kbp, state, state2] * occ * occ2

                                            dG_interact = g_frame_np[residue_num, residue2_num, state, state2] * occ * occ2
                                            protE_sum += dG_interact
                                            # inet.mod_interaction(residue_num, residue2_num, dG_interact, ph_number, md)
                                            # mod_interaction(inet, residue_num, residue2_num, dG_interact, ph_number, md)
                                            # index = nr_of_residues*residue_num - \
                                            #         int((residue_num+1.0)*residue_num/2.0) + residue2_num - residue_num - 1
                                            # index = nr_of_residues*residue_num - \
                                            #         ((residue_num+1)*residue_num)/2 + residue2_num - residue_num - 1
                                            # # print "%i | %i -> %i" % (residue_num, residue2_num, index)
                                            # all_states_interactions[md, ph_number, index] += dG_interact
                                        else:
                                            # dG_interact = g_frame[residue2_kbp, residue_kbp, state2, state] * occ * occ2
                                            dG_interact = g_frame_np[residue2_num, residue_num, state2, state] * occ * occ2
                                        protE_res_sum += dG_interact

                        protE_res[frame_nr][ph][residue_kbp] = protE_res_sum
                    protE[frame_nr][ph] = protE_sum


                kbp_results = KbpResult()
                kbp_results.g = g_frame_np
                kbp_results.occs = occs_np
                pkaint_tmp = {}
                for residue_local in pkaint_frame.keys():
                    pkaint_tmp[residue_local] = np.array(pkaint_frame[residue_local], dtype=np.float)
                kbp_results.pkaint = pkaint_tmp
                frame_kbp_results.add_job(kbp_results, frame_nr)

                # print time.time() - start


            f = open(pickle_filename, 'wb')
            pickle.dump((protE, protE_res, frame_kbp_results), f, protocol=2)
            f.close()

        md_protE.append(protE)
        md_protE_res.append(protE_res)
        md_kbp_results.append(frame_kbp_results)

    return md_protE, md_protE_res, inet, md_kbp_results


@cython.boundscheck(False) # turn of bounds-checking for entire function
cpdef calc_protonation_energy(kbp_result):
    """ Sets kbp_result.protonation_energy.

    @param kbp_result: KbpResult object
    @return: None
    """
    assert(isinstance(kbp_result, KbpResult))

    # start = time.time()

    cdef float occ, occ2
    cdef float dG_interact, dG
    cdef float pka_int
    cdef float RT_ln10_ph
    cdef float protE_sum, protE_res_sum
    cdef int state, state2
    cdef int nr_of_states, nr_of_states2
    cdef int residue_num, residue2_num
    cdef int ph_number
    cdef np.ndarray[np.float_t, ndim=4] g_matrix_np

    cdef int nr_of_residues

    R = 8.3144621                             # J/(K*mol)
    # kb = 1.3806488e-23                      # J/K
    # RT = R * 300.0 / 1000.0                 # kJ/mol
    RT_ln10 = R * 300.0 / 1000.0 * np.log(10) # kJ/mol

    # pickle_filename = jobfolder + 'prot_energies_detailed%s.pickle' % self.kbp_suffix
    # if os.path.exists(pickle_filename) and not restart:
    #     print("Reading prot_energies_detailed cPickle file for: %s" % jobfolder)

        # f = open(pickle_filename, 'rb')
        # (protE, protE_res, frame_kbp_results) = pickle.load(f)
        # f.close

    # else:
    # print "Calculating protonation energies for: " + kbp_result.descr.folder

    residue_list = kbp_result.descr.sorted_residue_list
    nr_of_residues = len(residue_list)
    ph_values = kbp_result.descr.ph_values
    g_matrix = kbp_result.g
    pkaint = kbp_result.pkaint
    occs = kbp_result.occs
    titratable_residues = kbp_result.descr.titratable_residues
    residue_list_ext = kbp_result.descr.residue_list_ext
    residue_list_ext_ref = kbp_result.descr.residue_list_ext



    max_nr_of_states = 0
    for (residue_num, residue_kbp, resname_kbp, nr_of_states) in residue_list_ext_ref:
        if nr_of_states > max_nr_of_states:
            max_nr_of_states = nr_of_states
    g_matrix_np = np.zeros((nr_of_residues, nr_of_residues, max_nr_of_states-1, max_nr_of_states-1))
    for (residue_num, residue_kbp, resname_kbp, nr_of_states) in residue_list_ext:
        for state in range(0, nr_of_states-1):
            for (residue2_num, residue2_kbp, resname2_kbp, nr_of_states2) in residue_list_ext[residue_num+1:]:
                for state2 in range(0, nr_of_states2-1):
                    g_matrix_np[residue_num, residue2_num, state, state2] \
                        = g_matrix[residue_num][state][residue2_num][state2]

    protE     = np.zeros(len(ph_values), dtype=np.float)
    protE_res = {}
    for (residue_num, residue_kbp, resname_kbp, nr_of_states) in residue_list_ext_ref:
        protE_res[residue_kbp] = np.zeros(len(ph_values), dtype=np.float)
    for ph_number, ph in enumerate(ph_values):
        RT_ln10_ph = RT_ln10 * ph

        protE_sum = 0.0
        for (residue_num, residue_kbp, resname_kbp, nr_of_states) in residue_list_ext_ref:

            # if residue_kbp not in ['DPP-19_A', 'DPP-21_A', 'DPP-40_A']:
            # if residue_kbp not in ['EPP-57_A']:
                # print residue_kbp + "skipped"
                # continue

            protE_res_sum = 0.0

            #####################
            ### INTRINSIC PKA ###
            #####################
            for state in range(0, nr_of_states-1):
                occ = occs[residue_num][state+1, ph_number]

                state_name = titratable_residues[resname_kbp][state+1]['name']
                # slope_sign ->  0 for R/0 state
                #            ->  1 for P   state
                #            -> -1 for D   state
                pka_int = pkaint[residue_kbp][state+1]
                if state_name == 'P':
                    dG = (pka_int + 1.0 * RT_ln10_ph) * occ
                elif state_name == 'D':
                    dG = (pka_int - 1.0 * RT_ln10_ph) * occ
                elif state_name == '0':
                    dG = pka_int * occ
                else:
                    error = "Sorry, state type %s is not yet implemented." % state_name
                    raise AssertionError(error)

                protE_sum     += dG
                protE_res_sum += dG

                # Interaction energy.
                # Iterate over all other residues.
                # if state_type != 'R' and occ > 0.001:
                if  occ > 0.001:
                    for (residue2_num, residue2_kbp, resname2_kbp, nr_of_states2) in residue_list_ext:
                        if residue_num == residue2_num:
                            continue

                        for state2 in  range(0, nr_of_states2-1):
                            occ2 = occs[residue2_num][state2+1, ph_number]

                            ###################
                            ### INTERACTION ###
                            ###################
                            if residue2_num > residue_num:
                                dG_interact = g_matrix_np[residue_num, residue2_num, state, state2] * occ * occ2
                                # dG_interact = g_matrix[residue_num][state][residue2_num][state2] * occ * occ2
                                protE_sum += dG_interact
                            else:
                                dG_interact = g_matrix_np[residue2_num, residue_num, state2, state] * occ * occ2
                                # dG_interact = g_matrix[residue2_num][state2][residue_num][state] * occ * occ2
                            protE_res_sum += dG_interact

            protE_res[residue_kbp][ph_number] = protE_res_sum
        protE[ph_number] = protE_sum

        # if ph_number > 2:
        #     break
    # print time.time() - start

        # f = open(pickle_filename, 'wb')
        # pickle.dump((protE, protE_res, frame_kbp_results), f, protocol=2)
        # f.close()

    kbp_result.prot_energy = protE
    kbp_result.residue_prot_energy = protE_res

    return



