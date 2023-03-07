# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 12:03:07 2012

@author: Tim Meyer
"""

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

from kbp2.kbp_results import *

def analyse_tension(kbp_result, sasa_factor, min_delta_e=-20.0):
    return analyse_tension_local(kbp_result, sasa_factor, min_delta_e=min_delta_e)

cpdef analyse_tension_local(kbp_result, sasa_factor, min_delta_e=-20.0):
    assert(isinstance(kbp_result, FrameKbpResults))

    # Average number of atoms in a sphere of 10 Angström (0.1 Atoms/A^3)
    # conf_energy_multiplier = 4.0/3.0*3.1415 * 10.0**3 * 0.1
    conf_energy_multiplier = 1.0

    ph_values = kbp_result.descr.ph_values
    residue_list_ext = kbp_result.descr.residue_list_ext
    md_weights = kbp_result.weights
    titratable_residues = kbp_result.descr.titratable_residues

    tension = np.zeros((len(residue_list_ext), len(ph_values)))
    best_state = np.zeros((len(residue_list_ext), len(ph_values)), dtype=int)
    best_md = np.zeros((len(residue_list_ext), len(ph_values)), dtype=int)

    for ph_number, ph in enumerate(ph_values):
        for (residue_num, residue_kbp, resname_kbp, nr_of_states) in residue_list_ext:
            # Get current energy of the residue and its protonation state.

            occs = kbp_result.combined_results.occs[residue_num][:,ph_number]
            current_state = np.argmax(occs)
            current_state_name = titratable_residues[resname_kbp][current_state]['name']
            current_state_name.replace('R', '0')
            current_md = np.argmax(md_weights[ph_number])

            current_res_prot_energy = 0.0
            # current_res_sasa_corr = 0.0
            # current_res_confE = 0.0
            for md_number, md_kbp_result in kbp_result:
                current_res_prot_energy += md_weights[ph_number][md_number] * \
                                           kbp_result.kbp_results[md_number].residue_prot_energy[residue_kbp][ph_number]
                # current_res_sasa_corr += md_weights[ph_number][md_number] * sasa_factor * \
                #                          kbp_result.kbp_results[md_number].residue_sasa[residue_num]
                # current_res_confE += md_weights[ph_number][md_number] * conf_energy_multiplier * \
                #                      np.sum(kbp_result.kbp_results[md_number].conf_energies['residue_energy'][residue_num])


            # Compare with the energy of the individual MDs
            # If protonation state differs, suggest this protonation state for the next MD
            best_alternative_state = None
            best_alternative_state_name = None
            best_alternative_delta_e = None
            best_alternative_md = None
            for md_number, md_kbp_result in kbp_result:
                occs = md_kbp_result.occs[residue_num][:,ph_number]
                md_state = np.argmax(occs)
                md_state_name = titratable_residues[resname_kbp][md_state]['name']
                md_state_name.replace('R', '0')

                res_prot_energy = md_kbp_result.residue_prot_energy[residue_kbp][ph_number]
                # res_sasa_corr = sasa_factor * md_kbp_result.residue_sasa[residue_num]
                # res_confE = conf_energy_multiplier * np.sum(md_kbp_result.conf_energies['residue_energy'][residue_num])

                if current_state != md_state:
                    if (current_state_name != md_state_name) or (residue_kbp in ['HSP']):
                        delta_e = res_prot_energy - current_res_prot_energy
                        # delta_e += res_sasa_corr - current_res_sasa_corr
                        # delta_e += res_confE - current_res_confE
                        if delta_e < min_delta_e:
                            if (best_alternative_delta_e is None) or (delta_e < best_alternative_delta_e):
                                best_alternative_state = md_state
                                best_alternative_state_name = md_state_name
                                best_alternative_delta_e = delta_e
                                best_alternative_md = md_number

            if best_alternative_state is not None:
                # print "%.2f - %s  %i (%s) -> %i (%s): %.2f" % (ph, residue_kbp, current_state, current_state_name, best_alternative_state, best_alternative_state_name, best_alternative_delta_e)
                tension[residue_num, ph_number] = best_alternative_delta_e
                best_state[residue_num, ph_number] = best_alternative_state
                best_md[residue_num, ph_number] = best_alternative_md
            else:
                # best_state[residue_num, ph_number] = current_state
                best_state[residue_num, ph_number] = -1
                best_md[residue_num, ph_number] = -1

    return tension, best_state, best_md



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
            # if residue_kbp not in ['NTE-1_A']:
                # print residue_kbp + "skipped"
                # continue
            # if resname_kbp in ['NTE']:
            #     continue

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

                # if resname_kbp in ['NTE']:
                #     dG *= 5.0

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

                            # scaling = 1.0
                            # sb_acid = ['EPP', 'DPP']
                            # # sb_base = ['ARG', 'LYS']
                            # sb_base = ['ARG']
                            # if ((resname_kbp in sb_acid) and (resname2_kbp in sb_base)) or \
                            #         ((resname_kbp in sb_base) and (resname2_kbp in sb_acid)):
                            #      scaling = 0.8



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



