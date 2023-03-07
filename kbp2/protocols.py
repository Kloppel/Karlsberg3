# coding=utf-8
import os
import kbp2


def h_min(charmm_pac, kbp2_settings, iteration):
    """
    Hydrogen minimization for protocol with ph=7
    """

    charmm_pac.backup_settings()

    info_to_log = ("Only hydrogens are minimzied! Iteration %i" % iteration)
    kbp2_settings.log(info_to_log)

    iteration_folder = charmm_pac.workdir
    workdir = iteration_folder+ 'modelling/'
    if not os.path.exists(workdir):
        os.mkdir(workdir)
    charmm_pac.workdir = workdir

    if not charmm_pac.charmm_instructions['minimize']:
        charmm_pac.force_energy_min()
    else:
        error = "Minimization instructions are not empty!"
        raise AssertionError(error)

    charmm_pac.run_charmm()

    charmm_pac.update_coords_with_modelled_structure()
    charmm_pac.restore_settings()

    #debug
    if not charmm_pac.charmm_instructions['minimize']:
        charmm_pac.charmm_instructions['minimize'] = []
    else:
        error = "Minimization instructions are not empty!"
        raise AssertionError(error)

    # # debug
    # minimization_checks = ['minimize', 'minimize_selections', 'gaps']
    # for min_check in minimization_checks:
    #     if charmm_pac.charmm_instructions[min_check]:
    #         charmm_pac.charmm_instructions[min_check] = []
    #     else:
    #         print min_check, charmm_pac.charmm_instructions[min_check]
    #         error = "Minimization instructions are not empty!"
    #         raise AssertionError(error)
    # if charmm_pac.charmm_config['missing_atom_ids']:
    #     charmm_pac.charmm_config['missing_atom_ids'] = []
    # else:
    #     print 'missing_atom_ids', charmm_pac.charmm_config['missing_atom_ids']
    #     error = "Minimization instructions are not empty!"
    #     raise AssertionError(error)


def open_sb_acids(charmm_pac, kbp2_settings, iteration):
    """
    Salt-bridge opening for protocol with ph=-10
    """

    info_to_log = ("Salt bridge opening acids! Iteration %i" % iteration)
    kbp2_settings.log(info_to_log)

    iteration_folder = charmm_pac.workdir
    workdir_mod = iteration_folder+ 'open_sb/'
    workdir_min = iteration_folder+ 'minimize_sb/'

    cutoff = kbp2_settings.sb_cutoff

    if iteration == 1:
        residues_in_salt_bridges = kbp2.modelling.open_salt_bridges(workdir_mod, workdir_min, charmm_pac,
                                                                        cutoff = cutoff, opening_acids = True, tyrosines = True)
        kbp2_settings.salt_bridges = list(residues_in_salt_bridges)
        kbp2_settings.log('Salt bridges Are below')
        kbp2_settings.log(residues_in_salt_bridges)
    else:
        workdir_mod = iteration_folder+ 'modelling/'
        kbp2_settings.log("Regular modelling -> Only hydrogens will be minimized in this iteration!")
        residues_in_salt_bridges = kbp2.modelling.open_salt_bridges(workdir_mod, workdir_min, charmm_pac)
        if residues_in_salt_bridges:
            error = "This should be empty list of salt bridges now cause they are not supposed to be found in iteration %i" % iteration
            kbp2_settings.log(error)
            raise AssertionError (error)


def open_sb_bases(charmm_pac, kbp2_settings, iteration):
    """
    Salt-bridge opening for protocol with ph=20
    """

    info_to_log = ("Salt bridge opening bases! Iteration %i" % iteration)
    kbp2_settings.log(info_to_log)

    iteration_folder = charmm_pac.workdir
    workdir_mod = iteration_folder+ 'open_sb/'
    workdir_min = iteration_folder+ 'minimize_sb/'

    cutoff = kbp2_settings.sb_cutoff

    if iteration == 1:
        residues_in_salt_bridges = kbp2.modelling.open_salt_bridges(workdir_mod, workdir_min, charmm_pac, cutoff = cutoff, opening_bases = True, tyrosines = True)
        kbp2_settings.salt_bridges = list(residues_in_salt_bridges)
        kbp2_settings.log('Salt bridges Are below')
        kbp2_settings.log(residues_in_salt_bridges)
    else:
        workdir_mod = iteration_folder+ 'modelling/'
        kbp2_settings.log("Regular modelling -> Only hydrogens will be minimized in this iteration!")
        residues_in_salt_bridges = kbp2.modelling.open_salt_bridges(workdir_mod, workdir_min, charmm_pac)
        if residues_in_salt_bridges:
            error = "This should be empty list of salt bridges now cause they are not supposed to be found in iteration %i" % iteration
            kbp2_settings.log(error)
            raise AssertionError (error)




#######################################################################################################################################################################
def open_sb_acids_md(charmm_pac, kbp2_settings, iteration):

    info_to_log = ("Salt bridge opening acids! Iteration %i" % iteration)
    # print info_to_log
    kbp2_settings.log(info_to_log)

    iteration_folder = charmm_pac.workdir
    workdir_mod = iteration_folder+ 'open_sb/'

    cutoff = kbp2_settings.sb_cutoff

    if iteration == 1:
        residues_in_salt_bridges = kbp2.modelling.open_salt_bridges_md(workdir_mod, charmm_pac, cutoff = cutoff, opening_acids = True, gbsw = True, tyrosines = True)
        kbp2_settings.salt_bridges = list(residues_in_salt_bridges)
        kbp2_settings.log('Salt bridges Are below')
        kbp2_settings.log(residues_in_salt_bridges)
    else:
        workdir_mod = iteration_folder+ 'modelling/'
        kbp2_settings.log("Regular modelling -> Only hydrogens will be minimized in this iteration!")
        residues_in_salt_bridges = kbp2.modelling.open_salt_bridges_md(workdir_mod, charmm_pac)
        if residues_in_salt_bridges:
            error = "This should be empty list of salt bridges now cause they are not supposed to be found in iteration %i" % iteration
            kbp2_settings.log(error)
            raise AssertionError (error)


def open_sb_bases_md(charmm_pac, kbp2_settings, iteration):

    info_to_log = ("Salt bridge opening bases! Iteration %i" % iteration)
    # print info_to_log
    kbp2_settings.log(info_to_log)

    iteration_folder = charmm_pac.workdir
    workdir_mod = iteration_folder+ 'open_sb/'

    cutoff = kbp2_settings.sb_cutoff

    if iteration == 1:
        residues_in_salt_bridges = kbp2.modelling.open_salt_bridges_md(workdir_mod, charmm_pac, cutoff = cutoff, opening_bases = True, gbsw = True, tyrosines = True)
        kbp2_settings.salt_bridges = list(residues_in_salt_bridges)
        kbp2_settings.log('Salt bridges Are below')
        kbp2_settings.log(residues_in_salt_bridges)
    else:
        workdir_mod = iteration_folder+ 'modelling/'
        kbp2_settings.log("Regular modelling -> Only hydrogens will be minimized in this iteration!")
        residues_in_salt_bridges = kbp2.modelling.open_salt_bridges_md(workdir_mod, charmm_pac)
        if residues_in_salt_bridges:
            error = "This should be an empty list of salt bridges now because they are not supposed to be found in iteration %i" % iteration
            kbp2_settings.log(error)
            raise AssertionError (error)


# def special_pac_ph51(charmm_pac, kbp2_settings, iteration):
#
#     '''Emanuele'''
#
#     info_to_log = ("Special pac! Iteration %i" % iteration)
#     kbp2_settings.log(info_to_log)
#
#     iteration_folder = charmm_pac.workdir
#     workdir = iteration_folder+ 'modelling/'
#
#     kbp2.modelling.pac_ph51(workdir, charmm_pac)
#
#
# def special_pac_ph52(charmm_pac, kbp2_settings, iteration):
#
#     '''Emanuele'''
#
#     info_to_log = ("Special pac! Iteration %i" % iteration)
#     kbp2_settings.log(info_to_log)
#
#     iteration_folder = charmm_pac.workdir
#     workdir = iteration_folder+ 'modelling_spec/'
#
#     kbp2.modelling.pac_ph52(workdir, charmm_pac)
#
# def special_pac_ph41(charmm_pac, kbp2_settings, iteration):
#
#     '''David'''
#
#     info_to_log = ("Special pac! Iteration %i" % iteration)
#     kbp2_settings.log(info_to_log)
#
#     iteration_folder = charmm_pac.workdir
#     workdir = iteration_folder+ 'modelling_spec/'
#
#     kbp2.modelling.pac_ph41(workdir, charmm_pac)
#
# def special_pac_ph42(charmm_pac, kbp2_settings, iteration):
#
#     '''David'''
#
#     info_to_log = ("Special pac! Iteration %i" % iteration)
#     kbp2_settings.log(info_to_log)
#
#     iteration_folder = charmm_pac.workdir
#     workdir = iteration_folder+ 'modelling_spec/'
#
#     kbp2.modelling.pac_ph42(workdir, charmm_pac)
#
# def special_pac_ph43(charmm_pac, kbp2_settings, iteration):
#
#     '''David'''
#
#     info_to_log = ("Special pac! Iteration %i" % iteration)
#     kbp2_settings.log(info_to_log)
#
#     iteration_folder = charmm_pac.workdir
#     workdir = iteration_folder+ 'modelling_spec/'
#
#     kbp2.modelling.pac_ph43(workdir, charmm_pac)
#
# def special_pac_ph44(charmm_pac, kbp2_settings, iteration):
#
#     '''David'''
#
#     info_to_log = ("Special pac! Iteration %i" % iteration)
#     kbp2_settings.log(info_to_log)
#
#     iteration_folder = charmm_pac.workdir
#     workdir = iteration_folder+ 'modelling_spec/'
#
#     kbp2.modelling.pac_ph44(workdir, charmm_pac)



