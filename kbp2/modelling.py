# coding=utf-8

import os
import re
import kbp2

def open_salt_bridges(workdir_mod, workdir_min, charmm_pac, cutoff = 4, opening_acids = False, opening_bases = False,
                      gbsw = False, tyrosines = True):
    """
    Modelling of the salt-bridge opening at ph=-10 or ph=20. Special patches are applied to acids or bases
    (see parameters below) in order to open a salt-bridge.

    :param cutoff: previously defined cutoff to determine acid/base pairs
    :param opening_acids: acids will be affected by the salt-bridge opening procedure
    :param opening_bases: bases will be affected by the salt-bridge opening procedure
    :param gbsw: usage of implicit solvent during modelling
    :param tyrosines: will tyrosines be included and regarded as salt-bridge/hydrogen bond partners
    """

    if opening_acids == True and opening_bases == True:
        raise AssertionError("Cannot choose two methods at once!")

    if not os.path.exists(workdir_mod):
        os.mkdir(workdir_mod)
    charmm_pac.workdir = workdir_mod

    dir = os.path.dirname(__file__)
    add_top = dir + '/additional_files/sb_pat.rtf'
    charmm_pac.top.append(add_top)
    charmm_pac.parse_top()

    charmm_pac_ssp = charmm_pac.structure.copy()

    if not charmm_pac_ssp.par_read:
        for par in charmm_pac.par:
            charmm_pac_ssp.read_par(par)
    if not charmm_pac_ssp.top_content:
        for top in charmm_pac.top:
            charmm_pac_ssp.read_top(top)

    list_of_residues_in_sb = []

    if opening_acids == False and opening_bases == False:
        # print "Regular modelling -> Only hydrogens will be minimized!"
        pass
    else:
        salt_struct = kbp2.file_parser.MDAnalysis_ssp(charmm_pac_ssp)
        if tyrosines:
            (residues_in_salt_bridges1, residues_not_in_salt_bridges1, salt_bridges1) = \
                salt_struct.get_residues_in_sb(tyr_acid = True, cutoff = cutoff)
            (residues_in_salt_bridges2, residues_not_in_salt_bridges2, salt_bridges2) = \
                salt_struct.get_residues_in_sb(tyr_base = True, cutoff = cutoff)
            for residue in residues_in_salt_bridges1:
                if residue in residues_in_salt_bridges2:
                    continue
                else:
                    residues_in_salt_bridges2.append(residue)
            # list_of_residues_in_sb = list(residues_in_salt_bridges2)
            list_of_residues_in_sb = residues_in_salt_bridges2
        else:
            (residues_in_salt_bridges, residues_not_in_salt_bridges, salt_bridges) = \
                salt_struct.get_residues_in_sb(cutoff = cutoff)
            # list_of_residues_in_sb = list(residues_in_salt_bridges)
            list_of_residues_in_sb = residues_in_salt_bridges

        if len(list_of_residues_in_sb) == 0:
            raise AssertionError("No salt bridges!")

    prev_titr_residue_dict = charmm_pac.get_titr_residue_dict()
    new_titr_residue_dict  = charmm_pac.get_titr_residue_dict()

    prev_titr_residues = charmm_pac.get_titr_residues()
    new_titr_residues = charmm_pac.get_titr_residues()

    template_state = {'charge' : 0,

                      'patch' : None,
                      'external_patches' : None,
                      'rename' : None,
                      'special' : None}


    resname = 'EPP'
    state = dict(template_state)
    state['charge'] = 1
    state['patch']  = 'EPPS'
    new_titr_residues[resname].append(state)

    resname = 'DPP'
    state = dict(template_state)
    state['charge'] = 1
    state['patch']  = 'DPPS'
    new_titr_residues[resname].append(state)

    resname = 'LYS'
    state = dict(template_state)
    state['charge'] = -1
    state['patch']  = 'LYSS'
    new_titr_residues[resname].append(state)

    resname = 'ARG'
    state = dict(template_state)
    state['charge'] = -1
    state['patch']  = 'ARGS'
    new_titr_residues[resname].append(state)

    resname = 'HSP'
    state = dict(template_state)
    state['charge'] = -1
    state['patch']  = 'HSPS'
    new_titr_residues[resname].append(state)

    resname = 'TYR'
    state = dict(template_state)
    state['charge'] = 1
    state['patch']  = 'TYRS'
    new_titr_residues[resname].append(state)

    #CTER NTER
    resname = 'CTE'
    state = dict(template_state)
    state['charge'] = 1
    state['patch']  = 'CTES'
    new_titr_residues[resname].append(state)

    resname = 'NTE'
    state = dict(template_state)
    state['charge'] = -1
    state['patch']  = 'NTES'
    new_titr_residues[resname].append(state)


    charmm_pac.set_titr_residues(new_titr_residues)
    charmm_pac.apply_titr_residue_dict(new_titr_residue_dict)

    charmm_pac.charmm_instructions['minimize'] = []

    if len(list_of_residues_in_sb) > 0:
        for residue in list_of_residues_in_sb:
            (resname, resid, segname) = re.split(r'[-_]', residue)
            resid = int(resid)

            if opening_acids:
                if resname == 'EPP':
                    charmm_pac.set_prot_residue(('EPP', resid, segname), patch='EPPS')
                if resname == 'DPP':
                    charmm_pac.set_prot_residue(('DPP', resid, segname), patch='DPPS')
                if resname == 'CTE':
                    charmm_pac.set_prot_residue(('CTE', resid, segname), patch='CTES')
                if resname == 'TYR':
                    charmm_pac.set_prot_residue(('TYR', resid, segname), patch='TYRS')

            if opening_bases:
                if resname == 'LYS':
                    charmm_pac.set_prot_residue(('LYS', resid, segname), patch='LYSS')
                if resname == 'HSP':
                    charmm_pac.set_prot_residue(('HSP', resid, segname), patch='HSPS')
                if resname == 'NTE':
                    charmm_pac.set_prot_residue(('NTE', resid, segname), patch='NTES')
                if resname == 'ARG':
                    charmm_pac.set_prot_residue(('ARG', resid, segname), patch='ARGS')

            residue_to_minimize = charmm_pac.structure.struct[segname][resid]
            charmm_pac.charmm_instructions['minimize'].append(residue_to_minimize)

    if gbsw:

        ### Definition of GB ###
        gbsw_commands = []
        gbsw_commands.append("stream \"/scratch/scratch/tmeyer/CHARMM_NAMD/toppar_36/gbsw/radius_gbsw.str\"")
        gbsw_commands.append("scalar wmain statistics select .not. type H* end")
        gbsw_commands.append("deXine check select (.not type H* ) .and. ( prop wmain .eq. 0.0 ) show end")
        gbsw_commands.append("if ?nsel ne 0  stop       !some heavy atom have a zero radius")
        gbsw_commands.append("GBSW conc 0.1 sgamma 0.03 GBenergy epsp 4 MOLSURF")

        ### adding GB ###
        command_block = " "
        for command in gbsw_commands:
            command_block += '\n'
            command_block += command

        charmm_pac.add_charmm_command(command_block, adj_task='hbuild')

    # print "Modelling salt bridge opening NOW!"
    charmm_pac.charmm_instructions['backbone_fixed'] = True
    charmm_pac.force_energy_min()
    charmm_pac.run_charmm()

    charmm_pac.apply_titr_residue_dict(prev_titr_residue_dict)
    charmm_pac.set_titr_residues(prev_titr_residues)

    charmm_pac.update_coords_with_modelled_structure()

    if len(list_of_residues_in_sb) > 0:
        # print "Minimization after salt bridge opening!"

        if not os.path.exists(workdir_min):
            os.mkdir(workdir_min)
        charmm_pac.workdir = workdir_min

        charmm_pac.charmm_instructions['backbone_fixed'] = True
        charmm_pac.run_charmm()
        charmm_pac.update_coords_with_modelled_structure()

    charmm_pac.charmm_instructions['minimize'] = []
    charmm_pac.charmm_instructions['backbone_fixed'] = False

    return list_of_residues_in_sb


def get_charmm_command_md(open_acids=False, open_bases=False):

    if open_acids:
        charmm_command_md ="""


!titratable (protonated) residue atoms

!tests comment

define EPPO select (resname EPP .and. type OE*) end
define EPPH select (resname EPP .and. type HE*) end
define DPPO select (resname DPP .and. type OD*) end
define DPPH select (resname DPP .and. type HD*) end
define CTERO select (type OT*) end

define TYRO select (type OH) end
define TYRH select (type HH) end

define LYSN select (type NZ) end
define LYSH select (resname LYS .and. type HZ3) end
define HISN select (resname HSP .and. (type ND1 .or. type NE2)) end
define HISH select (resname HSP .and. (type HD1 .or. type HE2)) end
define NTERN select (resid 1 .and. type N) end
define NTERH select (type HT*) end

!nontitratable residue atoms that are involved in SB
define ARGN select (type NH1 .or. type NH2) end
define ARGH select (resname ARG .and. type HH*) end
define SERO select (type OG) end
define CYSS select (type SG .and. .not. .bonded. type SG) end
define THRO select (type OG1) end

!groups to adjust to open SB
define ACIDres select (EPPO .or. DPPO .or. CTERO) end
define ACIDhyd select (EPPH .or. DPPH) end
define BASEres select (ARGN .or. LYSN .or. NTERN) end
define BASEhyd select (ARGH .or. LYSH .or. NTERH) end

!scalar charge set -0.20 sele ((DPPO .or. EPPO) .and. sb_resi) end      ! symm. charges for SB ACID oxygens
!scalar charge set 0.11 sele ((DPPH .or. EPPH).and. sb_resi) end        ! symm. charges for SB ACID hydrogens

!scalar charge set -0.80 sele (ARGN .and. sb_resi) end		           ! set charge for BASE atoms
!scalar charge set 0.46 sele (ARGH .and. sb_resi) end    		       ! set charge for BASE hydrogens

!scalar charge set -0.10 sele (TYRO) end                  ! set charge for TYR oxygen
!scalar charge set 0.43 sele (TYRH) end                   ! set charge for TYR hydrogen

shake bonh param

cons fix select .not. (sb_resi .or. resname TYR) .or. backbone end

dynamics start time 0.002 nstep 5000 -
   firstt 100.0 finalt 100.0 teminc 10.0 ihtfrq 100 -                  !temp control
   ieqfrq 2000 ichecw 1 twindl -5.0 twindh +5.0 iasors 0 -             !temp control
   atom cdie fshift vshift cutnb 14.0 ctofnb 12.0 -                    !nonbond
   inbfrq -1 ihbfrq -1

!COOR AXIS SELE ATOM A 20 CA END SELE ATOM A 20 CB END
!COOR ROTATE AXIS PHI 310.0 SELE RESID 20 .AND. .NOT. (backbone .and. type H) END
"""
    if open_bases:
        charmm_command_md = """
!titratable (deprotonated) residue atoms
define EPPO select (resname EPP .and. type OE*) end
define EPPH select (resname EPP .and. type HE*) end
define DPPO select (resname DPP .and. type OD*) end
define DPPH select (resname DPP .and. type HD*) end
define CTERO select (type OT*) end

define TYRO select (type OH) end
define TYRH select (type HH) end

define LYSN select (type NZ) end
define LYSH select (resname LYS .and. type HZ3) end
define HISN select (resname HSP .and. (type ND1 .or. type NE2)) end
define HISH select (resname HSP .and. (type HD1 .or. type HE2)) end
define NTERN select (resid 1 .and. type N) end
define NTERH select (type HT*) end

!nontitratable residue atoms that are involved in SB
define ARGN select (type NH1 .or. type NH2) end
define ARGH select (resname ARG .and. type HH*) end
define SERO select (type OG) end
define CYSS select (type SG .and. .not. .bonded. type SG) end
define THRO select (type OG1) end

!groups to adjust to open SB
define ACIDres select (EPPO .or. DPPO .or. CTERO) end
define ACIDhyd select (EPPH .or. DPPH) end
define BASEres select (ARGN .or. LYSN .or. NTERN) end
define BASEhyd select (ARGH .or. LYSH .or. NTERH) end

!scalar charge set -0.76 sele ((DPPO .or. EPPO) .and. sb_resi) end      ! symm. charges for SB ACID oxygens
!scalar charge set 0.00 sele ((DPPH .or. EPPH).and. sb_resi) end       ! symm. charges for SB ACID hydrogens

!scalar charge set -1.00 sele (ARGN .and. sb_resi) end		           ! set charge for BASE atoms
!scalar charge set 0.27 sele (ARGH .and. sb_resi) end		           ! set charge for BASE hydrogens

!scalar charge set -1.00 sele (TYRO) end                  ! set charge for TYR oxygen
!scalar charge set -0.10 sele (TYRH) end                  ! set charge for TYR hydrogen

shake bonh param

cons fix select .not. (sb_resi .or. resname TYR) .or. backbone end

dynamics start time 0.002 nstep 5000 -
   firstt 100.0 finalt 100.0 teminc 10.0 ihtfrq 100 -                  !temp control
   ieqfrq 2000 ichecw 1 twindl -5.0 twindh +5.0 iasors 0 -             !temp control
   atom cdie fshift vshift cutnb 14.0 ctofnb 12.0 -                    !nonbond
   inbfrq -1 ihbfrq -1

!COOR AXIS SELE ATOM A 20 CA END SELE ATOM A 20 CB END
!COOR ROTATE AXIS PHI 310.0 SELE RESID 20 .AND. .NOT. (backbone .and. type H) END
"""
    return charmm_command_md


def open_salt_bridges_md(workdir_mod, charmm_pac, cutoff = 4, opening_acids = False, opening_bases = False, gbsw = False, tyrosines = True):

    if opening_acids == True and opening_bases == True:
        raise AssertionError("Cannot choose two methods at once!")

    if not os.path.exists(workdir_mod):
        os.mkdir(workdir_mod)
    charmm_pac.workdir = workdir_mod

    #add_top = '/user/jdragelj/python/karlsberg/kbp2/additional_files/sb_pat.rtf'
    #charmm_pac.top.append(add_top)
    #charmm_pac.parse_top()

    charmm_pac_ssp = charmm_pac.structure.copy()

    if not charmm_pac_ssp.par_read:
        for par in charmm_pac.par:
            charmm_pac_ssp.read_par(par)
    if not charmm_pac_ssp.top_content:
        for top in charmm_pac.top:
            charmm_pac_ssp.read_top(top)

    list_of_residues_in_sb = []

    if opening_acids == False and opening_bases == False:
        # print "Regular modelling -> Only hydrogens will be minimized!"
        pass
    else:
        salt_struct = kbp2.file_parser.MDAnalysis_ssp(charmm_pac_ssp)
        if tyrosines:
            (residues_in_salt_bridges1, residues_not_in_salt_bridges1, salt_bridges1) = salt_struct.get_residues_in_sb(tyr_acid = True, cutoff = cutoff)
            (residues_in_salt_bridges2, residues_not_in_salt_bridges2, salt_bridges2) = salt_struct.get_residues_in_sb(tyr_base = True, cutoff = cutoff)
            for residue in residues_in_salt_bridges1:
                if residue in residues_in_salt_bridges2:
                    continue
                else:
                    residues_in_salt_bridges2.append(residue)
            # list_of_residues_in_sb = list(residues_in_salt_bridges2)
            list_of_residues_in_sb = residues_in_salt_bridges2
        else:
            (residues_in_salt_bridges, residues_not_in_salt_bridges, salt_bridges) = salt_struct.get_residues_in_sb(cutoff = cutoff)
            # list_of_residues_in_sb = list(residues_in_salt_bridges)
            list_of_residues_in_sb = residues_in_salt_bridges

        # print "Residues in salt bridges are: ", list_of_residues_in_sb

        if len(list_of_residues_in_sb) == 0:
            raise AssertionError("No salt bridges!")

    prev_titr_residue_dict = charmm_pac.get_titr_residue_dict()
    prev_titr_residues = charmm_pac.get_titr_residues()

    charmm_pac.backup_settings()

    if len(list_of_residues_in_sb) > 0:

        selection_all = 'define sb_resi select ('
        for i, residue_in_sb in enumerate(list_of_residues_in_sb):
            resname, resid, segname = re.split('[-_]', residue_in_sb)
            selection_res = " (segid %s .and. resid %s) " % (segname, resid)
            selection_all += selection_res
            if i+1 != len(list_of_residues_in_sb):
                selection_all += ' .or. -'
            else:
                selection_all += ') end'
            selection_all += '\n'

        if opening_acids:
            charmm_command_md = get_charmm_command_md(open_acids=True)
        if opening_bases:
            charmm_command_md = get_charmm_command_md(open_bases=True)

        charmm_pac.add_charmm_command(charmm_command_md, adj_task='minimize_modeled')
        charmm_pac.add_charmm_command(selection_all, adj_task='minimize_modeled')

    if gbsw:

        print "GBSW MD initialized in %s!" % charmm_pac.workdir

        ### Definition of GB ###
        gbsw_commands = []
        gbsw_commands.append("stream \"/scratch/scratch/tmeyer/CHARMM_NAMD/toppar_36/gbsw/radius_gbsw.str\"")
        gbsw_commands.append("scalar wmain statistics select .not. type H* end")
        gbsw_commands.append("define check select ( .not type H* ) .and. ( prop wmain .eq. 0.0 ) show end")
        gbsw_commands.append("if ?nsel ne 0  stop       !some heavy atom have a zero radius")
        gbsw_commands.append("GBSW conc 0.1 sgamma 0.03 GBenergy epsp 4 MOLSURF")

        ### adding GB ###
        command_block = " "
        for command in gbsw_commands:
            command_block += '\n'
            command_block += command

        charmm_pac.add_charmm_command(command_block, adj_task='hbuild')

    # print "Modelling salt bridge opening NOW!"
    charmm_pac.charmm_instructions['backbone_fixed'] = True
    charmm_pac.force_energy_min()
    charmm_pac.run_charmm()

    charmm_pac.restore_settings()

    #remove charmm commands from charmm tasks
    # for i, task in enumerate(charmm_pac.charmm_config['tasks']):
    #     if task == 'minimize_modeled':
    #         if len(charmm_pac.charmm_config['tasks']) > i+2 :
    #             del charmm_pac.charmm_config['tasks'][i+2]
    #             del charmm_pac.charmm_config['tasks'][i+1]


    charmm_pac.apply_titr_residue_dict(prev_titr_residue_dict)
    charmm_pac.set_titr_residues(prev_titr_residues)

    charmm_pac.update_coords_with_modelled_structure()
    charmm_pac.charmm_instructions['minimize'] = []
    charmm_pac.charmm_instructions['backbone_fixed'] = False


    return list_of_residues_in_sb


# def pac_ph51(workdir, charmm_pac):
#
#     '''Emanuele'''
#
#     charmm_pac.backup_settings()
#     if not os.path.exists(workdir):
#         os.mkdir(workdir)
#     charmm_pac.workdir = workdir
#
#     charmm_pac.backup_settings()
#     charmm_pac.charmm_instructions['minimize'] = []
#     charmm_pac.charmm_instructions['do_minimize'] = False
#
#     ##########################
#     ### CHANGE PROTONATION ###
#     ##########################
#
#     #placeholder
#     residue_tuple = ('DPP', 19, 'A')
#     state = 0
#     charmm_pac.set_prot_residue(residue_tuple, state=state)
#
#
#     spec_command = get_sel_commands()
#
#     charmm_pac.add_charmm_command(spec_command, adj_task='hbuild')
#     charmm_pac.run_charmm()
#     charmm_pac.update_coords_with_modelled_structure()
#     charmm_pac.restore_settings()
#
#
#     #after everything minimization should be empty -> for the next iteration
#     if not charmm_pac.charmm_instructions['minimize']:
#         charmm_pac.charmm_instructions['minimize'] = []
#     else:
#         error = "Minimization instructions are not empty!"
#         raise AssertionError(error)
#
#     return
#
# def pac_ph52(workdir, charmm_pac):
#
#     '''Emanuele'''
#
#     charmm_pac.backup_settings()
#     if not os.path.exists(workdir):
#         os.mkdir(workdir)
#     charmm_pac.workdir = workdir
#
#     charmm_pac.backup_settings()
#     charmm_pac.charmm_instructions['minimize'] = []
#     charmm_pac.charmm_instructions['do_minimize'] = False
#
#     ##########################
#     ### CHANGE PROTONATION ###
#     ##########################
#
#     #placeholder
#     residue_tuple = ('DPP', 21, 'A')
#     state = 1
#     charmm_pac.set_prot_residue(residue_tuple, state=state)
#
#     spec_command = get_sel_commands()
#
#     charmm_pac.add_charmm_command(spec_command, adj_task='hbuild')
#     charmm_pac.run_charmm()
#     charmm_pac.update_coords_with_modelled_structure()
#     charmm_pac.restore_settings()
#
#
#     #after everything minimization should be empty -> for the next iteration
#     if not charmm_pac.charmm_instructions['minimize']:
#         charmm_pac.charmm_instructions['minimize'] = []
#     else:
#         error = "Minimization instructions are not empty!"
#         raise AssertionError(error)
#
#     return
#
# def pac_ph41(workdir, charmm_pac):
#
#     '''David'''
#
#     charmm_pac.backup_settings()
#     if not os.path.exists(workdir):
#         os.mkdir(workdir)
#     charmm_pac.workdir = workdir
#
#     charmm_pac.backup_settings()
#     charmm_pac.charmm_instructions['minimize'] = []
#     charmm_pac.charmm_instructions['do_minimize'] = False
#
#     ##########################
#     ### CHANGE PROTONATION ###
#     ##########################
#
#     #placeholder
#     residue_tuple = ('DPP', 19, 'A')
#     state = 0
#     charmm_pac.set_prot_residue(residue_tuple, state=state)
#
#     spec_command = get_sel_commands()
#
#     charmm_pac.add_charmm_command(spec_command, adj_task='hbuild')
#     charmm_pac.run_charmm()
#     charmm_pac.update_coords_with_modelled_structure()
#     charmm_pac.restore_settings()
#
#
#     #after everything minimization should be empty -> for the next iteration
#     if not charmm_pac.charmm_instructions['minimize']:
#         charmm_pac.charmm_instructions['minimize'] = []
#     else:
#         error = "Minimization instructions are not empty!"
#         raise AssertionError(error)
#
#     return
#
# def pac_ph42(workdir, charmm_pac):
#
#     '''David'''
#
#     charmm_pac.backup_settings()
#     if not os.path.exists(workdir):
#         os.mkdir(workdir)
#     charmm_pac.workdir = workdir
#
#     charmm_pac.backup_settings()
#     charmm_pac.charmm_instructions['minimize'] = []
#     charmm_pac.charmm_instructions['do_minimize'] = False
#
#     ##########################
#     ### CHANGE PROTONATION ###
#     ##########################
#
#     #placeholder
#     residue_tuple = ('DPP', 21, 'A')
#     state = 1
#     charmm_pac.set_prot_residue(residue_tuple, state=state)
#
#     spec_command = get_sel_commands()
#
#     charmm_pac.add_charmm_command(spec_command, adj_task='hbuild')
#     charmm_pac.run_charmm()
#     charmm_pac.update_coords_with_modelled_structure()
#     charmm_pac.restore_settings()
#
#
#     #after everything minimization should be empty -> for the next iteration
#     if not charmm_pac.charmm_instructions['minimize']:
#         charmm_pac.charmm_instructions['minimize'] = []
#     else:
#         error = "Minimization instructions are not empty!"
#         raise AssertionError(error)
#
#     return
#
# def pac_ph43(workdir, charmm_pac):
#
#     '''David'''
#
#     charmm_pac.backup_settings()
#     if not os.path.exists(workdir):
#         os.mkdir(workdir)
#     charmm_pac.workdir = workdir
#
#     charmm_pac.backup_settings()
#     charmm_pac.charmm_instructions['minimize'] = []
#     charmm_pac.charmm_instructions['do_minimize'] = False
#
#     ##########################
#     ### CHANGE PROTONATION ###
#     ##########################
#
#     #placeholder
#     residue_tuple = ('DPP', 40, 'A')
#     state = 0
#     charmm_pac.set_prot_residue(residue_tuple, state=state)
#
#     spec_command = get_sel_commands()
#
#     charmm_pac.add_charmm_command(spec_command, adj_task='hbuild')
#     charmm_pac.run_charmm()
#     charmm_pac.update_coords_with_modelled_structure()
#     charmm_pac.restore_settings()
#
#
#     #after everything minimization should be empty -> for the next iteration
#     if not charmm_pac.charmm_instructions['minimize']:
#         charmm_pac.charmm_instructions['minimize'] = []
#     else:
#         error = "Minimization instructions are not empty!"
#         raise AssertionError(error)
#
#     return
#
# def pac_ph44(workdir, charmm_pac):
#
#     '''David'''
#
#     charmm_pac.backup_settings()
#     if not os.path.exists(workdir):
#         os.mkdir(workdir)
#     charmm_pac.workdir = workdir
#
#     charmm_pac.backup_settings()
#     charmm_pac.charmm_instructions['minimize'] = []
#     charmm_pac.charmm_instructions['do_minimize'] = False
#
#     ##########################
#     ### CHANGE PROTONATION ###
#     ##########################
#
#     #placeholder
#     residue_tuple = ('EPP', 43, 'A')
#     state = 1
#     charmm_pac.set_prot_residue(residue_tuple, state=state)
#
#     spec_command = get_sel_commands()
#
#     charmm_pac.add_charmm_command(spec_command, adj_task='hbuild')
#     charmm_pac.run_charmm()
#     charmm_pac.update_coords_with_modelled_structure()
#     charmm_pac.restore_settings()
#
#
#     #after everything minimization should be empty -> for the next iteration
#     if not charmm_pac.charmm_instructions['minimize']:
#         charmm_pac.charmm_instructions['minimize'] = []
#     else:
#         error = "Minimization instructions are not empty!"
#         raise AssertionError(error)
#
#     return
#
#
# def get_sel_commands():
#
#     spec_command = """
#
#
# !### Definition of GB ###
# stream "/scratch/scratch/tmeyer/CHARMM_NAMD/toppar_36/gbsw/radius_gbsw.str"
# scalar wmain statistics select .not. type H* end
# define check select (.not type H* ) .and. ( prop wmain .eq. 0.0 ) show end
# if ?nsel ne 0  stop       !some heavy atom have a zero radius
# GBSW conc 0.1 sgamma 0.03 GBenergy epsp 4 MOLSURF
#
#
# !define sel1 select .byres. ((segid A .and. resid 21) .around. 5) end
# !define bb select (type N .or. type C .or. type CA .or. type O) end
# !cons fix select .not. (sel1 .or. hydrogen) .or. bb end
#
# !define sel1 select (segid A .and. resid 21) .or. (segid A .and. resid 19) end
# !define bb select (type N .or. type C .or. type CA .or. type O) end
# !cons fix select .not. (sel1 .or. hydrogen) .or. bb end
#
# define sel1 select (segid A .and. (resid 21 .or. resid 19 .or. resid 40 .or. resid 43 .or. resid 35)) end
# define bb select (type N .or. type C .or. type CA .or. type O) end
# cons fix select .not. (sel1 .or. hydrogen) .or. bb end
#
#
# minimize sd nsteps 1000 tolg 0.1
# minimize abnr nsteps 5000 tolg 0.01
#
#
# """
#
#     return spec_command
