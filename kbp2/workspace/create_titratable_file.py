# -*- coding: utf-8 -*-
__author__ = 'Tim Meyer'


import sys
sys.path.append('/user/tmeyer/workspace/script/protein_toolbox')


from kbp2 import charmm


class titratable_residues(object):
    def __init__(self, top):
        """
            top : Single string or list of strings containing path to topology files.
        """
        # Structure of entries in residues: {<name> : [{'type' : R,
        #                                              'pka'   : 0.0,
        #                                              'patch' : None,
        #                                              'atoms' : {<name> : <charge>}
        #                                              'order' : [<atom names>]
        #                                              },
        #                                              ...]
        #                                   }
        self.residues = {}

        if type(top) == str:
            top = [top]
        # The charmm object is used to parse topology files.
        self.charmm = charmm.Charmm_manager(top=top, par=[''])

    def add_state(self, name, type, pka, patch=None, user_atom=None, nocharge=False):

        # Prepare new state.
        state = {'type'  : type,
                 'pka'   : pka,
                 'patch' : patch,
                 'atoms' : None,
                 'order' : []}

        atoms = {}
        order = []
        if user_atom is not None:
            ### If the user specifies atoms and their charges, they are used. ###
            for atom_name, atom_charge in user_atom.iteritems():
                if nocharge:
                    atoms[atom_name] = 0.0
                else:
                    atoms[atom_name] = atom_charge
                order.append(atom_name)
        else:
            ### Get charges from topology file and add them. ###
            # Copy atoms and their charge from 'name' entry in the topology files.
            atoms = {}
            order = []
            for atom in self.charmm.top_content['residues'][name]['atoms']:
                atom_name   = atom['name']
                atom_charge = atom['charge']
                if nocharge:
                    atom_charge = 0.0
                atoms[atom_name] = float(atom_charge)
                order.append(atom_name)
        if patch is not None:
            # Overwrite those atom charges that are specified in the patch.
            for atom in self.charmm.top_content['patches'][patch]['atoms']:
                atom_name   = atom['name']
                atom_charge = atom['charge']
                if nocharge:
                    atom_charge = 0.0
                atoms[atom_name] = float(atom_charge)
                # Add atom if they appear in the patch. This is not supported by Karlsberg+.
                # if not atom_name in order:
                #     order.append(atom_name)
                #     error = "The patch %s is adding an atom. This is not supported by Karlsberg+!"
                #     raise AssertionError(error)
        state['atoms'] = atoms
        state['order'] = order

        if self.residues.has_key(name):
            # Add state to existing residue.

            ### Checks: ###
            # Only one 'R' state is allowed.
            if type == 'R':
                for other_state in self.residues[name]:
                    if other_state['type'] == 'R':
                        error = 'State %s:%s can not be added. Only one reference state (type=R) is allowed!' \
                                % (name, type)
                        raise AssertionError(error)

            self.residues[name].append(state)
        else:
            # Create new residue and add state.
            self.residues[name] = [state]


    def purge_constant_charges(self, residue_names=None, exclude_ref=False):
        """
        Removes all atom entries for the specified residues, that are not changed within at least one state.

        residue_names : Single string, list of strings or None. Strings are the names of those residues
                        that are processed.
        """
        if residue_names is None:
            resnames = self.residues.keys()
        elif type(residue_names) is str:
            resnames = [residue_names]
        else:
            resnames = residue_names


        for resname in resnames:
            # Find atoms that do change.
            all_atoms     = {}
            changed_atoms = {}
            for state in self.residues[resname]:
                if exclude_ref:
                    if state['type'] == 'R':
                        continue
                for atom_name, atom_charge in state['atoms'].iteritems():
                    # atom_name   = atom['name']
                    # atom_charge = atom['charge']
                    if all_atoms.has_key(atom_name):
                        if abs(atom_charge - all_atoms[atom_name]) > 0.0001:
                            changed_atoms[atom_name] = True
                    else:
                        all_atoms[atom_name] = atom_charge
            # Delete all atoms, that are not changed.
            for state in self.residues[resname]:
                atoms = {}
                for atom_name, atom_charge in state['atoms'].iteritems():
                    # atom_name = atom['name']
                    if changed_atoms.has_key(atom_name):
                        # Copy the atom.
                        atoms[atom_name] = atom_charge
                # Overwrite atom list with filtered one.
                state['atoms'] = atoms

    def remove_atoms(self, atom_names):
        """
        Removes all atom entries for the specified in all residues.

        atom_names : List of strings. Strings are the atom names.
        """

        for resname in self.residues.keys():
            # Find atoms that do change.
            all_atoms     = {}
            kept_atoms = {}
            for state in self.residues[resname]:
                for atom_name, atom_charge in state['atoms'].iteritems():
                    # atom_name   = atom['name']
                    # atom_charge = atom['charge']
                    if all_atoms.has_key(atom_name):
                        if atom_name not in atom_names:
                            kept_atoms[atom_name] = True
                    else:
                        all_atoms[atom_name] = atom_charge
            # Delete all atoms, that are not changed.
            for state in self.residues[resname]:
                atoms = {}
                for atom_name, atom_charge in state['atoms'].iteritems():
                    # atom_name = atom['name']
                    if kept_atoms.has_key(atom_name):
                        # Copy the atom.
                        atoms[atom_name] = atom_charge
                # Overwrite atom list with filtered one.
                state['atoms'] = atoms

    def write_kbp_titratable_file(self, filename):
        text = '---\n'

        for resname in self.residues.keys():
            text += '%s:\n' % resname
            for state in self.residues[resname]:
                text += ' - sym: %s\n' % state['type']
                text += '   pka: %.2f\n' % state['pka']
                text += '   apc:\n'
                for atom_name in state['order']:
                    if not state['atoms'].has_key(atom_name):
                        continue
                    atom_charge = state['atoms'][atom_name]
                    charge_str = '%.3f' % atom_charge
                    if charge_str[-1] == '0':
                        # text += '    %s: %5.2f\n' % (atom_name, atom_charge)
                        text += '    %s: %.2f\n' % (atom_name, atom_charge)
                    else:
                        # text += '    %s: %6.3f\n' % (atom_name, atom_charge)
                        text += '    %s: %.3f\n' % (atom_name, atom_charge)
                if state['patch'] is not None:
                    text += '   patch: %s\n' % state['patch']

        f = open(filename, 'w')
        f.write(text)
        f.close()


if __name__ == '__main__':
    top = []
    # top.append("/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/top.inp")
    # top.append("/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/patches.rtf")
    top.append("/scratch/scratch/tmeyer/karlsbergplus/top.inp")
    top.append("/scratch/scratch/tmeyer/karlsbergplus/patches.rtf")

    # par = []
    # par.append("/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/par.inp")
    # par.append("/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/patches.prm")


    titr_res = titratable_residues(top)

    dummyRef = True

    if dummyRef:
        ref_state_pka = 100.0
        # ref_state_pka = 0.0
    else:
        ref_state_pka = 0.0

    if dummyRef:
        titr_res.add_state('DPP', 'R', 0.0, 'DPPN', nocharge=True)
        titr_res.add_state('DPP', '0', 0.0-ref_state_pka, 'DPP1')
    else:
        titr_res.add_state('DPP', 'R', 0.0, 'DPP1')
    titr_res.add_state('DPP', 'D', 4.0-ref_state_pka, None)
    # PROPKA:
    # titr_res.add_state('DPP', 'D', 3.8-ref_state_pka, None)
    # Fit:
    # titr_res.add_state('DPP', 'D', 3.85-ref_state_pka, None)
    titr_res.add_state('DPP', '0', 0.0-ref_state_pka, 'DPP2')


    if dummyRef:
        titr_res.add_state('EPP', 'R', 0.0, 'EPPN', nocharge=True)
        titr_res.add_state('EPP', '0', 0.0-ref_state_pka, 'EPP1')
    else:
        titr_res.add_state('EPP', 'R', 0.0, 'EPP1')
    titr_res.add_state('EPP', 'D', 4.4-ref_state_pka, None)
    # PROPKA:
    # titr_res.add_state('EPP', 'D', 4.5-ref_state_pka, None)
    # Fit:
    # titr_res.add_state('EPP', 'D', 4.6-ref_state_pka, None)
    titr_res.add_state('EPP', '0', 0.0-ref_state_pka, 'EPP2')


    if dummyRef:
        titr_res.add_state('HSP', 'R', 0.0, 'HSPN', nocharge=True)
        titr_res.add_state('HSP', '0',  0.0-ref_state_pka, 'HSPD')
    else:
        titr_res.add_state('HSP', 'R',  0.0, 'HSPD')
    titr_res.add_state('HSP', 'P', -7.0-ref_state_pka, None)
    # PROPKA:
    # titr_res.add_state('HSP', 'P', -6.5-ref_state_pka, None)
    # Fit:
    # titr_res.add_state('HSP', 'P', -7.4-ref_state_pka, None)
    titr_res.add_state('HSP', '0', -0.4-ref_state_pka, 'HSPE')

    if dummyRef:
        titr_res.add_state('ARG', 'R', 0.0, 'ARGN', nocharge=True)
        titr_res.add_state('ARG', '0',  0.0-ref_state_pka, 'ARGR')
    else:
        titr_res.add_state('ARG', 'R',  0.0, 'ARGR')
    titr_res.add_state('ARG', 'P', -12.0-ref_state_pka, None)
    # PROPKA:
    # titr_res.add_state('ARG', 'P', -12.0-ref_state_pka, None)
    # Fit: (from PROPKA)
    # titr_res.add_state('ARG', 'P', -12.0-ref_state_pka, None)

    if dummyRef:
        titr_res.add_state('LYS', 'R', 0.0, 'LYSN', nocharge=True)
        titr_res.add_state('LYS', '0',   0.0-ref_state_pka, 'LYSR')
    else:
        titr_res.add_state('LYS', 'R',   0.0, 'LYSR')
    titr_res.add_state('LYS', 'P', -10.4-ref_state_pka, None)
    # PROPKA:
    # titr_res.add_state('LYS', 'P', -10.5-ref_state_pka, None)
    # Fit:
    # titr_res.add_state('LYS', 'P', -10.3-ref_state_pka, None)

    if dummyRef:
        titr_res.add_state('TYR', 'R', 0.0, 'TYRN', nocharge=True)
        titr_res.add_state('TYR', '0', 0.0-ref_state_pka, None)
    else:
        titr_res.add_state('TYR', 'R', 0.0, None)
    titr_res.add_state('TYR', 'D', 9.6-ref_state_pka, 'TYRD')
    # PROPKA:
    # titr_res.add_state('TYR', 'D', 10.0-ref_state_pka, 'TYRD')
    # Fit:
    # titr_res.add_state('TYR', 'D', 9.5-ref_state_pka, 'TYRD')

    if dummyRef:
        titr_res.add_state('CYS', 'R', 0.0, 'CYSN', nocharge=True)
        titr_res.add_state('CYS', '0', 0.0-ref_state_pka, None)
    else:
        titr_res.add_state('CYS', 'R', 0.0, None)
    titr_res.add_state('CYS', 'D', 9.5-ref_state_pka, 'CYSD')
    # PROPKA:
    # titr_res.add_state('CYS', 'D', 9.0-ref_state_pka, 'CYSD')

    # titr_res.purge_constant_charges(exclude_ref=True)


    atoms = {'N' : -0.97,
             'HT1' : 0.22,
             'HT2': 0.22,
             'HT3' : 0.22,
             'CA' : 0.21,
             'C' : 0.51,
             'O' : -0.51}
             # 'HA' : 0.09}
    # atoms = {'N' : -0.97,
    #          'HT1' : 0.22,
    #          'HT2': 0.22,
    #          'HT3' : 0.22}
    if dummyRef:
        titr_res.add_state('NTE', 'R', 0.0, 'NTEN', atoms, nocharge=True)
        titr_res.add_state('NTE', '0',  0.0-ref_state_pka, 'NTEREF', atoms)
    else:
        titr_res.add_state('NTE', 'R',  0.0, 'NTEREF', atoms)
        # titr_res.add_state('NTE', '0',  0.0, 'NTEREF', atoms)

    atoms = {'N' : -0.30,
             'HT1' : 0.33,
             'HT2': 0.33,
             'HT3' : 0.33,
             'CA' : 0.21,
             'C' : 0.51,
             'O' : -0.51}
             # 'HA' : 0.10}
    titr_res.add_state('NTE', 'P', -7.5-ref_state_pka, None, atoms)
    # PROPKA:
    # titr_res.add_state('NTE', 'P', -8.0-ref_state_pka, None, atoms)
    # Fit:
    # titr_res.add_state('NTE', 'P', -8.4-ref_state_pka, None, atoms)



    atoms = {'C' : 0.34,
             'OT1' : -0.17,
             'OT2' : -0.17,
             'CA' : 0.07,
             'N' : -0.47,
             'HN' : 0.31}
             # 'HA' : 0.09}
    if dummyRef:
        titr_res.add_state('CTE', 'R', 0.0, 'CTEN', atoms, nocharge=True)
        titr_res.add_state('CTE', '0', 0.0-ref_state_pka, 'CTEREF', atoms)
    else:
        titr_res.add_state('CTE', 'R', 0.0, 'CTEREF', atoms)
        # titr_res.add_state('CTE', '0', 0.0, 'CTEREF', atoms)

    atoms = {'C' : 0.34,
             'OT1' : -0.67,
             'OT2' : -0.67,
             'CA' : 0.07,
             'N' : -0.47,
             'HN' : 0.31}
             # 'HA' : 0.09}
    titr_res.add_state('CTE', 'D', 3.8-ref_state_pka, None, atoms)
    # PROPKA:
    # titr_res.add_state('CTE', 'D', 3.2-ref_state_pka, None, atoms)
    # Fit:
    # titr_res.add_state('CTE', 'D', 3.6-ref_state_pka, None, atoms)

    # Remove backbone
    # titr_res.remove_atoms(['N', 'HN', 'CA', 'HA', 'C', 'O'])


    # Set HA or first atom to 0.001.
    if dummyRef:
        for name, residue in titr_res.residues.iteritems():
            first_atom = None
            for atom_name in residue[0]['order']:
                if residue[0]['atoms'].has_key(atom_name):
                    first_atom = atom_name
                    break
            if residue[0]['atoms'].has_key('HA') and name != 'NTE':
                residue[0]['atoms']['HA'] = 0.001
            else:
                residue[0]['atoms'][first_atom] = 0.001


    # titr_res.write_kbp_titratable_file('/user/tmeyer/titratable.yaml')
    # titr_res.write_kbp_titratable_file('/scratch/scratch/tmeyer/md_pka/runs/kbp_cache/2lzt/complete_titratable/titratable.yaml')

    # titr_res.write_kbp_titratable_file('/scratch/scratch/tmeyer/md_pka/runs/titratables/titratable_refsmall_shifted.yaml')
    # titr_res.write_kbp_titratable_file('/scratch/scratch/tmeyer/md_pka/runs/titratables/titratable_ter.yaml')

    titr_res.write_kbp_titratable_file('/scratch/scratch/tmeyer/md_pka/runs/titratables/titratable_refsmall_shifted_ter2_np.yaml')
    # titr_res.write_kbp_titratable_file('/scratch/scratch/tmeyer/md_pka/runs/titratables/titratable_refsmall_shifted_ter2_np_propka.yaml')
    # titr_res.write_kbp_titratable_file('/scratch/scratch/tmeyer/md_pka/runs/titratables/titratable_refsmall_shifted_ter2_np_fit.yaml')

    # titr_res.write_kbp_titratable_file('/scratch/scratch/tmeyer/md_pka/runs/kbp_cache/2lzt/ref_zero/titratable.yaml')
    # titr_res.write_kbp_titratable_file('/scratch/scratch/tmeyer/md_pka/runs/kbp_cache/2lzt/complete_titratable/titratable.yaml')
    # titr_res.write_kbp_titratable_file('/scratch/scratch/tmeyer/md_pka/runs/kbp_cache/2lzt/ref_zero_Rpkat/titratable.yaml')
    # titr_res.write_kbp_titratable_file('/scratch/scratch/tmeyer/md_pka/runs/kbp_cache/2lzt/ref_shifted/titratable.yaml')
    # titr_res.write_kbp_titratable_file('/scratch/scratch/tmeyer/md_pka/runs/titratables/titratable_refsmall_shifted_sidechains.yaml')
    # titr_res.write_kbp_titratable_file('/scratch/scratch/tmeyer/kbplus2_enere/titratable_refsmall_shifted_sidechains.yaml')
    # titr_res.write_kbp_titratable_file('/user/tmeyer/anke/titratable_refZero.yaml')


    #
    # for name, res in titr_res.residues.iteritems():
    #     print name
    #     for state in res:
    #         print state
    #     print