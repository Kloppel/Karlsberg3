# coding=utf-8
from collections import defaultdict
import re
import subprocess
import numpy as np

import sys,os,shutil

from file_parser import Simple_struct_parser

# sys.path.append('/user/tmeyer/workspace/script/protein_toolbox/cython')
# from all_vs_all_atoms import find_collisions

from cython.all_vs_all_atoms import find_collisions

import kbp_tools

def is_std_aa(name):
    """
    Returns 'one letter' code for each residue
    """
    std_aa = {'ARG' : 'R',\
              'HIS' : 'H',\
              'HSP' : 'H',\
              'HSD' : 'H',\
              'HSE' : 'H',\
              'LYS' : 'K',\
              'ASP' : 'D',\
              'GLU' : 'E',\
              'SER' : 'S',\
              'THR' : 'T',\
              'ASN' : 'N',\
              'GLN' : 'Q',\
              'CYS' : 'C',\
              'SEC' : 'U',\
              'GLY' : 'G',\
              'PRO' : 'P',\
              'ALA' : 'A',\
              'ILE' : 'I',\
              'LEU' : 'L',\
              'MET' : 'M',\
              'PHE' : 'F',\
              'TRP' : 'W',\
              'TYR' : 'Y',\
              'VAL' : 'V',\
              'EPP' : 'E',\
              'DPP' : 'D'\
            }

    if std_aa.has_key(name):
        return 1
    else:
        return 0

def is_na_base(name):
    """
    Returns 'one letter' code for each residue
    """
    std_base = {'ADE' : 'A',\
              'GUA' : 'G',\
              'CYT' : 'C',\
              'THY' : 'T',\
              'URA' : 'U',\
            }

    if std_base.has_key(name):
        return 1
    else:
        return 0


def get_std_titr_residues():
    """
    Returns a dictionary that contains information about standard titratable residues. The entries for 'pka',
    'atoms' and 'order' are not set.

    Returns a dict with the following structure:
    {<name> : [dict: <state parameters>]
    }

    with <state parameters>:
    <state parameters> = {'charge' : int: <total charge of the residue>,
                       'patch' : None|str: <name of the CHARMM patch required to create this state>,
                       'rename' : None|str: <new residue name>,
                       'external_patches' : [ (str: <CHARMM patch name>, list of (<resname>, <resid>, <segname>)) ],
                       'special' : None|str: <description used later to treat the residue in a special way>
                      }
    Remarks:
    - It is not allowed to set 'patch' and 'rename' at the same time.
    - 'external_patches': This keyword allows to specify patches that apply to other residues. Each list
    entry is a tuple with two entries, first the name if the patch and second a list of residues, that are
    used as parameter in the CHARMM patch statement. The residues are described are tuples in the format
    (<resid>, <segname>)
    e.g. The entry ['CRES' : [(12, 'A'), (34, 'B')]] will create to the following CHARMM command:
        PATCH CRES A 12 B 34
    - 'special': If set, this entry will cause a special treatment of the reside. Allowed keywords are:
              'nter': The residue is a N-terminal residue.
              'cter': The residue is a C-terminal residue.
    # @return: dict
    """
    template_state = {'charge' : 0,
                      'patch' : None,
                      'external_patches' : None,
                      'rename' : None,
                      'special' : None}

    titr_residues = {}

    resname = 'ASP'
    titr_residues[resname] = [dict(template_state) for i in range(2)]
    titr_residues[resname][0]['charge'] = -1
    titr_residues[resname][1]['charge'] = 0
    titr_residues[resname][1]['patch'] = 'ASPP'

    resname = 'GLU'
    titr_residues[resname] = [dict(template_state) for i in range(2)]
    titr_residues[resname][0]['charge'] = -1
    titr_residues[resname][1]['charge'] = 0
    titr_residues[resname][1]['patch'] = 'GLUP'

    resname = 'ARG'
    titr_residues[resname] = [dict(template_state) for i in range(2)]
    titr_residues[resname][0]['charge'] = 1
    titr_residues[resname][1]['charge'] = 0
    titr_residues[resname][1]['patch'] = 'ARGR'

    resname = 'LYS'
    titr_residues[resname] = [dict(template_state) for i in range(2)]
    titr_residues[resname][0]['charge'] = 1
    titr_residues[resname][1]['charge'] = 0
    titr_residues[resname][1]['patch'] = 'LYSR'

    resname = 'HIS'
    titr_residues[resname] = [dict(template_state) for i in range(3)]
    titr_residues[resname][0]['charge'] = 0
    titr_residues[resname][0]['rename'] = 'HSE'
    titr_residues[resname][1]['charge'] = 0
    titr_residues[resname][1]['rename'] = 'HSD'
    titr_residues[resname][2]['charge'] = 1
    titr_residues[resname][2]['rename'] = 'HSP'

    resname = 'TYR'
    titr_residues[resname] = [dict(template_state) for i in range(2)]
    titr_residues[resname][0]['charge'] = 0
    titr_residues[resname][1]['charge'] = -1
    titr_residues[resname][1]['patch'] = 'TYRD'

    resname = 'CYS'
    titr_residues[resname] = [dict(template_state) for i in range(2)]
    titr_residues[resname][0]['charge'] = 0
    titr_residues[resname][1]['charge'] = -1
    titr_residues[resname][1]['patch'] = 'CYSD'

    resname = 'NTE'
    titr_residues[resname] = [dict(template_state) for i in range(2)]
    titr_residues[resname][0]['charge'] = 1
    titr_residues[resname][0]['special'] = 'nter'
    titr_residues[resname][1]['charge'] = 0
    titr_residues[resname][1]['patch'] = 'NTEREF'
    titr_residues[resname][1]['special'] = 'nter'

    resname = 'CTE'
    titr_residues[resname] = [dict(template_state) for i in range(2)]
    titr_residues[resname][0]['charge'] = -1
    titr_residues[resname][0]['special'] = 'cter'
    titr_residues[resname][1]['charge'] = 0
    titr_residues[resname][1]['patch'] = 'CTEREF'
    titr_residues[resname][1]['special'] = 'cter'

    return titr_residues

# KNOWN ISSUES OF Charmm_manager CLASS:
#             1. In case of disulfide bridges, all must be either opened or closed
#             2. Does not carefully deal with Cys residues that are a part of a disulfide bridge
#                in case of modifying their charge (application of two patches in CHARMM)

class Charmm_manager:

    class TopologyError(Exception):
            def __init__(self, value):
                self.value = value
            def __str__(self):
                return repr(self.value)
    class GeneralError(Exception):
            def __init__(self, value):
                self.value = value
            def __str__(self):
                return repr(self.value)


    def __init__(self, workdir=None, pdb='',  charmm_bin='', top=[], par=[]):
        """
        workdir:    The folder where all files will be created and where Charmm will be executed.
        charmm_bin: The Charmm binary. Default: "charmm36"
        top:        One or more topology files. Can be a string or a list of strings.
                    Default: "/scratch/scratch/tmeyer/CHARMM_NAMD/toppar_36/top_all27_prot_na.rtf"
        par:        One or more parameter files. Can be a string or a list of strings.
                    Default: "/scratch/scratch/tmeyer/CHARMM_NAMD/toppar_36/par_all27_prot_na.prm"
        """

        self.top = []
        if not top:
            error = "Topology file is not provided"
            # raise AssertionError(error)
        else:
            if isinstance(top, str):
                top = [top]
            self.top = list(top)

        self.par = []
        if not par:
            error = "Parameter file is not provided"
            # raise AssertionError(error)
        else:
            if isinstance(par, str):
                par = [par]
            self.par = list(par)

        if not charmm_bin:
            self.charmm_bin = "charmm36b1_64"
        else:
            self.charmm_bin = charmm_bin

        if workdir is None:
            self.workdir = "/tmp/charmm_manager/"
        else:
            if workdir[-1] != '/':
                workdir += '/'
            self.workdir = workdir
        if not os.path.exists(self.workdir):
            os.mkdir(self.workdir)

        # Everything that has been done is commented here.
        self.log = []

        self.stream = ''

        # Will store input structure as Simple_struct_parser object.
        self.structure = None

        # Todo: delete this thing, it seems it is not used
        # Prefix for all files generated.
        # self.prefix = None

        # Status flags.
        self.structure_read      = False
        self.structures_checked  = False
        self.structures_prepared = False

        self.backuped_settings = None

        # Container for charmm output. Both .out and structure files.
        self.generated_files = {}

        # Stores profiles for charmm runs.
        self.charmm_script_profiles = {}


        ### Stores the content of the top files, as far as it is parsed. ###
        self.top_content = {}
        self.top_content['residues'] = {}
        # Structure of 'residues': {name : {name: <name>,
        #                                   charge: <charge>,
        #                                   atoms: [{name: <name>,
        #                                           type: <type>,
        #                                           charge: <charge>}]
        #                                   }}
        self.top_content['patches'] = {}
        # Structure of 'patches': {name : {name: <name>,
        #                                   charge: <charge>,
        #                                   atoms: [{name: <name>,
        #                                           type: <type>,
        #                                           charge: <charge>}]
        #                                   }}


        ### Decision strings that define user decision are stored here. ###
        self.decisions = {}
        # Descriptions of decisions that are decided, by the user or by using the default.
        self.decided   = []
        # Descriptions of decisions that have to be made.
        self.undecided = []


        ### Container for the consequences of decisions. ###
        self.charmm_instructions = {}

        # Actions concerting gaps:
        #  content: list of ('connect', start residue, end residue) or
        #                   ('fill', start residue, end residue, sequence string in three letter code)
        self.charmm_instructions['gaps'] = []

        #  content: list of (start residue, end residue) or residue. These residues will be minimized. residue are
        #  SimpleAtomAccessResidue objects.
        self.charmm_instructions['minimize'] = []

        # Decision whether a energy minimisation should be performed.
        # Only atoms from "self.charmm_instructions['minimize']" and all hydrogens are minimized.
        # content: False/True
        self.charmm_instructions['do_minimize'] = False

        self.charmm_instructions['backbone_fixed'] = False

        # Actions concerning termini:
        #  content: 'segname' => (nter residue, cter residue)
        self.charmm_instructions['ter'] = {}
        # List of segments, where the termini are patched manually, instead of via 'generate'.
        # Todo: Add better description
        self.charmm_instructions['patch_ter'] = {}
        # Caping of termini that are not a true teminus
        self.charmm_instructions['cap_termini'] = None
        # Extra functionality to completely ignore all termini -> should be used with care!!!
        self.charmm_instructions['ignore_termini'] = False

        # Actions concerning disulphide bridges:
        #  content: list of (first cystein residue, second cysteine residue)
        self.charmm_instructions['disu'] = []

        # Actions concerning missing atoms:
        #  content: string -> model or keep
        self.charmm_instructions['missing'] = 'model'

        # Patches that should be applied after the complete structure except for water chains is generated and before
        # coordinates are read. If any patch is defined, AUTOGGYATE ANGLE DIHEDRAL will be call afterwards.
        #  content: list of (str: <patch name>, list of (<segname>, <resid>, <atom-name>))
        self.charmm_instructions['patches'] = []

        # These patches that should be applied after AUTOGENERATE ANGLE DIHEDRAL command
        #  content: list of (str: <patch name>, list of  (<segname>, <resid>, <atom-name>))
        self.charmm_instructions['patches_no_autogen'] = []

        ### description
        self.charmm_instructions['minimize_selections'] = []

        ### The Charmm configuration. ###
        self.charmm_config = {}

        # Output files:
        self.charmm_config['output'] = ['pdb', 'crd', 'psf', 'xplor.psf']

        # IDs of missing atom:
        # Structure: [ (<segname>, <resid>, <atom-name>) ]
        self.charmm_config['missing_atom_ids'] = []

        # Should AUTOgenerate ANGLes DIHEdrals be called before the water chains are read?
        self.charmm_config['autogen'] = False

        # backup object without self.structure
        self.charmm_config['backup'] = None


        ### Task list for Charmm:
        self.charmm_config['tasks'] = []
        # Define sequence and read coordinates for not water chains from PDBs.
        self.charmm_config['tasks'].append('read_structure_sequence')
        # Model missing residues if requested in self.charmm_instructions['gaps'].
        self.charmm_config['tasks'].append('model_gaps')
        # Patch disulphide bridges mentioned in self.charmm_instructions['disu'].
        self.charmm_config['tasks'].append('patch_disu')
        # Patch termini manually if requested.
        self.charmm_config['tasks'].append('patch_ter')
        # Apply patches, that have been specified in self.charmm_instructions['patches']
        self.charmm_config['tasks'].append('patches')
        # Generate angles and dihedrals if necessary.
        self.charmm_config['tasks'].append('autogen')
        # Apply patches afetr AUTOGENERATE ANGLE DIHEDRAL
        self.charmm_config['tasks'].append('patches_no_autogen')
        # Define sequence and read coordinates for water chains from PDBs.
        self.charmm_config['tasks'].append('read_water_sequence')
        # Read the coordinates for non water chains.
        self.charmm_config['tasks'].append('read_coordinates')
        # Read the coordinates for water chains.
        self.charmm_config['tasks'].append('read_water_coordinates')
        # Model missing atoms if requested in self.charmm_instructions['missing'] using IC PARA / IC BUILD.
        self.charmm_config['tasks'].append('build_missing')
        # Run HBUILD command.
        self.charmm_config['tasks'].append('hbuild')
        # Minimize modeled atoms and residues (and hydrogens) if necessary.
        self.charmm_config['tasks'].append('minimize_modeled')
        # Write out the files defined in self.charmm_config['output'].
        self.charmm_config['tasks'].append('write_structure')

        # If a list entry in 'tasks' does not match with any key from above it is interpreted as

        # These instructions will be included in the input script.


        # A string containing the charmm output.
        self.charmm_output = ''

        # Will be set to True, if a charmm run has been completed successfully.
        self.charmm_normal_termination = False

        # Set by add_structure to the filename of the read structure without its suffix e.g. '2lzt' for 2lzt.pdb
        self.title = None

        # If this variable is not initialized, it is set to '<self.title>_out' when the CHARMM script is generated.
        # It is used as a prefix for the CHARMM output files.
        self.charmm_out_prefix = None

        # If True, the structures written out for charmm are removed after the successful CHARMM run.
        self.clean_up = False

        # Contains the names of files that are temporarily generated for charmm. They will be deleted after
        # running charmm. This does not include th input script!
        self.charmm_tmp_files = []

        # Contains a dict that contains information of what type of resides can change their protonation and
        # instructions of how this can be done in terms of CHARMM patch statements or renaming. It is set to default
        # values below.
        self.titr_residues = {}
        # This list is needed to find the corresponding entry in self.titr_residues, if the residue has been renamed.
        self.titr_restypes = {}
        # This dictionary stores the information of what protonation state is currently assigned to residues.
        # The structure is {(str: <resname>, int: <resid>, str: <segname>) : int: <state nr according to entry nr in self.titr_residues>}
        self.prot_state = {}


        ### First actions ###
        # Set self.titr_residues to default values.
        titr_residues = get_std_titr_residues()
        self.set_titr_residues(titr_residues)

        if top:
            self.parse_top()

        if pdb:
            self.add_structure(pdb)

    def copy(self, no_structure=False, no_backup=False):
        """
        Makes a copy of a Charmm_manager object
        """

        new = Charmm_manager()

        new.top = self.top
        new.par = self.par

        new.charmm_bin = self.charmm_bin

        new.workdir = self.workdir

        new.log = []
        new.log = list(self.log)

        new.stream = self.stream

        if no_backup:
            new.backuped_settings = None
        else:
            new.backuped_settings = self.backuped_settings.copy(no_structure=True, no_backup=True)

        if not no_structure:
            new.structure = self.structure.copy()
        else:
            new.structure = None


        new.structure_read = self.structure_read
        new.structures_checked = self.structures_checked
        new.structures_prepared = self.structures_prepared

        new.generated_files = {}
        new.generated_files = dict(self.generated_files)

        new.charmm_script_profiles = {}
        new.charmm_script_profiles = dict(self.charmm_script_profiles)

        new.top_content = {}

        new.top_content['residues'] = {}
        for name, properties in self.top_content['residues'].iteritems():
            new.top_content['residues'][name] = {}
            new.top_content['residues'][name]['name'] = properties['name']
            new.top_content['residues'][name]['charge'] = properties['charge']
            new.top_content['residues'][name]['atoms'] = []
            for atom in properties['atoms']:
                new.top_content['residues'][name]['atoms'].append(dict(atom))

        new.top_content['patches'] = {}
        for name, properties in self.top_content['patches'].iteritems():
            new.top_content['patches'][name] = {}
            new.top_content['patches'][name]['name'] = properties['name']
            new.top_content['patches'][name]['charge'] = properties['charge']
            new.top_content['patches'][name]['atoms'] = []
            for atom in properties['atoms']:
                new.top_content['patches'][name]['atoms'].append(dict(atom))

        new.decisions = {}
        for choice, options in self.decisions.iteritems():
            new.decisions[choice] = tuple(options)

        new.decided = []
        new.decided = list(self.decided)

        new.undecided = []
        new.undecided = list(self.undecided)

        new.charmm_instructions = {}

        new.charmm_instructions['gaps'] = []
        for gap in self.charmm_instructions['gaps']:
            new.charmm_instructions['gaps'].append(tuple(gap))

        new.charmm_instructions['minimize'] = []
        new.charmm_instructions['minimize'] = list(self.charmm_instructions['minimize'])

        new.charmm_instructions['do_minimize'] = self.charmm_instructions['do_minimize']

        new.charmm_instructions['backbone_fixed'] = self.charmm_instructions['backbone_fixed']

        new.charmm_instructions['ter'] = {}
        for chain, residues in self.charmm_instructions['ter'].iteritems():
            residue1, residue2 = residues
            if residue1 is not None and residue2 is not None:
                res1 = residue1.copy()
                res2 = residue2.copy()
            else:
                res1 = residue1
                res2 = residue2
            termini_tuple = (res1, res2)
            new.charmm_instructions['ter'][chain] = termini_tuple

        new.charmm_instructions['patch_ter'] = []
        new.charmm_instructions['patch_ter'] = list(self.charmm_instructions["patch_ter"])

        new.charmm_instructions['cap_termini'] = self.charmm_instructions['cap_termini']
        new.charmm_instructions['ignore_termini'] = self.charmm_instructions['ignore_termini']

        new.charmm_instructions['disu'] = []
        for disu in self.charmm_instructions['disu']:
            residue1, residue2 = disu
            res1 = residue1.copy()
            res2 = residue2.copy()
            disu_tuple = (res1, res2)
            new.charmm_instructions['disu'].append(disu_tuple)

        new.charmm_instructions['missing'] = self.charmm_instructions['missing']

        new.charmm_instructions['patches'] = []
        for patch_name, residue_list in self.charmm_instructions['patches']:
            if type(residue_list) is list:
                patch_residues = []
                for residue_tuple in residue_list:
                    new_tuple = tuple(residue_tuple)
                    patch_residues.append(new_tuple)
            else:
                patch_residues = tuple(residue_list)
            patch = (patch_name, patch_residues)
            new.charmm_instructions['patches'].append(patch)

        new.charmm_instructions['patches_no_autogen'] = []
        for patch_name, residue_list in self.charmm_instructions['patches_no_autogen']:
            # print patch_name
            # print residue_list
            if type(residue_list) is list:
                patch_residues = []
                for residue_tuple in residue_list:
                    new_tuple = tuple(residue_tuple)
                    patch_residues.append(new_tuple)
            else:
                patch_residues = tuple(residue_list)
            patch = (patch_name, patch_residues)
            new.charmm_instructions['patches_no_autogen'].append(patch)

        new.charmm_instructions['minimize_selections'] = list(self.charmm_instructions['minimize_selections'])

        new.charmm_config = {}
        new.charmm_config['output'] = []
        new.charmm_config['output'] = list(self.charmm_config['output'])

        new.charmm_config['missing_atom_ids'] = []
        for missing_id in self.charmm_config['missing_atom_ids']:
            new_missing_id = tuple(missing_id)
            new.charmm_config['missing_atom_ids'].append(new_missing_id)

        new.charmm_config['autogen'] = self.charmm_config['autogen']

        new.charmm_config['tasks'] = []
        for task in self.charmm_config['tasks']:
            new.charmm_config['tasks'].append(task)

        new.charmm_output = self.charmm_output

        new.charmm_normal_termination = self.charmm_normal_termination

        new.title = self.title

        new.charmm_out_prefix = self.charmm_out_prefix

        new.clean_up = self.clean_up

        new.charmm_tmp_files = []
        for file in self.charmm_tmp_files:
            new.charmm_tmp_files.append(file)

        new.titr_residues = {}
        for residue, descr in self.titr_residues.iteritems():
            states = []
            for state in descr:
                newstate = dict(state)
                states.append(newstate)
            new.titr_residues[residue] = states

        new.titr_restypes = {}
        new.titr_restypes = dict(self.titr_restypes)

        new.prot_state = {}
        for residue_tuple, state in self.prot_state.iteritems():
            new_tuple = tuple(residue_tuple)
            new.prot_state[new_tuple] = state

        # new.parse_top()

        # if pdb:
        #     new.add_structure(pdb)

        return new

    def backup_settings(self):

        self.backuped_settings = self.copy(no_structure=True, no_backup=True)

    def restore_settings(self):

        charmm_out_prefix = self.charmm_out_prefix
        backup_settings_copy = self.backuped_settings.copy(no_structure=True, no_backup=True)
        structure = self.structure

        self.__dict__ = self.backuped_settings.__dict__
        self.backuped_settings = backup_settings_copy
        self.structure = structure
        self.charmm_out_prefix = charmm_out_prefix

    def update_coords_with_modelled_structure(self, allow_hydrogen_missmatch=False):

        modelled_structure = self.get_modelled_structure()
        self.transfer_coordinates(modelled_structure, allow_hydrogen_missmatch=allow_hydrogen_missmatch)

    def update_charges_with_modelled_structure(self, allow_hydrogen_missmatch=False):

        modelled_structure = self.get_modelled_structure()
        self.transfer_charges(modelled_structure, allow_hydrogen_missmatch=allow_hydrogen_missmatch)

    def add_patch(self, patch_name, residue_list, no_autogen=False):
        """
        Adds a new patch of user's choice.

        @patch_name: name of the patch
        @residue_list: list of residues tuples, structure is:
        [residue: tuple:(str: <resname>, int: <resid>, str: <segname>)]

        It will be added to the patch list for the specific object of the Charmm_manager Class.
        """
        if type(residue_list) is not list:
            residue_list = [residue_list]
        residue_list_tuples = []
        for residue_descr in residue_list:
            if type(residue_descr) is not tuple:
                resname, resid, segname = re.split('[-_]', residue_descr)
                residue_tuple = (resname, int(resid), segname)
                residue_list_tuples.append(residue_tuple)
            else:
                residue_list_tuples.append(residue_descr)
        if not no_autogen:
            self.charmm_instructions['patches'].append((patch_name, residue_list_tuples))
        else:
            self.charmm_instructions['patches_no_autogen'].append((patch_name, residue_list_tuples))

    def add_charmm_command (self, command, adj_task=None):
        """
        Adds a new command to the charmm script after some task from the list of tasks

        @command: str - charmm command of users choice
        @adj_task:  str - task after which the command will be done e.g. after the adj_task = 'patches'

        """
        if adj_task is not None:
            for i, entry in enumerate(self.charmm_config['tasks']):
                if entry == adj_task:
                    self.charmm_config['tasks'].insert(i+1, command)
                    break
        else:
            raise AssertionError('Not a complete input: State a task or a line!')

    def set_titr_residues(self, titr_residues):
        """
        Stores the definition of titratable residues.

        @param titr_residues: For definition of titr_residues, see function 'get_std_titr_residues'.
        @return: None
        """

        titr_restypes = {}
        for resname, titr_residue in titr_residues.iteritems():
            titr_restypes[resname] = resname
            for state_descr in titr_residue:
                if state_descr['rename'] is not None:
                    alternative_name = state_descr['rename']
                    titr_restypes[alternative_name] = resname

        self.titr_residues = titr_residues

        # This list is needed to find the corresponding entry in self.titr_residues, if the residue has been renamed.
        self.titr_restypes = titr_restypes

    def get_titr_residues(self):
        """
        Creates and returns of copy of the titratable residues dictionary.
        {}
        Returns a copy of self.titr_residues
        """
        titr_residues = {}
        for resname, states in self.titr_residues.iteritems():
            titr_residues[resname] = []
            for state in states:
                copied_state = dict(state)
                titr_residues[resname].append(copied_state)
        return titr_residues


    def add_decision(self, decision, option):
        """
        Adds a decision to handle a specific problem regarding the structure.
        Decisions and problem can occur from different sources.
        @decision:str
        @option:str

        e.g. decision: 'disu_bridges', option: 'closed' -> In this structure salt bridges will be closed!
        """
        self.__check_status('add_decision')


#        reg = re.compile(r'^(.+?)__.+')
#        reg_m = reg.match(problem)
#        problem_type = reg_m.groups()[0]
        problem_type = decision.split('__')[0]
        #                                  Was it used?
        self.decisions[decision] = (option, False, problem_type)

    def __get_decision(self, decision):
        if self.decisions.has_key(decision):
            (option, used, problem_type) = self.decisions[decision]
            self.decisions[decision] = (option, True, problem_type)
            return option
        else:
            return None

    def is_allowed(self, fkt_name):
        """
        Check whether the function call of 'fkt_name' is allowed or not. Nothing will be added to the log.
        """
        (permission, message) = self.__check_status(fkt_name, just_asking=True)
        return permission

    def __log(self, message):
        self.log.append(message)

    def __warn(self, message):
        self.log("WARNING: " + message)


    def __check_status(self, fkt_name, just_asking=False):
        """
        Check if the call of the function 'fkt_name' is allowed.

        Keep track of activities in this class. Every function in Charmm_manager that changes the internal
        stored data should call this function with its own name as 'fkt_name' before it starts changing anything. With
        just_asking=False it is assumed that this function has been called by 'fkt_name' and the execution of
        'fkt_name' will be logged if permission is granted. If permission is denied the calling class should stop
        immediately.

        Returns the tuple (permission, message).
        permission: True if the function call is allowed.
        message:    Is '' if permission is granted and contains a string with the error message if not.
        """

        message = None
        permission = False

        if fkt_name == 'add_structure':
            permission = True
            message = ''
            self.structure_read = True
            if self.structures_checked:
                if self.structures_prepared:
                    self.__warn("__check_status: New structure added although the structures were prepared already. " +\
                                "Structures must be checked and prepared again.")
                else:
                    self.__warn("__check_status: New structure added although the structures were checked already. " +\
                                "Structures must be checked again.")
            self.structures_checked  = False
            self.structures_prepared = False
        elif fkt_name == 'parse_top':
            permission = True
            message = ''
        elif fkt_name == 'add_decision':
            permission = True
            message = ''
            if self.structures_checked:
                permission = False
                message = "__check_status: The function 'check_structures' has been called already. Calling the function "\
                          "'add_decision' is not allowed any more."
                if not just_asking:
                    raise self.GeneralError(message)
        elif fkt_name == 'check_structures':
            permission = True
            message = ''
            if not self.structure_read:
                permission = False
                message = "__check_status: A structure must be added with 'add_structure' before 'check_structures' " +\
                          "can be called."
                if not just_asking:
                    raise self.GeneralError(message)
            self.structures_checked = True

        elif fkt_name == 'run_charmm':
            permission = True
            message = ''
            if not self.top:
                permission = False
                message = "__check_status: Topology file is missing! Topology file must be provided before 'run_charmm' "\
                          "can be called."
                if not just_asking:
                    raise self.GeneralError(message)
            if not self.par:
                permission = False
                message = "__check_status: Parameter file is missing! Parameter file must be provided before 'run_charmm' "\
                          "can be called."
                if not just_asking:
                    raise self.GeneralError(message)

            if not self.top_content:
                permission = False
                message = "__check_status: No topology content! Parse topology file before calling 'run_charmm'"
                if not just_asking:
                    raise self.GeneralError(message)

        if message is None:
            message = "__check_status: Function '%s' is unknown." % fkt_name
            if not just_asking:
                raise self.GeneralError(message)

        return (permission, message)

    def parse_top(self):
        """
        Parses the topology files.
        So far only names and charges are read.
        """

        (permission, message) = self.__check_status('parse_top')
        if not permission:
            print(message)
            return -1

        for filename in self.top:
            with open(filename, 'r') as f:
                last_residue = None
                last_entry_is_patch = False
                for line in f:
                    # Consider that the line is commented out or empty.
                    if not line.strip() or line.strip()[0] == '!':
                        continue
                    # Consider that a part of the line is commented out.
                    line = line.split('!')[0]
                    entries = line.split()
                    if len(entries) == 0:
                        continue
                    if entries[0] in ['RESI', 'PRES', 'END']:
                        if last_residue is not None:
                            # End of residue entry found. Store previous residue.
                            name = last_residue['name']
                            if last_entry_is_patch:
                                self.top_content['patches'][name] = last_residue
                            else:
                                self.top_content['residues'][name] = last_residue
                            last_residue = None

                    if entries[0] in ['RESI', 'PRES']:
                        # New residue found.
                        if entries[0] == 'PRES':
                            last_entry_is_patch = True
                        else:
                            last_entry_is_patch = False
                        last_residue = {}
                        last_residue['atoms'] = []
                        last_residue['name'] = entries[1]
                        if len(entries) > 2:
                            last_residue['charge'] = int(float(entries[2]))
                        else:
                            last_residue['charge'] = None
                    elif entries[0] == 'ATOM':
                        # New ATOM entry found.
                        if last_residue is not None:
                            # The new ATOM entry belongs to a RESI entry. Adding it to the last found residue.
                            new_atom = {}
                            new_atom['name']   = entries[1]
                            new_atom['type']   = entries[2]
                            new_atom['charge'] = entries[3]
                            last_residue['atoms'].append(new_atom)


    def add_structure(self, structure):
        """
        Set the structure.

        Parameter:
        structure: string with path to a pdb/crd file or a Simple_struct_parser object.
        """

        name = ''

        if isinstance(structure, str):
            pdb = Simple_struct_parser()
            suffix = structure.split('.')[-1]
            if suffix.upper().find('PDB') > -1:
                pdb.read_pdb(structure)
            elif suffix.upper().find('CRD') > -1:
                pdb.read_crd(structure)
            else:
                error = "Could not guess file type of '%s' or the file type is not supported." % structure
                raise AssertionError(error)
            name = structure
            # Generate prefix
            end_of_path = name.rfind('/') + 1
            start_of_ending = name.rfind('.')
            title = name[end_of_path:start_of_ending]
        elif isinstance(structure, Simple_struct_parser):
            pdb = structure
            name = pdb.file_log[0] + ' (Simple_struct_parser object) '

            # Generate prefix
            end_of_path = name.rfind('/') + 1
            start_of_ending = pdb.file_log[0].rfind('.')
            title = pdb.file_log[0][end_of_path:start_of_ending]

        else:
            print("ERROR in Charmm_manager::add_structure. Parameter is not of type string or Simple_struct_parser: ")
            print("-> " + str(type(structure)))
            print()
            return -1

        (permission, message) = self.__check_status('add_structure')
        if not permission:
            print(message)
            return -1

        pdb.remove_altloc()

        self.structure = pdb

        self.__log("Structure '%s' added." % name)

        self.title = title

        # Todo: This line make the __handle_icodes function obsolete.
        self.structure.resolve_icodes()

        return 1

    def transfer_coordinates(self, new_structure, allow_hydrogen_missmatch=False):
        """Transfers the coordinates from new_structure to the stored structure. For each atom in structure, there
        must be an corresponding atom in structure.

        new_structure: Simple_structure_parser object
        """

        structure = self.structure
        for new_atom in new_structure.atoms:
            resid = new_atom['resid']
            segname = new_atom['segname']
            name = new_atom['name']
            if (structure.struct.has_key(segname)) and (structure.struct[segname].has_key(resid)) and \
                    (structure.struct[segname][resid].has_key(name)):
                current_atom = structure.struct[segname][resid][name]
                current_atom['coord'] = np.array(new_atom['coord'])
            else:
                name = new_atom['name']
                element = name[0]
                if not (element == 'H' and allow_hydrogen_missmatch):
                    error = "Atom %s in resid %i in segment %s does not exist in the current structure." \
                            % (name, resid, segname)
                    raise AssertionError(error)

    def transfer_charges(self, new_structure, allow_hydrogen_missmatch=False):
        """Transfers the charges from new_structure to the stored structure. For each atom in structure, there
        must be an corresponding atom in structure.

        new_structure: Simple_structure_parser object
        """

        structure = self.structure
        for new_atom in new_structure.atoms:
            resid = new_atom['resid']
            segname = new_atom['segname']
            name = new_atom['name']
            if (structure.struct.has_key(segname)) and (structure.struct[segname].has_key(resid)) and \
                    (structure.struct[segname][resid].has_key(name)):
                current_atom = structure.struct[segname][resid][name]
                current_atom['charge'] = new_atom['charge']
            else:
                name = new_atom['name']
                element = name[0]
                if not (element == 'H' and allow_hydrogen_missmatch):
                    error = "Atom %s in resid %i in segment %s does not exist in the current structure." \
                            % (name, resid, segname)
                    raise AssertionError(error)



    # Add a sequence that will follow first_res. HIS will be renamed to HSD.
    # Parameters:
    # first_res : <simple_atom_access_residue object>
    # sequence  : Sequence to add in space separated three letter code.
    def add_sequence(self, first_res, sequence, is_last_residue=False):
        pdb = self.structure

        temlate_atom = first_res.get_first_atom().copy()
        temlate_atom['icode'] = ' '
        temlate_atom['coord'] = np.array( [9999, 9999, 9999] )
        temlate_atom['mass'] = None
        temlate_atom['vdw'] = None
        temlate_atom['occupancy'] = 0.0
        temlate_atom['tempFactor'] = 0.0
        temlate_atom['index'] = 0

        if not is_last_residue:
            resid = first_res.resid + 1
        else:
            resid = first_res.resid - len(sequence.split(' '))

        if not is_last_residue:
            atom_index = 0
            for i, atm in enumerate(pdb.atoms):
                for atm_ref in first_res.iter_atoms():
                    if atm is atm_ref:
                        atom_index = i
            atom_index += 1
        else:
            atom_index = 0
            stop = False
            for i, atm in enumerate(pdb.atoms):
                for atm_ref in first_res.iter_atoms():
                    if atm is atm_ref:
                        atom_index = i
                        stop = True
                        break
                if stop:
                    break



        for resname in sequence.split(' '):
            if resname == 'HIS':
                resname = 'HSD'
            temlate_atom['resname'] = resname
            temlate_atom['resid']   = resid
            # From topology file:
            atom_list = self.top_content['residues'][resname]['atoms']
            for atm in atom_list:
                temlate_atom['name']   = atm['name']
                temlate_atom['type']   = atm['type']
                temlate_atom['charge'] = atm['charge']

                temlate_atom['element'] = atm['name'][0]
                pdb.atoms.insert(atom_index, temlate_atom.copy())
                atom_index += 1

            resid += 1

        pdb.create_struct()

        pdb.reassign_indices()

        # for i, atm in enumerate(pdb.atoms):
        #     atm.__list_position_in_host = i


        return 1

    def check_structures(self, quiet=False):
        """
        generelles Verhalten: default Aktion wird ausgeführt, Alternativen werden kommentiert.

        -> unknwon residues
          -> delete (default)
          -> keep

        -> fix segment names
          -> hetatm(s) in atom chain:
            -> keep if connected to atoms (recommended)
            -> keep
            -> move to new chain
          -> hetatm only chain:
            -> move all monomers in one chain / move all polymers in seperate chains (recommended)
            -> move all in one chain
            -> do nothing

        -> non-std AA / hetatm grouped by chain (several choices possible if not stated otherwise)
          -> keep all (default / no other options are allowed)
          -> delete all (no other options are allowed)
          -> delete some
          -> delete chain X
          -> delete residue Y
          -> delete type
          -> delete water


        -> gaps
        -> missing residues (if header available)

        -> Disu patching

        -> charge pattern

        -> print information
        """

        (permission, message) = self.__check_status('check_structures')
        if not permission:
            print(message)
            return -1

        self.__handle_charmm_renaming()

        self.__handle_chain_assignment()

        self.__handle_icodes()

        #self.__check_structure()

        self.__handle_protonation()

        self.__handle_unknown()

        self.__handle_gaps()

        self.__handle_termini()

        # __handle_termini must be called first.
        self.__handle_charmm_atom_names()

        self.__handle_disu_bridges()

        self.__find_missing_atoms()


        if not quiet:
            if len(self.decided) > 0:
                print("\n### The following changes have been applied to the structure. ###\n")
                for decided in self.decided:
                    print(decided + '\n')

            if len(self.undecided) > 0:
                print("\n### There are undecided decisions for whose no default behaviour exists. ###\n")
                for undecided in self.undecided:
                    print(undecided + '\n')

        return True

    def __handle_charmm_renaming(self):
        # renaming scheme:
        rename = {
            # 'CA'  : 'CAL',
            # 'THP' : 'THY',
            'HOH' : 'TIP3',
            'MSE' : 'MET',
            # 'CL'  : 'CLA',
            # 'ZN'  : 'ZN2'
        }

        resnames_to_rename = {}
        pdb = self.structure
        for segment in pdb.struct.iter_segments():
            for residue in segment.iter_residues():
                resname = residue.resname
                if rename.has_key(resname):
                    replacement   = rename[resname]
                    decision_name = "rename__%s_%s" % (resname, replacement)

                    decision = self.__get_decision(decision_name)
                    # Set default
                    if decision is None:
                        decision = 'rename'
                    if decision == 'rename':
                        for atom in residue.iter_atoms():
                            atom['resname'] = replacement


                    # Has this resname been handled already?
                    if not resnames_to_rename.has_key(resname):
                        resnames_to_rename[resname] = replacement

                        text =  "# Residue of type " + str(resname) + " must be renamed in order to be recognised by Charmm.\n"
                        text += "  Decision: " + decision_name + "\n"
                        text += "  Options:  rename - Rename to " + replacement + " (default)\n"
                        text += "            keep   - Do nothing.\n"
                        text += "  Chosen:   " + decision

                        self.decided.append(text)

        pdb.create_struct()

        return True

    def __handle_chain_assignment(self):
        pdb = self.structure
        for segment in pdb.struct.iter_segments():

            # Check that all waters are in separate chains.
            waters = []
            is_pure_water_chain = True
            for residue in segment.iter_residues():
                if residue.resname == 'TIP3':
                    waters.append(residue)
                else:
                    is_pure_water_chain = False
            # If the segment contains water and non water residues, move the waters to a separate chain.
            if waters and not is_pure_water_chain:
                # Determine new segname for the waters.
                for i in ['',1,2,3,4,5,6,7,8,9]:
                    if not pdb.struct.has_key('HOH' + str(i)):
                        new_segname = 'HOH' + str(i)
                        break
                else:
                    for i in range(10,100):
                        if not pdb.struct.has_key('W' + str(i)):
                            new_segname = 'W' + str(i)
                            break
                for water in waters:
                    for atm in water.iter_atoms():
                        atm['segname'] = new_segname
                        atm['chain'] = ' '

                pdb.create_struct()





    def __handle_protonation(self):
        histidines = []
        acids      = []
        bases      = []
        defaults = 0

        pdb = self.structure
        for segment in pdb.struct.iter_segments():
            for residue in segment.iter_residues():
                resname = residue.resname

                #if resname in ('ASP GLU'):
                #if resname in ('ARG LYS'):
                if resname == 'HIS':
                    resid   = residue.resid
                    segname = residue.segname

                    #decision_name = "prot__%s-%i_%s" % (resname, residue.resid, residue.segname)
                    decision_name = "prot__%s" % (residue)
                    decision = self.__get_decision(decision_name)
                    if decision is None:
                        # Try chainid
                        decision_name = "prot__%s-%i_%s" % (resname, residue.resid, residue.segname[0])
                        decision = self.__get_decision(decision_name)

                    # Set default
                    if decision is None:
                        decision = 'HSE'
                        defaults += 1

                    histidines.append( (residue, decision) )

                    # for atom in residue.iter_atoms():
                    #     atom['resname'] = decision

                    residue_tuple = (resname, resid, segname)

                    self.set_prot_residue(residue_tuple, rename=decision)



                #if self.top_content['known_residues'].has_key(resname):
                #if abs(self.top_content['known_residues'][resname]) > 0.01:

        if len(histidines) > 0:
            if defaults > 0:
                insert = "The default setting will be used for %i of them." % defaults
            else:
                insert = "All residues have been set by the user."
            text =  "# There are %s HIS residues in the structure. %s\n" % (len(histidines), insert)
            text += "  Decision: " + "prot__<RESNAME>-<RESID>_<SEGNAME/CHAIN>" + "\n"
            text += "  Options:  HSE - Set to HSE. (default)\n"
            text += "            HSD - Set to HSD.\n"
            text += "            HSP - Set to HSP.\n"
            text += "  Chosen:  "
            c = 0
            for residue, decision in histidines:
                if c > 4:
                    text += "\n           "
                    c = 0
                text += (" " + str(residue) + ":" + decision).ljust(16)
                c += 1

            self.decided.append(text)


        # pdb.create_struct()

        return True

    def __handle_unknown(self):
        print_message = False
        unknown_resnames = defaultdict(int)

        pdb = self.structure
        del_list = []
        for segment in pdb.struct.iter_segments():
            for residue in segment.iter_residues():
                resname = residue.resname

                assert resname != 'HIS', "Histidine with undefined protonation state found! (%s-%i_%s)" \
                                         % (resname, residue.resid, residue.segname)

                if not self.top_content['residues'].has_key(resname):
                    print_message = True
                    unknown_resnames[resname] += 1
                    decision_name = "delete_unknown"
                    decision = self.__get_decision(decision_name)
                    # Set default
                    if decision is None:
                        decision = 'yes'
                    if decision == 'yes':
                        for atm in residue.iter_atoms():
                            del_list.append(atm)
        if del_list:
            pdb.del_atom(del_list)

        # pdb.create_struct()

        if print_message:
            resname_text = ""
            for resname, count in unknown_resnames.iteritems():
                resname_text += " %ix%s" % (count, resname)
            resname_text = resname_text.strip()

            text =  "# Some residues have been found that are not known to Charmm. (%s)\n" % resname_text
            text += "  Decision: " + decision_name + "\n"
            text += "  Options:  yes - Delete the residues. (default)\n"
            text += "            no  - Keep them.\n"
            text += "  Chosen:   " + decision

            self.decided.append(text)


        return True

    def __handle_gaps(self):
        # Set the default behaviour.
        decision_name = "gap__model_gaps"
        decision = self.__get_decision(decision_name)
        if decision is None:
            default_decision = 'model/keep'
            #default_decision = 'keep'
        else:
            default_decision = decision



        pdb = self.structure
        # Look at each chain separately.
        for segment in pdb.struct.iter_segments():
            segname = segment.get_segname()

            # Look for gaps in the residue numbering.
            last_resid   = None
            last_resname = None
            for residue in segment.iter_residues():
                resid   = residue.resid
                resname = residue.resname
                if residue.hetatm or resname in ['TIP3', 'HOH', 'POPC', 'POPE']:
                    continue
                if (last_resid is None) or (resid == last_resid + 1):
                    last_resid   = resid
                    last_resname = resname
                    last_residue = residue
                else:
                    true_gap = True
                    gap_start = last_residue
                    gap_end   = residue
                    missing_residues = last_resid - resid

                    decision_name = "gap__%s-%i_%s__%s-%i_%s" \
                            % (last_resname, last_resid, segname,  resname, resid, segname)
                    decision = self.__get_decision(decision_name)

                    # Check if the gap is just a gap in the numbering and not in the structure.
                    for atom_start in gap_start.iter_atoms():
                        if atom_start['name'][0] == 'H':
                            continue
                        coord = atom_start['coord']
                        for atom_end in gap_end.iter_atoms():
                            if atom_start['name'][0] == 'H':
                                continue
                            coord_comp = atom_end['coord']
                            dist = coord_comp - coord
                            if 0.5 < np.sqrt(np.dot( dist, dist )) < 1.7:
                                true_gap = False


                    #  Try to extract the missing sequence from the header.
                    missing_sequence = []
                    chainid = segment.chainid
                    for seq in pdb.header_info['missing']:
                        if (seq[0][1]  == gap_start.chainid) and (seq[0][2]  == last_resid + 1) and\
                           (seq[-1][1] == gap_end.chainid)   and (seq[-1][2] == resid - 1):
                            missing_sequence = seq
                            missing_sequence_str = ""
                            for rn, c, rid in seq:
                                missing_sequence_str += ' ' + rn
                            missing_sequence_str = missing_sequence_str.strip()
                            break


                    # Set default behaviour.
                    if decision is None:
                        if true_gap:
                            decision = default_decision
                            # decision = 'model/keep'
                            #if is_std_aa(resname) and is_std_aa(last_resname):
                            #    if missing_residues < 10:
                            #        decision = '...'

                        else:
                            # Charmm is taking care of that itself.
#                            decision = 'connect'
                            pass

                    # Todo: more feedback, on what will be actually done. Maybe give a selection commant to inspect modelled residues?
                    guessed_sequence_accepted = False
                    if decision == 'connect':
                        # Todo: Warning if the gap is a true gap.
                        self.charmm_instructions['gaps'].append( ('connect', gap_start, gap_end) )
                    elif decision == 'keep':
                        pass
                    elif decision == 'model/keep':
                        if 0 < len(missing_sequence) < 15:
                            self.charmm_instructions['gaps'].append( ('fill', gap_start, gap_end, missing_sequence_str) )
                            self.add_sequence(gap_start, missing_sequence_str)
                            guessed_sequence_accepted = True
                    elif decision is not None:
                        # Todo: Warning if gap size is not equal to sequence length.
                        sequence = decision
                        self.charmm_instructions['gaps'].append( ('fill', gap_start, gap_end, sequence) )
                        self.add_sequence(gap_start, sequence)

                    if true_gap:
                        text =  "# A gap has been found in the resid numbering AND in the structure.\n"
                        if guessed_sequence_accepted:
                            text += "  Suggestion for the missing sequence from the header:\n"
                            text += " "
                            c = 1
                            for rn in missing_sequence_str.split(' '):
                                if c > 15:
                                    text += "\n "
                                    c = 0
                                text += " " + rn
                                c += 1
                            text += "\n"
                        text += "  Decision: " + decision_name + "\n"
                        text += "  Options:  keep       - Do nothing, the gap remains in the structure.\n"
                        text += "            model/keep - Try to extract the missing sequence from the PDB header.\n"
                        text += "                         If possible and if the length of the missing sequence is\n"
                        text += "                         smaller than 10, the gap will be modelled. Otherwise the gap \n"
                        text += "                         will be kept.(default)\n"
                        text += "            <sequence> - Model the missing residues with the given sequence. \n"
                        text += "                         The sequence must be space separated three letter codes. \n"
                        text += "  Chosen:   " + decision
                        self.decided.append(text)
#                    else:
#                        text =  "# A gap has been found in the resid numbering, BUT NOT in the structure.\n"
#                        text += "  Decision: " + decision_name + "\n"
#                        text += "  Options:  connect - Link the two residues by applying a patch in Charmm. (default)\n"
#                        text += "            keep    - Do nothing, the gap remains in the structure.\n"
#                        text += "  Chosen:   " + decision
#                        self.decided.append(text)



                    last_resid   = resid
                    last_resname = resname

        return True

    def __handle_termini(self):
        decision_name = 'cap_termini'
        text = "# Capping of termini will be decided here. Should they be caped?\n"
        text += "  Decision: " + decision_name + "\n"
        text += "  Options:  cap - Apply patches in Charmm to cap termini if not true termini. (default)\n"
        text += "            dont_cap  - Do nothing, apply regular termini patches!"
        self.decided.append(text)
        return self.find_termini()

    def find_termini(self):
        pdb = self.structure
        # Look at each chain separately.
        # decision_name = 'cap_termini'
        # text = "# Capping of termini will be decided here. Should they be caped?\n"
        # text += "  Decision: " + decision_name + "\n"
        # text += "  Options:  cap - Apply patches in Charmm to cap termini if necessary. (default)\n"
        # text += "            dont_cap  - Do nothing, apply regular termini patches! \n"

        for segment in pdb.struct.iter_segments():
            segname = segment.get_segname()

            # Find first and last standard amino acid in this segment.
            nter = None
            cter = None

            f_ter = None
            t_ter = None
            for residue in segment.iter_residues():
                if is_std_aa(residue.resname):

                    # has_std_aa = True
                    if nter is None:
                        nter = residue
                    else:
                        cter = residue

                elif is_na_base(residue.resname):
                    if f_ter is None:
                        f_ter = residue
                    else:
                        t_ter = residue

            # Check if the residues are real N- and C-Termini. For N Termini resid=1 is used as check.
            # if (nter is not None) and ((nter.resid != 1) or ('CAY' in nter)):
            #     nter = None
            # if (cter is not None) and not ((('OXT' in cter) or ('OT1' in cter) and not ('HT1' in cter))):
            #     cter = None


            if (f_ter is not None) and (t_ter is not None):
                self.charmm_instructions['ter'][segname] = (f_ter, t_ter)
            else:
                self.charmm_instructions['ter'][segname] = (nter, cter)


            decision_name = 'cap_termini'
            decision = self.__get_decision(decision_name)
            # Set default
            if decision is None:
                decision = 'cap'

            # text = "Termini capping - Chosen decision:   " + decision + " for segname %s \n" %segname
            text = "  Chosen:   " + decision + " for segname %s" %segname

            self.decided.append(text)

            if decision == 'cap':
                self.charmm_instructions['cap_termini'] = True
            if decision == 'dont_cap':
                self.charmm_instructions['cap_termini'] = False
            if decision == 'ignore_termini':
                print 'WARNING: This is inappropriate modelling and it should be used only for testing!!!'
                self.charmm_instructions['ignore_termini'] = True

        # To set termini patching to None and avoid atom renaming:
        # self.charmm_instructions['ter'][segname] = (None, None)

        # To set termini patching to None:
        # self.charmm_instructions['patch_ter'] = {'A' : None}

        return True

    def __handle_charmm_atom_names(self):
        pdb = self.structure
        # Look at each chain separately.
        for segment in pdb.struct.iter_segments():
            segname = segment.get_segname()

            # Change C-terminus atom names, if it is a true terminus. True termini have an OXT atom.
            is_ter = False
            cter = self.charmm_instructions['ter'][segname][1]
            if cter is not None:
                for atom in cter.iter_atoms():
                    if atom['name'] == 'OXT':
                        atom['name'] = 'OT1'
                    if atom['name'] == 'O':
                        atom['name'] = 'OT2'
                # for atom in cter.iter_atoms():
                #     if atom['name'] == 'OXT':
                #         atom['name'] = 'OT2'
                #         is_ter = True
                # if is_ter:
                #     for atom in cter.iter_atoms():
                #         if atom['name'] == 'O':
                #             atom['name'] = 'OT1'

            for residue in segment.iter_residues():
                if residue.resname == 'ILE':
                    for atom in residue.iter_atoms():
                        if atom['name'] == 'CD1':
                            atom['name'] = 'CD'
                if residue.resname == 'TIP3':
                    for atom in residue.iter_atoms():
                        if atom['name'] == 'O':
                            atom['name'] = 'OH2'
                if residue.resname == 'CAL':
                    for atom in residue.iter_atoms():
                        if atom['name'] == 'CA':
                            atom['name'] = 'CAL'
                # For MSE to MET conversion.
                if residue.resname == 'MET':
                    for atom in residue.iter_atoms():
                        if atom['name'] == 'SE':
                            atom['name'] = 'SD'
                if residue.resname == 'CLA':
                    for atom in residue.iter_atoms():
                        if atom['name'] == 'CL':
                            atom['name'] = 'CLA'

        self.structure.create_struct()

        return True

    def __handle_disu_bridges(self):
        # Find all CYS residues
        cys_list = []
        pdb = self.structure
        for segment in pdb.struct.iter_segments():
            for residue in segment.iter_residues():
                if residue.resname == 'CYS':
                    cys_list.append( residue )


        # Find disulphide bridges. Will fail if SG atom is not resolved.
        disu_bridges = []
        for i, cys_residue in enumerate(cys_list):
            if not cys_residue.has_key('SG'):
                continue
            sg_atom_coord = cys_residue['SG']['coord']
            for cys_residue_comp in cys_list[i+1:]:
                if not cys_residue_comp.has_key('SG'):
                    continue
                sg_atom_comp_coord = cys_residue_comp['SG']['coord']
                delta = sg_atom_coord - sg_atom_comp_coord
                if 1.8 < np.sqrt(np.dot( delta, delta )) < 2.2:
                    disu_bridges.append( (cys_residue, cys_residue_comp) )

        if disu_bridges:
            decision_name = 'disu_bridges'
            decision = self.__get_decision(decision_name)
            # Set default
            if decision is None:
                decision = 'closed'

            text =  "# There are %i disulphide bridges in the structure. Should they be closed?\n" % len(disu_bridges)
            text += "  Decision: " + decision_name + "\n"
            text += "  Options:  closed - Apply patches in Charmm to connect the CYS pairs. (default)\n"
            text += "            open   - Do nothing, CYS pairs will NOT be bound! \n"
            text += "  Chosen:   " + decision

            self.decided.append(text)

            if decision == 'closed':
                self.charmm_instructions['disu'] = disu_bridges


        return True


    def __handle_icodes(self):
        # Todo: store the information about renamed residues.

        # # Find residues with icodes.
        # pdb = self.structure.struct
        # icode_residues = []
        # icode_residues_sequ = []
        # last_residue = None
        # for residue in pdb.iter_residues():
        #     if residue.icode != ' ':
        #         if last_residue is None:
        #             icode_residues_sequ.append(residue)
        #         else:
        #             if last_residue.resid == residue.resid:
        #                 icode_residues_sequ.append(residue)
        #             else:
        #                 icode_residues.append(icode_residues_sequ)
        #                 icode_residues_sequ = [residue]
        #         last_residue = residue
        #
        # if len(icode_residues_sequ) > 0:
        #     icode_residues.append(icode_residues_sequ)

        # if len(icode_residues) > 0:
        #
        #     # WARNING: Can change the resid numbering!
        #
        #     decision_name = 'icode_renaming'
        #     decision = self.__get_decision(decision_name)
        #     # Set default
        #     if decision is None:
        #         decision = 'shift'
        #
        #     text  = "# There are %i inserted sequences in the structure (residues with icode entries). The residues in\n" % len(icode_residues)
        #     text += "  the corresponding chains have to be renumbered. How should this be done?\n"
        #     for icode_residues_sequ in icode_residues:
        #         first = icode_residues_sequ[0]
        #         last  = icode_residues_sequ[-1]
        #         text += "    %s_%i[%s]_%s -> %s_%i[%s]_%s\n" % (first.resname, first.resid, first.icode, first.segname,\
        #                                                 last.resname, last.resid, last.icode, last.segname, )
        #     text += "  Decision: " + decision_name + "\n"
        #     text += "  Options:  shift - All residues that follow a residue with icode are increased by one. (default)\n"
        #     text += "  Chosen:   " + decision
        #
        #     self.decided.append(text)
        #
        #     if decision == 'shift':
        #         self.structure.resolve_icodes()

        if self.structure.icodes > 0:
            # Todo: It would be more useful to really count the residues. structure.icode containes the number of atoms
            icode_atoms = self.structure.icodes
            # decision_name = 'icode_renaming'
            # decision = self.__get_decision(decision_name)
            # # Set default
            # if decision is None:
            #     decision = 'shift'

            text  = "# There are %i inserted atoms in the structure (atoms with icode entries). The corresponding " \
                    "chains had to be renumbered.\n" % icode_atoms
            # text += "  Decision: " + decision_name + "\n"
            # text += "  Options:  shift - All residues that follow a residue with icode are increased by one. (default)\n"
            # text += "  Chosen:   " + decision

            self.decided.append(text)

            # Has been done after reading the structure.
            # if decision == 'shift':
            #     self.structure.resolve_icodes()



    def __check_structure(self):
        # Find overlapping atoms.
        pdb = self.structure
        atom_coord_list = [atom['coord'] for atom in pdb.atoms]
        atom_coord_list = np.array(atom_coord_list, dtype=np.double)
        overlapping_atoms = find_collisions(atom_coord_list, 0.5)

        # Get more context for colliding atoms.
        protein_collisions = []
        ligand_collisions = []
        protein_ligand_collisions = []
        for index1, index2 in overlapping_atoms:
            atom1 = pdb.get_atom_ref(index1)
            atom2 = pdb.get_atom_ref(index2)

            # Ignore atoms with undefined coordinates.
            if atom1['coord'][0] > 998:
                continue

            if atom1['hetatm'] and atom2['hetatm']:
                ligand_collisions.append( (atom1, atom2) )
            elif not atom1['hetatm'] and not atom2['hetatm']:
                protein_collisions.append( (atom1, atom2) )
            elif (not atom1['hetatm'] and atom2['hetatm']):
                # Protein atom is stored at first position.
                protein_ligand_collisions.append( (atom1, atom2) )
            elif (atom1['hetatm'] and not atom2['hetatm']):
                # Protein atom is stored at first position.
                protein_ligand_collisions.append( (atom2, atom1) )

        undecided_decisions = False

        ### 1) Protein collisions ###
        if len(protein_collisions) > 0:
            decision_name = 'protein_collisions'
            decision = self.__get_decision(decision_name)
            # Set default
            if decision is None:
                decision = '--None--'

            text =  u"# There are %i collisions in the structure between protein atoms (distance < 0.5\u212B):\n" % len(protein_collisions)
            for first_atom, second_atom in protein_collisions:
                text += "  %s-%i_%s:%s   <-->   %s-%i_%s:%s\n" % \
                        (first_atom['resname'],first_atom['resid'],first_atom['segname'],first_atom['name'],\
                         second_atom['resname'],second_atom['resid'],second_atom['segname'],second_atom['name'])
            text += "  Decision: " + decision_name + "\n"
            text += "  Options:  delete - Delete the atoms. The first on of each pair is conserved.\n"
            text += "            keep   - Do nothing. \n"
            text += "  Chosen:   " + decision

            if decision == '--None--':
                self.undecided.append(text)
                undecided_decisions = True
            else:
                self.decided.append(text)
                # Apply decision
                if decision == 'delete':
                    del_list = []
                    for first_atom, second_atom in protein_collisions:
                        del_list.append(second_atom)
                    pdb.del_atom(del_list)

        ### 2) Protein-Ligand collisions ###
        if len(protein_ligand_collisions) > 0:
            decision_name = 'protein_ligand_collisions'
            decision = self.__get_decision(decision_name)
            # Set default
            if decision is None:
                decision = '--None--'

            text =  u"# There are %i collisions in the structure between protein and ligand atoms (distance < 0.5\u212B):\n" % len(protein_ligand_collisions)
            for first_atom, second_atom in protein_ligand_collisions:
                text += "  %s-%i_%s:%s <---> %s-%i_%s:%s\n" % \
                        (first_atom['resname'],first_atom['resid'],first_atom['segname'],first_atom['name'],\
                         second_atom['resname'],second_atom['resid'],second_atom['segname'],second_atom['name'])
            text += "  Decision: " + decision_name + "\n"
            text += "  Options:  delete - Delete the atoms of the ligand.\n"
            text += "            keep   - Do nothing. \n"
            text += "  Chosen:   " + decision

            if decision == '--None--':
                self.undecided.append(text)
                undecided_decisions = True
            else:
                self.decided.append(text)
                # Apply decision
                if decision == 'delete':
                    del_list = []
                    for first_atom, second_atom in protein_ligand_collisions:
                        del_list.append(second_atom)
                    pdb.del_atom(del_list)

        ### 3) Ligand collisions ###
        if len(ligand_collisions) > 0:
            decision_name = 'ligand_collisions'
            decision = self.__get_decision(decision_name)
            # Set default
            if decision is None:
                decision = '--None--'

            text = u"# There are %i collisions in the structure between ligand atoms (distance < 0.5\u212B):\n" % len(ligand_collisions)
            for first_atom, second_atom in ligand_collisions:
                text += "  %s-%i_%s:%s <---> %s-%i_%s:%s\n" % \
                        (first_atom['resname'],first_atom['resid'],first_atom['segname'],first_atom['name'],\
                         second_atom['resname'],second_atom['resid'],second_atom['segname'],second_atom['name'])
            text += "  Decision: " + decision_name + "\n"
            text += "  Options:  delete - Delete the atoms. The first on of each pair is conserved.\n"
            text += "            keep   - Do nothing. \n"
            text += "  Chosen:   " + decision

            if decision == '--None--':
                self.undecided.append(text)
                undecided_decisions = True
            else:
                self.decided.append(text)
                # Apply decision
                if decision == 'delete':
                    del_list = []
                    for first_atom, second_atom in ligand_collisions:
                        del_list.append(second_atom)
                    pdb.del_atom(del_list)
                # Todo: Check if the whole molecule collides -> Then it is save to delete one molecule.


        if undecided_decisions:
            return False
        else:
            return True


    def __find_missing_atoms(self):
        # Check that no atom is missing.
        pdb = self.structure

        missing_hydrogens_water = 0
        missing_hydrogens_protein = 0
        missing_hydrogens_ligands = 0
        missing_atoms_protein = 0
        missing_atoms_ligands = 0
        missing_atoms_protein_list = {}
        missing_atoms_ligands_list  = {}
        # A list of the ID of all atoms that are missing.
        missing_atoms_list = []

        # Generate a list of CYS residues involved in disulfid bridges. This is
        # needed later to check for missing hydrogens.
        disu_cys_list = []
        for (cys1, cys2) in self.charmm_instructions['disu']:
            disu_cys_list.append(cys1)
            disu_cys_list.append(cys2)

        for seg in pdb.struct.iter_segments():
            for res in seg.iter_residues():
                resname = res.resname
                atom_list = self.top_content['residues'][resname]['atoms']
                # Todo: add termini atoms here to the atom_list
                # ter = self.charmm_instructions['ter'][seg.segname]
                # if res == ter[0]:
                #     atom_list += self.top_content['residues'][resname]['atoms']

                for atm_ref in atom_list:
                    if not res.has_key( atm_ref['name'] ) or (9999 in res[atm_ref['name']]['coord']):
                        # Count missing hydrogens seperatly.
                        if atm_ref['name'][0] == 'H':
                            if res.resname == 'TIP3':
                                # It is a water hydrogen.
                                missing_hydrogens_water += 1
                            elif res.hetatm:
                                # It is a ligand hydrogen.
                                missing_hydrogens_ligands += 1
                            else:
                                # It is a protein hydrogen.
                                # pass

                                # The HN atom in the N-terminal residue is deleted by the NTER patch.
                                nter_res = self.charmm_instructions['ter'][seg.segname][0]
                                if nter_res is not None:
                                    if res == nter_res and atm_ref['name'] == 'HN':
                                        continue
                                    if res in disu_cys_list and atm_ref['name'] == 'HG1':
                                        continue
                                    missing_hydrogens_protein += 1

                            continue

                        # The O-atom in CTER residue will be renamed through a patch in charmm, so it is ok, that this
                        # one is missing here.
                        if atm_ref['name'] == 'O':
                            ter = self.charmm_instructions['ter'][seg.segname]
                            if res in ter:
                                continue

                        #Nucleic Bases additons
                        if atm_ref['name'] == 'P':
                            ter = self.charmm_instructions['ter'][seg.segname]
                            if res in ter:
                                continue
                        if atm_ref['name'] == 'O1P':
                            ter = self.charmm_instructions['ter'][seg.segname]
                            if res in ter:
                                continue
                        if atm_ref['name'] == 'O2P':
                            ter = self.charmm_instructions['ter'][seg.segname]
                            if res in ter:
                                continue

                        #common patching in case of NA is deoxy patch, therefore
                        if atm_ref['name'] == "O2'":
                            deoxy = False
                            for patch, patch_residues in self.charmm_instructions['patches']:
                                if (patch == 'DEO1') or (patch == 'DEO2'):
                                    patch_resi = patch_residues[0]
                                    tuple = (res.resname, res.resid, res.segname)
                                    if tuple == patch_resi:
                                        if not deoxy:
                                            deoxy = True

                            if deoxy:
                                continue


                        # Count missing atoms.
                        missing_atoms_list.append( (res.segname, res.resid, atm_ref['name']) )
                        res_id = str(res)
                        if res.hetatm:
                            # It is a ligand atom.
                            missing_atoms_ligands += 1
                            if not res_id in missing_atoms_ligands_list:
                                missing_atoms_ligands_list[ res_id ] = []
                            missing_atoms_ligands_list[ res_id ].append( atm_ref['name'] )
                        else:
                            # It is a protein atom.
                            missing_atoms_protein += 1
                            if not res_id in missing_atoms_protein_list:
                                missing_atoms_protein_list[ res_id ] = []
                            missing_atoms_protein_list[ res_id ].append( atm_ref['name'] )

        if missing_atoms_protein or missing_atoms_ligands:

            # OLD VERSION -> WORKING WELL
            # decision_name = 'missing_atoms'
            # # decision = self.__get_decision(decision_name)
            # # The 'keep' choice is not available at the moment.
            # decision = 'model'
            #
            # # Set default
            # if decision is None:
            #     decision = 'model'
            #
            # text = u"# There are %i protein and %i ligand atoms missing. Should they be modelled?:\n" % \
            #             ( missing_atoms_protein, missing_atoms_ligands )
            # if missing_atoms_protein:
            #     first = True
            #     # text1 = "  Missing protein atoms: "
            #     text1 = "  Missing atoms: "
            #     text += text1
            #     for res_id, atm_list in missing_atoms_protein_list.iteritems():
            #         text2 = "%s: " % res_id
            #         if not first:
            #             text2 = len(text1)*' ' + text2
            #         first = False
            #         for atm_name in atm_list:
            #             text2 += "%3s " % atm_name
            #         text += text2 + '\n'
            # if missing_atoms_ligands:
            #     first = True
            #     text1 = "  Missing ligand atoms: "
            #     text += text1
            #     for res_id, atm_list in missing_atoms_ligands_list.iteritems():
            #         text2 = "%s: " % res_id
            #         if not first:
            #             text2 = len(text1)*' ' + text2
            #         first = False
            #         for atm_name in atm_list:
            #             text2 += "%3s " % atm_name
            #         text += text2 + '\n'
            #
            #
            # text += "  Decision: " + decision_name + "\n"
            # text += "  Options:  model - Model all missing atoms. (default)\n"
            # text += "            [keep  - Do not model any (non hydrogen) atoms.] option currently not available\n"
            # text += "  Chosen:   " + decision
            #
            #
            # self.charmm_instructions['missing'] = decision
            # self.charmm_instructions['do_minimize'] = True
            #
            # self.charmm_config['missing_atom_ids'] = missing_atoms_list
            #
            # self.decided.append(text)

            decision_name = 'missing_atoms'
            decision = self.__get_decision(decision_name)
            # The 'keep' choice is not available at the moment.

            # Set default
            if decision is None:
                decision = 'model'

            text = u"# There are %i protein and %i ligand atoms missing. Should they be modelled?:\n" % \
                        ( missing_atoms_protein, missing_atoms_ligands )

            if decision == 'model':
                if missing_atoms_protein:
                    first = True
                    # text1 = "  Missing protein atoms: "
                    text1 = "  Missing atoms: "
                    text += text1
                    for res_id, atm_list in missing_atoms_protein_list.iteritems():
                        text2 = "%s: " % res_id
                        if not first:
                            text2 = len(text1)*' ' + text2
                        first = False
                        for atm_name in atm_list:
                            text2 += "%3s " % atm_name
                        text += text2 + '\n'
                if missing_atoms_ligands:
                    first = True
                    text1 = "  Missing ligand atoms: "
                    text += text1
                    for res_id, atm_list in missing_atoms_ligands_list.iteritems():
                        text2 = "%s: " % res_id
                        if not first:
                            text2 = len(text1)*' ' + text2
                        first = False
                        for atm_name in atm_list:
                            text2 += "%3s " % atm_name
                        text += text2 + '\n'
                self.charmm_instructions['do_minimize'] = True
                self.charmm_config['missing_atom_ids'] = missing_atoms_list

            elif decision == 'keep':
                self.charmm_instructions['do_minimize'] = False
            else:
                error = 'Unknown decision choice.'
                raise AssertionError(error)

            text += "  Decision: " + decision_name + "\n"
            text += "  Options:  model - Model all missing atoms. (default)\n"
            text += "            [keep  - Do not model any (non hydrogen) atoms.] option currently not available\n"
            text += "  Chosen:   " + decision

            self.charmm_instructions['missing'] = decision
            self.decided.append(text)


        if missing_hydrogens_water or missing_hydrogens_protein or missing_hydrogens_ligands:

            decision_name = 'missing_hydrogens'
            # decision = self.__get_decision(decision_name)
            # The 'keep' choice is not available at the moment.
            decision = 'model'

            # Set default
            if decision is None:
                decision = 'model'

            text = u"# There are %i protein, %i ligand and %i water hydrogens missing. Should they be modelled?:\n" % \
                        ( missing_hydrogens_protein, missing_hydrogens_ligands, missing_hydrogens_water )

            text += "  Decision: " + decision_name + "\n"
            text += "  Options:  model - Model and optimize all hydrogens. (default)\n"
            text += "            [keep  - Do not model any hydrogen.] option currently not available\n"
            text += "  Chosen:   " + decision

            if decision == 'model':
                self.charmm_instructions['do_minimize'] = True
            self.charmm_instructions['missing'] = decision

            self.decided.append(text)

        # Todo: Check for missing termini

#        print missing_hydrogens_water
#        print missing_hydrogens_protein
#        print missing_hydrogens_ligands
#        print missing_atoms_protein
#        print missing_atoms_ligands
#        pprint(missing_atoms_protein_list)
#        pprint(missing_atoms_ligands_list)

        return True

    def __generate_charmm_input_script(self):
        if self.charmm_out_prefix is None:
            output_prefix = self.title + '_out'
            self.charmm_out_prefix = output_prefix
        else:
            output_prefix = self.charmm_out_prefix

        pdb = self.structure

        # Write the input files for charmm and store them, so they can be deleted after the charmm run.
        self.charmm_tmp_files = []
        for seg in pdb.struct.iter_segments():
            filename = self.workdir + self.title + '_in_' + seg.segname + '.crd'
            self.structure.write_crd(filename, segname=seg.segname)
            self.charmm_tmp_files.append(filename)


        #####################
        ### Write header. ###
        #####################
        script = "! Input script for Charmm generated by 'Charmm_manager'.\n\n"

        script += 'dimension chsize 999999'
        script += '\n'

        ##########################################
        ### Read topology and parameter files. ###
        ##########################################
        append = ''
        for top in self.top:
#            script += "read rtf card name \""  + top + '\" \n' append
            script += 'read rtf card name "%s" %s\n' % (top, append)
            append = 'append'
        append = ''
        for par in self.par:
#            script += "read para card name \"" + par + '\"\n'
            script += 'read para card name "%s" %s\n' % (par, append)
            append = 'append'
        script += '\n'

        #jd_test
        if self.stream:
            script += 'stream %s' % (self.stream)
        script += '\n'


        first_segment_coor = True
        for task in self.charmm_config['tasks']:
            #################################################################################
            ### Read sequence, generate segment and read coordinates for not water chains ###
            #################################################################################
            if task == 'read_structure_sequence':
                # first_segment = True
                for seg in pdb.struct.iter_segments():
                    if seg.get_first_atom()['resname'] == 'TIP3':
                        continue

                    # Read sequence from crd file.
                    script += "read sequence coor name \"" + self.title + '_in_' + seg.segname + '.crd" resid \n'

#                    # Generate setup.
#                    script += 'generate setup ' + seg.segname
#                    if seg.get_first_atom()['resname'] == 'TIP3':
#                        script += ' noangle nodihedral\n'
#                    else:
#                        script += '\n'

                    # Generate setup.
                    if seg.segname in self.charmm_instructions['patch_ter']:
                        # Termini are patched explicitly later.
                        script += 'generate setup %s first none last none' % seg.segname + '\n'

                    elif self.charmm_instructions['ignore_termini']:
                        print 'WARNING: Chain %s will have no termini patches, be AWARE!'
                        script += 'generate setup %s' % seg.segname + '\n'
                    else:
                        nter, cter = self.charmm_instructions['ter'][seg.segname]

                        if (nter is None):
                            nter_patch = 'NONE'
                        elif is_na_base(nter.resname):
                                nter_patch = '5TER'

                        elif ((nter.resid != 1) or ('CAY' in nter)) and self.charmm_instructions['cap_termini']:
                            if nter.resname == 'PRO':
                                nter_patch = 'ACP'
                            else:
                                nter_patch = 'ACE'
                        else:
                            if nter.resname == 'GLY':
                                nter_patch = 'GLYP'
                            elif nter.resname == 'PRO':
                                nter_patch = 'PROP'
                            else:
                                nter_patch = 'NTER'

                        if cter is None:
                            cter_patch = 'NONE'
                        elif is_na_base(cter.resname):
                                cter_patch = '3TER'
                        elif (not ((('OXT' in cter) or ('OT1' in cter) and not ('HT1' in cter)))) and self.charmm_instructions['cap_termini']:
                            if cter.resname == 'PRO':
                                cter_patch = ''
                            else:
                                cter_patch = 'CT1'
                        else:
                            if cter.resname == 'PRO':
                                cter_patch = ''
                            else:
                                cter_patch = 'CTER'

                        if cter_patch != '':
                            script += 'generate setup %s first %s last %s\n' % (seg.segname, nter_patch, cter_patch)
                        else:
                            script += 'generate setup %s first %s\n' % (seg.segname, nter_patch)

                    # # Read coordinates from crd file.
                    # script += "read coor card name \"" + self.title + '_in_' + seg.segname + '.crd"'
                    # if not first_segment:
                    #     script += ' append\n'
                    # else:
                    #     script += '\n'
                    #     first_segment = False

                    script += '\n'

            #############################################################################
            ### Read sequence, generate segment and read coordinates for water chains ###
            #############################################################################
            elif task == 'read_water_sequence':
                for seg in pdb.struct.iter_segments():
                    if seg.get_first_atom()['resname'] != 'TIP3':
                        continue

                    # Read sequence from crd file.
                    script += "read sequence coor name \"" + self.title + '_in_' + seg.segname + '.crd" resid \n'

                    # Generate setup.
                    script += 'generate setup ' + seg.segname + ' noangle nodihedral\n'

                    # Read coordinates from crd file.
                    # script += "read coor card name \"" + self.title + '_in_' + seg.segname + '.crd"'
                    # script += ' append\n'

                    script += '\n'

            #######################################
            ### Read coordinates for all chains ###
            #######################################
            elif task == 'read_coordinates':
                # This is now a "global" variable in this function.
                # first_segment_coor = True
                for seg in pdb.struct.iter_segments():
                    if seg.get_first_atom()['resname'] == 'TIP3':
                        continue

                    # Read coordinates from crd file.
                    script += "read coor card name \"" + self.title + '_in_' + seg.segname + '.crd"'
                    if not first_segment_coor:
                        script += ' resid\n'
                    else:
                        script += '\n'
                        first_segment_coor = False
                script += '\n'

            elif task == 'read_water_coordinates':
                # This is now a "global" variable in this function.
                # first_segment_coor = True
                for seg in pdb.struct.iter_segments():
                    if seg.get_first_atom()['resname'] != 'TIP3':
                        continue

                    # Read coordinates from crd file.
                    script += "read coor card name \"" + self.title + '_in_' + seg.segname + '.crd"'
                    if not first_segment_coor:
                        script += ' resid\n'
                    else:
                        script += '\n'
                        first_segment_coor = False
                script += '\n'

            ##################################
            ### Model gaps in the sequence ###
            ##################################
            elif task == 'model_gaps':
                for gap in self.charmm_instructions['gaps']:
                    if gap[0] == 'connect':
                        start_res = gap[1]
                        end_res   = gap[2]
                        patch_name = 'LINK' # ToDo: Test this feature!
                        id1 = '%s %s' % ( start_res.segname, str(start_res.resid) )
                        id2 = '%s %s' % ( end_res.segname,   str(end_res.resid) )
                        script += "patch %s %s %s\n" % (patch_name, id1, id2)
                        script += '\n'

            ###########################
            ### Build missing atoms ###
            ###########################
            elif task == 'build_missing':
                if self.charmm_instructions['missing'] == 'model':
                    script += 'ic para\n'
                    script += 'ic build\n'
                    script += '\n'

            #################################
            ### Patch disulphide bridges. ###
            #################################
            elif task == 'patch_disu':
                for (cys1, cys2) in self.charmm_instructions['disu']:
                    id1 = '%s %s' % ( cys1.segname, str(cys1.resid) )
                    id2 = '%s %s' % ( cys2.segname, str(cys2.resid) )
                    script += 'patch disu %s %s\n' % (id1, id2)
                    self.charmm_config['autogen'] = True
                if self.charmm_instructions['disu']:
                    script += '\n'

            ######################
            ### Patch termini. ###
            ######################
            elif task == 'patch_ter':
                for segname in self.charmm_instructions['ter'].keys():
                    if segname in self.charmm_instructions['patch_ter']:
                        (nter, cter) = self.charmm_instructions['ter'][segname]
                        id_nter = '%s %s' % ( nter.segname, str(nter.resid) )
                        id_cter = '%s %s' % ( cter.segname, str(cter.resid) )
                        script += 'patch nter %s\n' % id_nter
                        script += 'patch cter %s\n' % id_cter
                        script += '\n'
                        # Call 'AUTOgenerate ANGLes DIHEdrals' before water chains are created.
                        self.charmm_config['autogen'] = True

            ###############
            ### Patches ###
            ###############
            elif task == 'patches':
                for patch_name, residue_list in self.charmm_instructions['patches']:
                    residue_descr = ''
                    if type(residue_list) is not list:
                        residue_list = [residue_list]
                    for residue_tuple in residue_list:
                        if len(residue_tuple) == 2:
                            resid, segname = residue_tuple
                        elif len(residue_tuple) == 3:
                            resname, resid, segname = residue_tuple
                        else:
                            error = "Entry %s charmm_instructions['patches'] has not the correct format. It should be" \
                                    "tuple of length 2 or 3"
                            raise(AssertionError(error))
                        residue_descr += ' %s %i' % (segname, resid)
                    script += 'patch %s %s setup\n' % (patch_name, residue_descr)

                if self.charmm_instructions['patches']:
                    script += '\n'
                    # Call 'AUTOgenerate ANGLes DIHEdrals' before water chains are created.
                    self.charmm_config['autogen'] = True

            #####################################
            ### AUTOgenerate ANGLes DIHEdrals ###
            #####################################
            elif task == 'autogen':
                if self.charmm_config['autogen']:
                    script += 'AUTOgenerate ANGLes DIHEdrals\n'
                    script += '\n'

            ##################################
            ### Patches after autogenerate ###
            ##################################
            elif task == 'patches_no_autogen':
                for patch_name, residue_list in self.charmm_instructions['patches_no_autogen']:
                    residue_descr = ''
                    if type(residue_list) is not list:
                        residue_list = [residue_list]
                    for residue_tuple in residue_list:
                        if len(residue_tuple) == 2:
                            resid, segname = residue_tuple
                        elif len(residue_tuple) == 3:
                            resname, resid, segname = residue_tuple
                        else:
                            error = "Entry %s charmm_instructions['patches'] has not the correct format. It should be" \
                                    "tuple of length 2 or 3"
                            raise(AssertionError(error))
                        residue_descr += ' %s %i' % (segname, resid)
                    script += 'patch %s %s setup\n' % (patch_name, residue_descr)

            ##############
            ### HBUILD ###
            ##############
            elif task == 'hbuild':
                script += 'hbuild\n'
                script += '\n'

            ##############################################
            ### Minimize modelled parts of the protein ###
            ##############################################
            elif task == 'minimize_modeled':
                if self.charmm_instructions['do_minimize']:
                    # Create a selection for all modelled single atoms.
                    atm_sel = ''
                    sel_def = ''
                    first_one = True
                    c = 1
                    c_sel = 1
                    c_sel_def = 1
                    c_sel_def_linelength = 1
                    for resname, resid, name in self.charmm_config['missing_atom_ids']:
                        # if not first_one:
                        #     atm_sel += " .or. "
                        # else:
                        #     first_one = False
                        if len(sel_def) > 0:
                            sel_def += " .or. "
                        c += 1
                        c_sel += 1

                        # atm_id = "%s_%s_%s" % (resname, resid, name)
                        # def_str += "define %s select segid %s .and. resid %s .and. type %s end\n" % (atm_id, resname, resid, name)
                        if c <= 2:
                            # atm_sel += "(segid %s .and. resid %s .and. type %s)" % (resname, resid, name)
                            sel_def += "(segid %s .and. resid %s .and. type %s)" % (resname, resid, name)
                            # atm_sel += "atm_id"
                        else:
                            # atm_sel += "(segid %s .and. resid %s .and. type %s) -\n               " % (resname, resid, name)
                            sel_def += "(segid %s .and. resid %s .and. type %s) -\n              " % (resname, resid, name)
                            # atm_sel += "atm_id -\n               "
                            c = 1

                        if c_sel > 40:
                            script += "define sel_%i select %s end\n" % (c_sel_def, sel_def)
                            if len(atm_sel) > 0:
                                atm_sel +=  " .or. sel_%i" % c_sel_def
                            else:
                                atm_sel +=  "sel_%i" % c_sel_def
                            sel_def = ''
                            c_sel_def += 1
                            c_sel = 1
                            c_sel_def_linelength += 1
                            if c_sel_def_linelength > 5:
                                atm_sel += " -\n                  "
                                c_sel_def_linelength = 0

                    if len(sel_def) > 0:
                        script += "define sel_%i select %s end\n" % (c_sel_def, sel_def)
                        if len(atm_sel) > 0:
                            atm_sel +=  " .or. sel_%i" % c_sel_def
                        else:
                            atm_sel +=  "sel_%i" % c_sel_def

                    # Add selections for all modelled gap residues.
                    for gap in self.charmm_instructions['gaps']:
                        if gap[0] != 'fill':
                            continue
                        start_resid = gap[1].resid
                        end_resid   = gap[2].resid
                        segname     = gap[1].segname
                        #if not first_one:
                        if len(atm_sel) > 0:
                            atm_sel += " .or. "
                        else:
                            first_one = False
                        c += 1
                        if c <= 2:
                            atm_sel += "(segid %s .and. resid %i:%i)" % (segname, start_resid+1, end_resid-1)
                        else:
                            atm_sel += "(segid %s .and. resid %i:%i) -\n               " % (segname, start_resid+1, end_resid-1)
                            c = 1

                    # Add selections for specified residues.
                    for residues in self.charmm_instructions['minimize']:
                        if type(residues) is list:
                            start_resid, end_resid = residues
                            residue = None
                        else:
                            residue = residues

                        if not first_one:
                            atm_sel += " .or. "
                        else:
                            first_one = False
                        c += 1
                        if c <= 2:
                            if residues is not None:
                                atm_sel += "(segid %s .and. resid %i)" \
                                           % (residue.segname, residue.resid)
                            else:
                                atm_sel += "(segid %s .and. resid %i:%i)" \
                                           % (start_resid.segname, start_resid.resid, end_resid.resid)
                        else:
                            if residues is not None:
                                atm_sel += "(segid %s .and. resid %i) -\n               " \
                                           % (residue.segname, residue.resid)
                            else:
                                atm_sel += "(segid %s .and. resid %i:%i) -\n               " \
                                           % (start_resid.segname, start_resid.resid, end_resid.resid)
                            c = 1

                    for selstr in self.charmm_instructions['minimize_selections']:
                        if len(atm_sel) > 0:
                            atm_sel += ' .or. (%s) ' % selstr
                            # first_one = False
                        else:
                            atm_sel += ' (%s) ' % selstr


                    if not atm_sel:
                        atm_sel = "none"

                    if self.charmm_instructions['backbone_fixed']:
                        script += "define backbone select type n .or. type ca .or. type c end\n"
                        script += "cons fix select .not. ( %s .or. hydrogens ) .or. backbone end\n" % atm_sel
                    else:
                        script += "cons fix select .not. ( %s .or. hydrogens )   end\n" % atm_sel
                    script += "minimize sd nsteps 1000 tolg 0.1\n"
                    script += "minimize abnr nsteps 5000 tolg 0.01\n"
                    script += "\n"




            ############################
            ### Write out structures ###
            ############################
            elif task == 'write_structure':
                script += "ioformat exten"  + '\n'
                for file_type in self.charmm_config['output']:
                    if file_type == 'psf':
                        script += "write psf card name \"" + output_prefix + '.psf"\n'

                    if file_type == 'xplor.psf':
                        script += "write psf xplor card name \"" + output_prefix + '.xplor.psf"\n'

                    if file_type == 'crd':
                        script += "write coord card name \"" + output_prefix + '.crd"\n'

                    if file_type == 'pdb':
                        script += "write coord pdb name \"" + output_prefix + '.pdb"\n'
                script += "\n"

#                self.charmm_config['output'] = ['pdb', 'crd', 'psf', 'xplor.psf']
#                    write coord pdb name "output.pdb"
            elif task == '':
                pass
            else:
                script += task
                script += "\n"

        script += "stop\n"

        return script

    def __remove_tmp_files(self):
        for fn in self.charmm_tmp_files:
            os.remove(fn)

    def run_charmm(self, submit=False, dolly=False):

        (permission, message) = self.__check_status('run_charmm')

        decision_str = ''
        if len(self.decided) > 0:
            decision_str += "\n### The following changes have been applied to the structure. ###\n"
            for decided in self.decided:
                decision_str += decided + '\n\n'
        if len(self.undecided) > 0:
            decision_str += "\n### There are undecided decisions for whose no default behaviour exists. ###\n"
            for undecided in self.undecided:
                decision_str += undecided + '\n\n'
        decision_filename = self.workdir +  'modelling_decision.dat'
        f = open(decision_filename, 'w')
        f.write(decision_str)
        f.close()

        script = self.__generate_charmm_input_script()

        script_name = self.title + '_charmm.inp'
        filename = self.workdir + script_name

        f = open(filename, 'w')
        f.write(script)
        f.close()


        from subprocess import call

        output_filename = self.workdir +  self.title + '_charmm.out'

        if submit:
            # print 'CHARMM job will be submited to Dolly Cluster. Waiting for completion is ON! Never trust a sheep.'
            # submit_charmm(self.workdir, script_name, self.title + '_charmm.out', wait_to_complete=True)
            submit_charmm(self.workdir, script_name, self.title + '_charmm.out', wait_to_complete=False)
            # submit_charmm(self.workdir, script_name, self.title + '_charmm.out')
            return
        elif dolly:
            commands  = "tcsh;"
            # commands  = "csh;"
            commands  = "cd " + self.workdir + ";"
            # commands += '/usr/local/dolly/bin/charmm36' + " -input " + script_name + " > %s" % output_filename
            commands += 'charmm36' + " < " + script_name + " > %s" % output_filename
            call(commands, shell=True)
        else:
            commands  = "tcsh;"
            commands  = "cd " + self.workdir + ";"
            # commands += '/usr/local/dolly/bin/charmm36' + " -input " + script_name + " > %s" % output_filename
            commands += self.charmm_bin + " < " + script_name + " > %s" % output_filename
            call(commands, shell=True)


        f_out = open(output_filename, 'r')
        self.charmm_normal_termination = False
        for next_line in f_out:
            self.charmm_output += next_line
            if 'NORMAL TERMINATION BY NORMAL STOP' in next_line:
                self.charmm_normal_termination = True
        f_out.close()


        # commands  = "cd " + self.workdir + "\n"
        # commands += self.charmm_bin + " -input " + script_name + "\n"
        #
        # commands += "exit\n"
        # process = subprocess.Popen("tcsh",\
        #                 shell=True,\
        #                 stdin=subprocess.PIPE,\
        #                 stdout=subprocess.PIPE,\
        #                 stderr=subprocess.PIPE,\
        #                 )
        # process.stdin.write(commands)
        #
        # output_filename = self.workdir +  self.title + '_charmm.out'
        # f_out = open(output_filename, 'w')
        # while True:
        #     next_line = process.stdout.readline()
        #
        #     if not next_line:
        #         break
        #
        #     self.charmm_output += next_line
        #     f_out.write(next_line)
        #
        #     if next_line.find('NORMAL TERMINATION BY NORMAL STOP'):
        #         self.charmm_normal_termination = True
        # f_out.close()

        process = None


        if self.charmm_normal_termination:
            if self.clean_up:
               self.__remove_tmp_files()

            return 1
        else:
            return 0

    def set_run_charmm(self):

        '''Instead of run, this fucntion just sets the run!'''

        (permission, message) = self.__check_status('run_charmm')

        decision_str = ''
        if len(self.decided) > 0:
            decision_str += "\n### The following changes have been applied to the structure. ###\n"
            for decided in self.decided:
                decision_str += decided + '\n\n'
        if len(self.undecided) > 0:
            decision_str += "\n### There are undecided decisions for whose no default behaviour exists. ###\n"
            for undecided in self.undecided:
                decision_str += undecided + '\n\n'
        decision_filename = self.workdir +  'modelling_decision.dat'
        f = open(decision_filename, 'w')
        f.write(decision_str)
        f.close()

        script = self.__generate_charmm_input_script()

        script_name = self.title + '_charmm.inp'
        filename = self.workdir + script_name

        f = open(filename, 'w')
        f.write(script)
        f.close()



    def get_modelled_structure(self):
        """
        Returns an object of Simple_structure_parser Class.

        """
        # charmm_out_title =
        output_filename_base = self.workdir + self.title + '_out'
        crd_filename = self.workdir + self.charmm_out_prefix + '.crd'
        xplor_psf_filename = self.workdir + self.charmm_out_prefix + '.xplor.psf'

        modelled_structure = Simple_struct_parser()
        modelled_structure.read_crd(crd_filename)
        modelled_structure.read_xplor_psf(xplor_psf_filename)

        return modelled_structure

    def apply_titr_residue_dict(self, titr_residue_dict):
        """Applies the protonation states as specified in titr_residue_dict to the structure.

        @param: titr_residue_dict: Dict created by self.get_titr_residue_dict.
        """

        # Apply terminal patches at the end, since N-terminus patch may change atoms types.
        resnames_to_delayed = ['NTE', 'CTE']
        delayed_residues = []
        for residue_tuple, state in titr_residue_dict.iteritems():
            if residue_tuple[0] in resnames_to_delayed:
                delayed_residues.append((residue_tuple, state))
            else:
                self.set_prot_residue(residue_tuple, state=state)
        for residue_tuple, state in delayed_residues:
            self.set_prot_residue(residue_tuple, state=state)


    def set_prot_residue(self, residue, state=None, charge=None, patch='undefined', rename='undefined'):
        """
        Changes the state of of the specified residue. Only one of the parameters state, charge, patch and rename is
        allowed.
        One of these three parameters has to be specified.
        The changes are documented in self.prot_state and all instructions are applied. If the residue is not in its
        original protonation state, all previous applied changes are reversed.

        @param residue: Description of the residue. Can be (<resname>, <resid>, <segname>) or
        string <resname>-<resid>_<segname>
        @param state: int: State nr according to entry nr in self.titr_residues.
        @param charge: int: Selects the first state for this residue type that matches the given charge.
        @param patch: None|str: Selects the state according to the applied patch. Can be None.
        @param rename: None|str: Selects the state according to the string applied for renameing. Can be None.
        @return: None
        """
        if type(residue) == str:
            entries = re.split(r'[-_]', residue)
            resname, resid, segname = entries
            residue_tuple = (resname, int(resid), segname)
        elif type(residue) == tuple:
            residue_tuple = residue

        resname, resid, segname = residue_tuple

        if not self.titr_restypes.has_key(resname):
            error = "No information about how to change protonation of residue %s available!" % resname
            raise AssertionError(error)
        restype = self.titr_restypes[resname]

        if (state is None) and (charge is None) and (patch == 'undefined') and (rename == 'undefined'):
            just_reverse = True
        else:
            just_reverse = False
            new_state = self.__find_state(restype, state, charge, patch, rename)

        # Get the current protonation state.
        if self.prot_state.has_key(residue_tuple):
            old_state = self.prot_state[residue_tuple]
        else:
            old_state = None

        if not just_reverse:
            if old_state == new_state:
                # Nothing to do
                return

        if old_state is not None:
            # Reverse previous changes
            self.__reverse_prot_state(residue_tuple)

        if not just_reverse:
            # Apply new state
            self.__apply_prot_state(residue_tuple, new_state)

        return


    def __find_state(self, restype, state=None, charge=None, patch='undefined', rename='undefined'):
        """ Finds the state in self.titr_residues[restype] that matches one of the keywords criteria given by the
        parameter state, charge, patch or rename. Only one of the parameter state, charge, patch and rename is allowed
        to be set.
        @return int: The matching state.
        """

        #Checks if the residue belongs to the ones that are titratable
        if not self.titr_restypes.has_key(restype):
            error = "No information about how to change protonation of residue %s available!" % restype
            raise AssertionError(error)

        restype_ref = self.titr_restypes[restype]

        # Find the new protonation state.
        parameters_set = 0
        if state is not None:
            parameters_set += 1
            if not 0 <= state < len(self.titr_residues[restype_ref]):
                error = "Residue %s has no state %i!" % (str(restype), state)
                raise AssertionError(error)
            new_state = state
        elif charge is not None:
            parameters_set += 1
            for i, state_descr in enumerate(self.titr_residues[restype_ref]):
                if state_descr['charge'] == charge:
                    new_state = i
                    break
            else:
                error = "Residue %s has no state with charge %i!" % (str(restype), charge)
                raise AssertionError(error)
        elif patch != 'undefined':
            parameters_set += 1
            for i, state_descr in enumerate(self.titr_residues[restype_ref]):
                if state_descr['patch'] == patch:
                    new_state = i
                    break
            else:
                error = "Residue %s has no state with patch %s!" % (str(restype), patch)
                raise AssertionError(error)
        elif rename != 'undefined':
            parameters_set += 1
            for i, state_descr in enumerate(self.titr_residues[restype_ref]):
                if state_descr['rename'] == rename:
                    new_state = i
                    break
            else:
                error = "Residue %s has no state with renaming statement %s!" % (str(restype), rename)
                raise AssertionError(error)
        if parameters_set > 1:
            error = "More than one of the parameters 'state', 'charge' and 'patch' has been specified. Only one is " \
                    "allowed!"
            raise(AssertionError(error))
        elif parameters_set == 0:
            error = "Non of the parameters 'state', 'charge' and 'patch' has been specified. One is required!"
            raise(AssertionError(error))

        return new_state

    def __apply_prot_state(self, residue_tuple, state):
        """

        @param residue_tuple: tuple: (<resname>, <resid>, <segname>)
        @param state: int: state nr according to entry nr in self.titr_residues
        @return:
        """
        resname, resid, segname = residue_tuple
        restype = self.titr_restypes[resname]
        state_descr = self.titr_residues[restype][state]

        if state_descr['rename'] is not None:
            self.structure.struct[segname][resid].rename(state_descr['rename'])
        if state_descr['patch'] is not None:
            patch_name = state_descr['patch']
            self.charmm_instructions['patches'].append((patch_name, residue_tuple))
        if state_descr['external_patches'] is not None:
            for patch_name, residue_list in state_descr['external_patches']:
                self.charmm_instructions['patches'].append((patch_name, list(residue_list)))

        # special?

        # Save new state
        self.prot_state[residue_tuple] = state

        return

    def __reverse_prot_state(self, residue_tuple):
        """

        @param residue_tuple: tuple: (<resname>, <resid>, <segname>)
        @return:
        """
        resname, resid, segname = residue_tuple
        restype = self.titr_restypes[resname]
        current_state = self.prot_state[residue_tuple]
        state_descr = self.titr_residues[restype][current_state]

        if state_descr['rename'] is not None:
            # Rename the residue_tuple to its original name.
            self.structure.struct[segname][resid].rename(restype)
        if state_descr['patch'] is not None:
            # Remove previous applied patches from self.charmm_instructions['patches']
            patch_name = state_descr['patch']
            # COMMENT:
            # There could be a try / except ValueError statement here. But if the patch is not in the list
            # there is probably something wrong. Or can this situation occur?
            # It is the same for state_descr['external_patches'] below.
            self.charmm_instructions['patches'].remove((patch_name, residue_tuple))
        if state_descr['external_patches'] is not None:
            for patch_name, residue_list in state_descr['external_patches']:
                self.charmm_instructions['patches'].remove((patch_name, residue_list))

        # Remove state entry
        self.prot_state[residue_tuple] = None

        return


    def get_titr_residue_dict(self, all=False):
        """Returns a dictionary that contains an entry for every residue in the structure that is titratable, according
         to self.titr_residues.

         @returns: dict: {(<resname>, <resid>, <segname>) : None}
        """
        cys_list = []
        for disu_bridge in self.charmm_instructions['disu']:
            for cys in disu_bridge:
                residue_descr = '%s-%i_%s' % (cys.resname, cys.resid, cys.segname)
                cys_list.append(residue_descr)

        titr_residue_dict = {}
        for residue in self.structure.struct.iter_residues():
            resname = residue.resname
            if resname in self.titr_restypes:
                restype = self.titr_restypes[resname]
                resid = residue.resid
                segname = residue.segname

                residue_descr = '%s-%i_%s' % (resname, int(resid), segname)

                if resname == 'CYS':
                    if residue_descr in cys_list:
                        continue

                residue_tuple = (restype, resid, segname)

                if residue_tuple in self.prot_state:
                    state = self.prot_state[residue_tuple]
                    titr_residue_dict[residue_tuple] = state
                else:
                    if all:
                        titr_residue_dict[residue_tuple] = 0
                    else:
                        titr_residue_dict[residue_tuple] = None

        # Add termini
        self.find_termini()

        for segname, (nter, cter) in self.charmm_instructions['ter'].iteritems():
            for restype, residue in zip(['NTE', 'CTE'], [nter, cter]):

                if residue is None:
                    continue
                if restype == 'NTE':
                    if (nter.resid != 1) or ('CAY' in nter):
                        continue
                if restype == 'CTE':
                    if not (('OXT' in cter) or (('OT1' in cter) and not ('HT1' in cter))):
                        continue
                if is_na_base(residue.resname):
                    continue

                resid = residue.resid
                residue_tuple = (restype, resid, segname)
                if residue_tuple in self.prot_state:
                    state = self.prot_state[residue_tuple]
                    titr_residue_dict[residue_tuple] = state
                else:
                    if all:
                        titr_residue_dict[residue_tuple] = 0
                    else:
                        titr_residue_dict[residue_tuple] = None

        return titr_residue_dict

    @staticmethod
    def copy_titr_residue_dict(titr_residue_dict):
        """
        Returns a full copy of residue titr_residue_dict
        This should be used in case of modification of this dictionary if the original should be kept.
        """
        titr_residue_dict_copy = {}
        for residue_tuple, state in titr_residue_dict.iteritems():
            titr_residue_dict_copy[tuple(residue_tuple)] = state
        return titr_residue_dict_copy


    def mod_titr_residue_dict(self, titr_residue_dict, restype, state=None, charge=None, patch='undefined',
                              rename='undefined'):
        """Modifies titr_residue_dict that is a dict created by self.get_titr_residue_dict(). All entries that
        match the given resname set to the state nr, as defined by the keywords state, charge, patch or rename and
        self.titr_residues. Only one of the keywords state, charge, patch and rename is allowed.

        @param: dict: Created by self.get_titr_residue_dict()
        @param restype: Residue name. There must be an entry with the same name in self.titr_restypes
        @param state: int: State nr according to entry nr in self.titr_residues.
        @param charge: int: Selects the first state for this residue type that matches the given charge.
        @param patch: None|str: Selects the state according to the applied patch. Can be None.
        @param rename: None|str: Selects the state according to the string applied for renaming. Can be None.
        @return: None
        """

        restype_ref = self.titr_restypes[restype]

        new_state = self.__find_state(restype_ref, state, charge, patch, rename)

        for residue_tuple in titr_residue_dict.keys():
            resname, resid, segname = residue_tuple
            restype_comp = self.titr_restypes[resname]
            if restype_comp == restype_ref:
                titr_residue_dict[residue_tuple] = new_state

        return

    def mod_titr_residue_dict_single(self, titr_residue_dict, residue_tuple, state=None, charge=None, patch='undefined',
                              rename='undefined'):
        """Modifies titr_residue_dict that is a dict created by self.get_titr_residue_dict(). The entry that
        matches the given resname is set to the state nr, as defined by the keywords state, charge, patch or rename and
        self.titr_residues. Only one of the keywords state, charge, patch and rename is allowed.

        @param: dict: Created by self.get_titr_residue_dict()
        @param residue_tuple: Tuple in the format (<resname>, <resid>, <segname>)
        @param state: int: State nr according to entry nr in self.titr_residues.
        @param charge: int: Selects the first state for this residue type that matches the given charge.
        @param patch: None|str: Selects the state according to the applied patch. Can be None.
        @param rename: None|str: Selects the state according to the string applied for renaming. Can be None.
        @return: None
        """

        resname, resid, segname = residue_tuple
        restype = self.titr_restypes[resname]
        new_state = self.__find_state(restype, state, charge, patch, rename)
        residue_tuple_ref = (restype, resid, segname)
        titr_residue_dict[residue_tuple] = new_state

        return

    def mod_titr_residue_dict_from_titrable_yaml(self, titr_residue_dict, prot_state, titratable_yaml, restrict=None,
                                                 residue_list=None):
        """Modifies titr_residue_dict that is a dict created by self.get_titr_residue_dict by transferring the
        information from 'prot_state'. 'prot_state' is a definition of the protonation states, created by
        kbp_tools.determine_md_protonation_pattern.

        @param: dict: Created by self.get_titr_residue_dict()
        @param prot_state: States defined by kbp_tools.determine_md_protonation_pattern.
        @param titratable_yaml: state definition created by kbp_tools.parse_titratable_yaml
        @param: Todo: restrict_restype, restrict_charge = restrict
        @param residue_list: list of residues described by residue tuple: (resname, resid, segname). Only these
                             residues are modified.

        @return: None
        """
        if restrict is not None:
            restrict_restype, restrict_charge = restrict

        for residue_descr, state_yaml in prot_state.iteritems():
            resname, resid, segname = re.split(r'[-_]', residue_descr)
            resid = int(resid)
            resname_std = kbp_tools.get_real_resname(resname)
            restype = self.titr_restypes[resname_std]
            residue_tuple = (restype, resid, segname)

            if titratable_yaml[resname][state_yaml].has_key('patch'):
                patch = titratable_yaml[resname][state_yaml]['patch']
            else:
                patch = None


            # There are some differences in the patches/renaming used for the MD and by Karlsberg+.
            rename = None
            patch_md = patch
            # 1) Check for patches that are done by renaming.
            if patch == 'HSPD':
                rename = 'HSD'
                patch_md = None
            elif patch == 'HSPE':
                rename = 'HSE'
                patch_md = None
            elif (patch == None) and (resname == 'HSP'):
                rename = 'HSP'
                patch_md = None
            # 2) Check for differences in the used patches.
            elif patch_md is not None:
                # patch_md =    patch.replace('LYSP', 'LYSR')
                # patch_md = patch_md.replace('ARGP', 'ARGR')
                patch_md = patch_md.replace('DPP1', 'ASPP')
                patch_md = patch_md.replace('DPP2', 'ASPP')
                patch_md = patch_md.replace('EPP1', 'GLUP')
                patch_md = patch_md.replace('EPP2', 'GLUP')

            if rename is not None:
                state = self.__find_state(restype, rename=rename)
            else:
                state = self.__find_state(restype, patch=patch_md)

            # state_definition = self.titr_residues[restype][state]
            # if state_definition['patch'] is None \
            #         and state_definition['rename'] is None \
            #         and 'external_patches' is None:

            # if state == 0:
                # This residue is in its reference state.
                # titr_residue_dict[residue_tuple] = None
            # else:
            #     titr_residue_dict[residue_tuple] = state

            if restrict is not None:
                if restrict_restype != restype:
                    continue
                if restrict_charge is not None:
                    if self.titr_residues[restype][state]['charge'] != restrict_charge:
                        continue

            if residue_list is not None:
                if not (resname_std, resid, segname) in residue_list:
                    continue

            titr_residue_dict[residue_tuple] = state

        return

    def get_titr_residue_list(self, titr_residue_dict):
        """Returns a list of residues in titr_residue_dict.
        @param titr_residue_dict: dict created by self.get_titr_residue_dict()
        @return: residue_list: list: [str: <resname>_<resid>_<segname> e.g. 'ARG-24_A']
        """
        cys_list = []
        for disu_bridge in self.charmm_instructions['disu']:
            for cys in disu_bridge:
                residue_descr = '%s-%i_%s' % (cys.resname, cys.resid, cys.segname)
                cys_list.append(residue_descr)

        residue_list = []
        for residue_tuple in titr_residue_dict.keys():
            resname, resid, segname = residue_tuple
            residue_descr = '%s-%i_%s' % (resname, resid, segname)
            if resname == 'CYS':
                if residue_descr in cys_list:
                    continue
            residue_list.append(residue_descr)

        return residue_list

    def rename_residue(self, residue_tuple, new_resname, patch_translation_dict={}):
        """Renames the specified residue. Changes are applied to:
                - The Simple_structure_parser object self.structure
                - self.charmm_instructions['patches']
                - self.charmm_instructions['minimize']
                - self.charmm_instructions['ter']
                - self.charmm_instructions['disu']

                - self.prot_state
                  If the residues protonation state has been modified using the 'set_prot_residue' function, then there
                  must be an entry in self.titr_residues, for the new resname.

        @param residue_tuple: residue tuple: (resname, resid, segname)
        @param new_resname: The new residue name.
        @return:
        """
        resname, resid, segname = residue_tuple
        new_residue_tuple = new_resname, resid, segname

        # Change residue names in self.structure
        self.structure[segname][resid].rename(new_resname)

        # Change self.prot_state
        if residue_tuple in self.prot_state:
            # Check that there is a entry for the new resname
            if not new_resname in self.titr_residues:
                error = "Protonation state of residue %s has been modified previously, but there is no information in "\
                        "'self.titr_residues' about how this protonation state can be created for the new resname %s." \
                        % (str(residue_tuple), new_resname)
                raise AssertionError(error)

            prot_state = self.prot_state[residue_tuple]

            self.__reverse_prot_state(residue_tuple)
            self.__apply_prot_state(new_residue_tuple, prot_state)

            self.prot_state[new_residue_tuple] = prot_state
            self.prot_state.pop(residue_tuple)

        # Check self.charmm_instructions['patches']
        short_residue_tuple = (resid, segname)
        for patch_index, (patch_name, residue_tuple_patch) in enumerate(self.charmm_instructions['patches']):
            matching_patch = False
            if len(residue_tuple_patch) == 2:
                if residue_tuple_patch == short_residue_tuple:
                    matching_patch = True
            elif len(residue_tuple_patch) == 3:
                if residue_tuple_patch == residue_tuple:
                    matching_patch = True
            if matching_patch:
                if patch_name not in patch_translation_dict:
                    error = "Residue %s was patched with %s, but there is no information in " \
                            "'patch_translation_dict' about how this patch should be renamed for the new " \
                            "resname %s." \
                            % (str(residue_tuple), patch_name, new_resname)
                    raise AssertionError(error)
                else:
                    new_patch_name = patch_translation_dict[patch_name]
                    self.charmm_instructions['patches'][patch_index] = (new_patch_name, residue_tuple)

        # Todo: Check the following entries
        # - self.charmm_instructions['minimize']
        # - self.charmm_instructions['ter']
        # - self.charmm_instructions['disu'] -> crash




    # def change_prot(self, resid, chain_seg=None, ter_type=None):
    #     pdb = self.structure
    #
    #     res = pdb.struct[chain_seg][resid]
    #     if ter_type is None:
    #         resname = res.resname
    #     else:
    #         resname = ter_type
    #     if resname == 'TYR':
    #         patch_name = 'TYRD'
    #     elif resname == 'ARG':
    #         # ARGR for implicit deprotonation.
    #         patch_name = 'ARGR'
    #         # patch_name = 'RPP1'
    #     elif resname == 'LYS':
    #         # Implicit deprotonation.
    #         patch_name = 'LYSR'
    #     elif resname == 'ASP':
    #         patch_name = 'ASPP'
    #     elif resname == 'GLU':
    #         patch_name = 'GLUP'
    #     elif resname == 'CYS':
    #         patch_name = 'CYSD'
    #     elif resname == 'NTE':
    #         patch_name = 'NTEREF'
    #     elif resname == 'CTE':
    #         patch_name = 'CTEREF'
    #     elif resname == 'DPP':
    #         patch_name = 'DPP2'
    #     elif resname == 'EPP':
    #         patch_name = 'EPP2'
    #     else:
    #         self.GeneralError("Error in 'change_prot': Unknown resname %s" % resname)
    #
    #     charmm_id = '%s %s' % ( res.segname, str(res.resid) )
    #     charmm_command = 'patch %s %s\n' % (patch_name, charmm_id)
    #
    #     for i, entry in enumerate(self.charmm_config['tasks']):
    #         if entry == 'patch_disu':
    #             self.charmm_config['tasks'].insert(i+1, charmm_command)
    #             break
    #
    #     self.charmm_config['autogen'] = True
    #
    #     # Todo: Es ermöglichen Termini zu patchen.

    # def change_type_prot(self, resname, exclude_list=[], charged=None):
    #     pdb = self.structure
    #     recreate_struct = False
    #
    #     # For Cysteins make sure that they are not part of a disulfid bridge,
    #     if resname == 'CYS':
    #         disu_cys_list = []
    #         for (cys1, cys2) in self.charmm_instructions['disu']:
    #             disu_cys_list.append(cys1)
    #             disu_cys_list.append(cys2)
    #
    #     if (resname not in ['HIS', 'CTE', 'NTE']) and (charged != True):
    #         for seg in pdb.struct.iter_segments():
    #             for res in seg.iter_residues():
    #                 if resname == res.resname:
    #                     resid = res.resid
    #                     segname = seg.segname
    #                     if resname == 'CYS'and res in disu_cys_list:
    #                         continue
    #                     else:
    #                         residue = "%s-%i_%s" % (resname, resid, segname)
    #                         if not residue in exclude_list:
    #                             self.change_prot(resid, segname)
    #                         # else:
    #                         #    print "Skipping residue " + residue
    #     elif resname == 'CTE':
    #         for seg in pdb.struct.iter_segments():
    #             (nter, cter) = self.charmm_instructions['ter'][seg.segname]
    #             if cter is not None:
    #                 self.change_prot(cter.resid, seg.segname, ter_type='CTE')
    #     elif resname == 'NTE':
    #         for seg in pdb.struct.iter_segments():
    #             (nter, cter) = self.charmm_instructions['ter'][seg.segname]
    #             if nter is not None:
    #                 self.change_prot(nter.resid, seg.segname, ter_type='NTE')
    #     elif resname == 'HIS':
    #         exclude_list_his = []
    #         for residue in exclude_list:
    #             st = residue.find('-')
    #             en = residue.find('_')
    #             resname = residue[0:st]
    #             resid = int(residue[st + 1: en])
    #             chain_seg = residue[en+1:]
    #             if resname in ['HSD', 'HSE', 'HIS', 'HSP']:
    #                 exclude_list_his.append("%s_%s" % (resid, chain_seg))
    #
    #         for atom in pdb.struct.iter_atoms():
    #             # HIS is considered as neutral Histidine, since the default behaviour is to make it neutral.
    #             if atom['resname'] in ['HSD', 'HSE', 'HIS']:
    #                 hist_id = "%s_%s" % (atom['resid'], atom['segname'])
    #                 if hist_id in exclude_list_his:
    #                     continue
    #
    #                 if (charged == True) or (charged is None):
    #                     atom['resname'] = 'HSP'
    #                     recreate_struct = True
    #                 elif charged == False:
    #                     # Don't change anything.
    #                     pass
    #             elif atom['resname'] == 'HSP':
    #                 hist_id = "%s_%s" % (atom['resid'], atom['segname'])
    #                 if hist_id in exclude_list_his:
    #                     continue
    #
    #                 if (charged == False) or (charged is None):
    #                     # HSD is the Karlsberg+ reference state, so it is choosen here as a default for neutral Histidine.
    #                     atom['resname'] = 'HSD'
    #                     recreate_struct = True
    #                 elif charged == True:
    #                     # Don't change anything.
    #                     pass
    #
    #
    #
    #     if recreate_struct:
    #         pdb.create_struct()

    # def convert_kbp_residues_epp_dpp(self):
    #     # Renames all EPP/DPP in ASP/GLU and renames the atoms HD1<->HD2 in Asp and HE1<->HE2 in Glu to make sure,
    #     # that the ASPP/GLUP patches create the (neutral) reference state of Karlsberg+.
    #
    #     pdb = self.structure
    #     for seg in pdb.struct.iter_segments():
    #         for res in seg.iter_residues():
    #             if res.resname == 'DPP':
    #                 new_resname = 'ASP'
    #                 atom_to_switch1 = 'HD1'
    #                 atom_to_switch2 = 'HD2'
    #             elif res.resname == 'EPP':
    #                 new_resname = 'GLU'
    #                 atom_to_switch1 = 'HE1'
    #                 atom_to_switch2 = 'HE2'
    #             else:
    #                 continue
    #
    #             atm1 = res[atom_to_switch1]
    #             atm2 = res[atom_to_switch2]
    #             atm1['name'] = atom_to_switch2
    #             atm2['name'] = atom_to_switch1
    #
    #             for atom in res.iter_atoms():
    #                 atom['resname'] = new_resname
    #
    #     pdb.create_struct()

    # def mutate_residue(self, resname, resid, chain_seg, new_resname):
    #     # Change reisdue names in the structure.
    #     pdb = self.structure
    #     seg = pdb.struct[chain_seg]
    #     residue = None
    #     for res in seg.iter_residues():
    #         if res.resid == resid and res.resname == resname:
    #             for atom in res.iter_atoms():
    #                 atom['resname'] = new_resname
    #             residue = res
    #
    #     pdb.create_struct()
    #
    #
    #     # Make sure the changed residues will be minimized.
    #     self.charmm_instructions['minimize'].append((residue, residue))
    #
    #     # Todo: Make sure, that this residues atoms are not added to the list for minimization when looking for missing atoms!

    def force_energy_min(self):
        # Hydrogens will be minimized, if this action has not already been triggered by "check_structures".
        self.charmm_instructions['do_minimize'] = True


def submit_charmm(workdir, input_filename, output_filename, wait_to_complete=False):

    """Function that submits a charmm job to the cluster via shell script."""

    charmms_file = workdir + 'charmm_run.sh'
    file = open(charmms_file, 'w')
    script = """#!/bin/tcsh

cd  %s """ % workdir + """
charmm_c39a2q < %s > %s """ %(input_filename, output_filename)+ """

        """

    script += '\n'
    file.write(script)
    file.close()
    st = os.stat(charmms_file)
    import stat
    os.chmod(charmms_file, st.st_mode | stat.S_IEXEC)

    # print('Submiting CHARMM job. The success of it must be determined by the user')
    shell = subprocess.Popen('csh\n',\
            stdin=subprocess.PIPE,\
            stdout=subprocess.PIPE,\
            stderr=subprocess.PIPE,\
            shell=True\
            )

    qsub_parameter = '-q D63.q'
    # qsub_parameter = '-q D61.q,D62.q,D63.q,D64.q'
    shell.stdin.write('qsub %s -o %s -e %s %scharmm_run.sh\n' % (qsub_parameter, workdir, workdir, workdir))
    shell.stdin.write('exit\n')
    print('Submited CHARMM job.')


    if wait_to_complete:
        terminate = False
        while not terminate:
            #todo: Implement better communication with the Grid engine, how to know the job outcome!
            if os.path.exists(workdir + output_filename):
                f_out = open(workdir + output_filename, 'r')
                for next_line in f_out:
                    if 'NORMAL TERMINATION BY NORMAL STOP' in next_line:
                        terminate = True
                    if 'ABNORMAL TERMINATION' in next_line:
                        terminate = True


def calc_gbsw_energy(structure, top, par, workdir='/tmp/gbsw_energy/', make_neutral=True, clean_up=True):
    """
        Uses the Charmm_manager object to calculate the energy of a protein with GBSW.

        Parameters:
        structure:  String with path to a pdb file or a Simple_struct_parser object.
        top:        List containing one or several path to the required topologie files.
        par:        List containing one or several path to the required parameter files.
        workdir:    A folder where the input and output files will be stored. if it exists, it must be empty.
                    Default: '/tmp/gbsw_energy'
        clean_up:   All created files and the workdir are deleted after the calculation has finished.
                    Default: True
    """

    if workdir[-1] != '/':
        workdir += '/'

    if not os.path.exists(workdir):
        os.mkdir(workdir)
    else:
        files_in_workdir = os.listdir(workdir)
        if files_in_workdir:
            error = "There are files in the workdir folder! The folder must be empty: " + workdir
            raise AssertionError(error)


    c = Charmm_manager(workdir=workdir, top=top, par=par)

    c.add_structure(structure)

    c.check_structures(quiet=True)

    if make_neutral:
        # The Residues DPP and EPP are a bit problematic.They are not recognized by gbsw. Also the reference state
        # for Karlsberg+ is created by a patch that places the hydrogen on the opposite side then the patch
        # for GLU/ASP (ASPP, GLUP).

        c.convert_kbp_residues_epp_dpp()





        # Choose the Karlsberg+ reference Histidine. Without this line HSP from Kalrsberg+ pqr would be interpreted as
        # charged HIS.
        c.change_type_prot('HIS')

        # c.change_type_prot('TYR')
        c.change_type_prot('ARG')
        c.change_type_prot('LYS')
        c.change_type_prot('ASP')
        c.change_type_prot('GLU')
        # c.change_type_prot('CYS')
        c.change_type_prot('CTE')
        c.change_type_prot('NTE')


    build_missing_index = c.charmm_config['tasks'].index('build_missing')
    c.charmm_config['tasks'].pop(build_missing_index)

    # Using the gbsw parameter file makes a little difference!
    # /scratch/scratch/tmeyer/CHARMM_NAMD/toppar_36/gbsw/par_all22_prot_gbsw.inp
    gbsw_commands = []
    gbsw_commands.append(" ")
    # gbsw_commands.append("ENERgy")
    gbsw_commands.append("stream \"/scratch/scratch/tmeyer/CHARMM_NAMD/toppar_36/gbsw/radius_gbsw.str\"")
    gbsw_commands.append("scalar wmain statistics select .not. type H* end")
    gbsw_commands.append("define check select (.not type H* ) .and. ( prop wmain .eq. 0.0 ) show end")
    gbsw_commands.append("if ?nsel ne 0  stop       !some heavy atom have a zero radius")
    gbsw_commands.append(" ")
    # gbsw_commands.append("ENERgy EPS 4.0")
    # CHARMM defaults due to: http://www.charmmtutorial.org/index.php/The_Energy_Function
    # ENERgy NBXMod  5 ATOM CDIEl SHIFt VATOm VDIStance VSWItch -
    #    CUTNb 14.0 CTOFnb 12.0 CTONnb 10.0 EPS 1.0 E14Fac 1.0 WMIN 1.5

    # Example 1 in the charmm docs
    #  !To perform a single-point energy calculation with infinite cutoffs:
    # gbsw_commands.append("GBSW sw 0.3 sgamma 0.03 dgp 1.5 GBenergy  epsp 4")
    # From mmtbs tutorial: http://mmtsb.org/workshops/mmtsb-ctbp_2006/Tutorials/GBSW_Tutorial/GBSW_Tutorial.html
    # gbsw_commands.append("gbsw sgamma 0.005 nang 50 dgp 1.5 GBenergy epsp 4")
    gbsw_commands.append("GBSW conc 0.1 sgamma 0.03 GBenergy epsp 4 MOLSURF")
    # removed since it is the default value:
    # dgp 1.5
    # sw 0.3
    gbsw_commands.append("ENERgy EPS 4.0")
    # gbsw_commands.append("ENERgy")
    gbsw_commands.append(" ")
    gbsw_commands.append(" ")
    insert_pos = c.charmm_config['tasks'].index('minimize_modeled')
    for cmd in reversed(gbsw_commands):
        c.charmm_config['tasks'].insert(insert_pos, cmd)
    # insert_pos = c.charmm_config['tasks'].index('write_structure')
    # c.charmm_config['tasks'].insert(insert_pos, 'ENERgy')


    ### From the CHARMM documentation ###
    # 2.  Choice of SW
    # In prinicple, one can choose any SW. However, it should be noted that
    # GBSW calculations take more time as SW increases. As default, SW=0.3
    # is recommended for the smooth boundary and SW=0.2 for the molecular
    # surface.
    #
    # The nonpolar solvation contribution is considered only when non-zero
    # SGAMMA is issued.  Note that the dimension is kcal/(molxA^2), and 0.01
    # to 0.04 might be suitable for SGAMMA.
    #
    # SW	[0.3]	 : half of smoothing length in Ang.
    #           (default value is changed to 0.2 when MOLSURF is issued.)
    # AA0	[aa0(sw)]: coefficient for the Coulomb Field Approximation term
    # AA1	[aa1(sw)]: coefficient for the correction term
    #           (optimized default values for aa0(sw) and aa1(sw)
    #           are given below)
    #
    # MOLSURF	[FALSE]  : approximation to PB with molecular surface
    # GBENER  [FALSE]  : calculate and print the solvation energy
    #           (No cutoff is used for GB electrostatic solvation energy.)
    # ROTINV  [FALSE]  : rotationally invariant numerical quadrature procedure
    #
    # NANG	[38]	 : number of angular integration points
    # NRAD 	[0]	 : number of radial integration points
    #                    (default value means the use of optimized 24 radial
    #            integration points)
    # RMAX	[20.0]	 : maximum distance for radial integration in Ang.
    # DGP	[1.5]	 : grid spacing for lookup table in Ang.
    # RBUFFER	[0.0]	 : buffer length for lookup table in Ang.
    #
    # EPSP	[1.0] 	 : dielectric constant of both protein and reference state
    # EPSW	[80.0]	 : solvent dielectric constant
    # CONC	[0.0] 	 : salt concentration in M
    # TEMP	[300.0]  : temperature in K (only necessary with CONC)
    # SGAMMA	[0.0]	 : nonpolar surface tension coefficients in kcal/(molxA^2)
    #
    # TMEMB	[0.0]	 : thickness of low-dielectric membrane slab centered
    #                    at Z=0 (in Ang.)
    # MSW	[sw]	 : half of membrane switching length in Ang.
    #
    # IGBFRQ	[1]	 : updating frequency of effective Born radii


    # c.force_energy_min()


    c.run_charmm()

    # f= open(workdir + "charmm.out", 'w')
    # f.writelines(c.charmm_output)
    # f.close()

      # Electrostatic solvation energy   =      -585.23072 [kcal/mol]
      #
      #
      # Surface area                     =      8421.69589 [Angs**2]
      # Nonpolar solvation energy        =        42.10848 [kcal/mol]
      #
      # Total solvation energy           =      -543.12224 [kcal/mol]

    # ENER ENR:  Eval#     ENERgy      Delta-E         GRMS
    # ENER INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals    IMPRopers
    # ENER CROSS:           CMAPs        PMF1D        PMF2D        PRIMO
    # ENER EXTERN:        VDWaals         ELEC       HBONds          ASP         USER
    # ENER PBEQ:             PBnp       PBelec        GBEnr         GSBP
    #  ----------       ---------    ---------    ---------    ---------    ---------
    # ENER>        0  -4247.89978   1857.21370     16.84889
    # ENER INTERN>      524.46169    596.14538     73.19139    583.84434     56.32630
    # ENER CROSS>      -198.49359      0.00000      0.00000      0.00000
    # ENER EXTERN>     -544.13800  -3482.02358      0.00000    256.96442      0.00000
    # ENER PBEQ>          0.00000      0.00000  -2114.17813      0.00000
    #  ----------       ---------    ---------    ---------    ---------    ---------

    energies_total = []
    energies_CMAPs = []
    energies_ELEC  = []
    energies_GBEnr = []
    energies_tot_solv = []
    energies_elec_solv = []
    energies_nonpolar_solv = []

    for line in c.charmm_output.split('\n'):
        entries = line.split()
        if entries and entries[0] == 'ENER>':
            # 4.184kJ/mol= 1.0kcal/mol.
            # Convert to kJ/mol
            e = float(entries[2]) * 4.184
            energies_total.append(e)
        if entries and entries[0] == 'ENER' and entries[1] == 'CROSS>':
            # Convert to kJ/mol
            e = float(entries[2]) * 4.184
            energies_CMAPs.append(e)
        if entries and entries[0] == 'ENER'and entries[1] == 'EXTERN>':
            # Convert to kJ/mol
            e = float(entries[3]) * 4.184
            energies_ELEC.append(e)
        if entries and entries[0] == 'ENER' and entries[1] == 'PBEQ>' :
            # Convert to kJ/mol
            e = float(entries[4]) * 4.184
            energies_GBEnr.append(e)

        reg = re.compile(r'^\s*Total solvation energy[\s=]+([-.\d]+) \[kcal/mol\]$')
        # reg = re.compile(r'^\s*Nonpolar solvation energy[\s=]+([-.\d]+) \[kcal/mol\]$')
        # reg = re.compile(r'^\s*Electrostatic solvation energy[\s=]+([-.\d]+) \[kcal/mol\]$')
        reg_m = reg.match(line)
        if reg_m is not None:
            e = 1.0 * float(reg_m.groups()[0]) * 4.184
            energies_tot_solv.append(e)

        reg = re.compile(r'^\s*Electrostatic solvation energy[\s=]+([-.\d]+) \[kcal/mol\]$')
        # reg = re.compile(r'^\s*Electrostatic solvation energy[\s=]+([-.\d]+) \[kcal/mol\]$')
        reg_m = reg.match(line)
        if reg_m is not None:
            e = 1.0 * float(reg_m.groups()[0]) * 4.184
            energies_elec_solv.append(e)

        reg = re.compile(r'^\s*Nonpolar solvation energy[\s=]+([-.\d]+) \[kcal/mol\]$')
        reg_m = reg.match(line)
        if reg_m is not None:
            e = 1.0 * float(reg_m.groups()[0]) * 4.184
            energies_nonpolar_solv.append(e)



    # for e in energies:
    #     print(e)
    if clean_up:
        shutil.rmtree(workdir)

    if not energies_CMAPs:
        energies_CMAPs.append(0)

    # print energies_ELEC
    # return (energies_total[0], energies_CMAPs[0], energies_ELEC[0], energies_GBEnr[0], energies_tot_solv[0])
    return (energies_elec_solv[0], energies_nonpolar_solv[0], energies_ELEC[0], energies_total[0])


if __name__ == "__main__":

    pdb        = "/scratch/scratch/tmeyer/md_pka/mod_bio/1a2p_m.pdb1"
    # workdir    = "/scratch/scratch/tmeyer/tmp/modelling/"
    workdir    = "/scratch/scratch/tmeyer/md_pka/mod_bio/1a2p_m_charmm/"
    charmm_bin = "charmm36"

    top = []
    # top.append("/scratch/scratch/tmeyer/CHARMM_NAMD/toppar_36/top_all27_prot_na.rtf")
    top.append("/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/top.inp")
    top.append("/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/patches.rtf")
    par = []
    # par.append("/scratch/scratch/tmeyer/CHARMM_NAMD/toppar_36/par_all27_prot_na.prm")
    par.append("/scratch/scratch/tmeyer/CHARMM_NAMD/toppar_36/par_all22_prot.inp")
    par.append("/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/patches.prm")


    charmm = Charmm_manager(workdir, pdb, charmm_bin, top, par)
    third_residue = charmm.structure.struct['A'][3]
    charmm.add_sequence(third_residue, 'ALA GLN', True)
    charmm.check_structures(quiet=False)
    # charmm.charmm_instructions['do_minimize'] = False
    charmm.run_charmm()

    # # Add the new residue to dict of titratable residues.
    # titr_residues = charmm.get_titr_residues()
    # template_state = {'charge' : 0,
    #                   'patch' : None,
    #                   'external_patches' : None,
    #                   'rename' : None,
    #                   'special' : None}
    # resname = 'EPP'
    # titr_residues[resname] = [dict(template_state) for i in range(3)]
    # titr_residues[resname][0]['charge'] = -1
    # titr_residues[resname][1]['charge'] = 0
    # titr_residues[resname][1]['patch'] = 'EPP2'
    # titr_residues[resname][2]['charge'] = 0
    # titr_residues[resname][2]['patch'] = 'EPP1'
    # charmm.set_titr_residues(titr_residues)
    #
    #
    # residue_to_rename = ('GLU', 7, 'A')
    # patch_translation_dict = {'GLUP' : 'EPP2'}
    #
    # # charmm.set_prot_residue(residue_to_rename, charge=0)
    # charmm.add_patch('GLUP', residue_to_rename)
    #
    # charmm.rename_residue(residue_to_rename, 'EPP', patch_translation_dict)
    # charmm.run_charmm()

    # modelled_structure = charmm.get_modelled_structure()


# def set_dolly_run_charmm(folder, dolly_number, charmm_structures):
#
#     tcsh_script_name = folder + "run%s.sh" % dolly_number
#     tcsh_script = """\n"""
#     for i, charmm_struct in enumerate(charmm_structures):
#
#          ### Generate shell script ###
#         tcsh_script += """\n"""
#         tcsh_script += """#!/bin/tcsh"""
#         tcsh_script += """
# set JOBNAME%i=%s """ % (i,charmm_struct.title)    + """
# source /scratch/scratch/jdragelj/mfes_settings.csh
# cd %s """ % charmm_struct.workdir + """
# %s < %s_charmm.inp > %s_charmm.out """ %  (charmm_struct.charmm_bin, charmm_struct.title, charmm_struct.title) + """
# mfes-0.3b.x86_64 -i config.in -c
# rm %sexclusion.stl """ % charmm_struct.workdir + """
# rm %sexclusion.vol """ % charmm_struct.workdir + """
# """
#         tcsh_script += """\n"""
#
#     f = open(tcsh_script_name, 'a')
#     f.write(tcsh_script)
#     f.close()
#
#     # Add right to execute file.
#     import stat
#     st = os.stat(tcsh_script_name)
#     os.chmod(tcsh_script_name, st.st_mode | stat.S_IEXEC)



def get_charmm_energy(workdir, inter=False, KJ=True):

    '''
    Gets a list of energies from output of CHARMM run

    Example output:

    # ENER ENR:  Eval#     ENERgy      Delta-E         GRMS
    # ENER INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals    IMPRopers
    # ENER CROSS:           CMAPs        PMF1D        PMF2D        PRIMO
    # ENER EXTERN:        VDWaals         ELEC       HBONds          ASP         USER
    #  ----------       ---------    ---------    ---------    ---------    ---------
    # ENER>        0   1522.52937      0.00000     18.03137
    # ENER INTERN>      544.02813   1124.93499    169.23668    713.96413     66.85994
    # ENER CROSS>      -176.10529      0.00000      0.00000      0.00000
    # ENER EXTERN>     -468.39449   -451.99472      0.00000      0.00000      0.00000
    #  ----------       ---------    ---------    ---------    ---------    ---------
    '''

    if workdir[-1] != '/':
        workdir += '/'

    charmm_energies = {}

    files = os.listdir(workdir)
    for file in files:
        if '_charmm.out' in file:
            charmm_out_filename = workdir + file


    f = open(charmm_out_filename)
    for line in f:
        ener_keyword = 'ENER ENR:'
        if inter:
            ener_keyword = 'INTE ENR:'
        if ener_keyword in line:
            types_and_names = []
            types_and_values = []
            while not '----------' in line:
                line = line.strip('\n')
                types_and_names.append(line.split(':'))
                line = f.next()
            line = f.next()
            while not '----------' in line:
                types_and_values.append(line.split('>'))
                line = f.next()
            break


    for (types_descr, names), (types, values) in zip(types_and_names, types_and_values):
        names_list = names.split()
        values_list = values.split()
        for name, value in zip(names_list, values_list):
            value = float(value)
            # Convert to kJ/mol
            if KJ:
                value *= 4.184
            if not name in charmm_energies:
                charmm_energies[name] = [value]
            else:
                charmm_energies[name].append(value)
        f.close()

    return charmm_energies


def get_charmm_sasa(workdir):
    '''
    :param workdir: Charmm otuput directory, function will find a file itself
    :return: Returns total sasa value only

    example:
 CHARMM>    coor surf select segid FLUR end acce rpro 1.4
 SELRPN>     42 atoms have been selected out of  37698
 SURFAC: Lennard-Jones radii values being used
 SURFAC: RPRObe=   1.40000
 SURFAC: Analytic surface area method used

 CHARMM>    scalar wmain stat
 Statistics for37698 selected atoms:
       minimum =    0.00000      maximum =    56.2405     weight =    37698.0
       average =   1.699390E-02  variance=   0.689882     total  =    640.636
    '''


    if workdir[-1] != '/':
        workdir += '/'

    charmm_sasa = None

    files = os.listdir(workdir)
    for file in files:
        if '_charmm.out' in file:
            charmm_out_filename = workdir + file

    f = open(charmm_out_filename)
    extract = False
    for line in f:
        ener_keyword = 'SCALAR WMAIN STAT'
        if ener_keyword in line:
            extract = True
        if extract:
            if 'average' in line:
                if 'variance' in line:
                    if 'total' in line:
                        line_data = line.split()
                        charmm_sasa = line_data[7]
    if charmm_sasa is None:
        print workdir
        raise AssertionError('SASA not extracted!')

    return charmm_sasa





    # energies_to_plot= ['ENERgy', 'ELEC', 'VDWaals', 'BONDs', 'ANGLes', 'DIHEdrals', 'IMPRopers']
    # from matplotlib import pyplot as plt
    # plt.figure()
    # for energy_to_plot in energies_to_plot:
    #     y = charmm_energies[energy_to_plot]
    #     plt.plot(y)
    # plt.legend(energies_to_plot)#
    # plt.show()

    # ENER ENR:  Eval#     ENERgy      Delta-E         GRMS
    # ENER INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals    IMPRopers
    # ENER CROSS:           CMAPs        PMF1D        PMF2D        PRIMO
    # ENER EXTERN:        VDWaals         ELEC       HBONds          ASP         USER
    #  ----------       ---------    ---------    ---------    ---------    ---------
    # ENER>        0   1522.52937      0.00000     18.03137
    # ENER INTERN>      544.02813   1124.93499    169.23668    713.96413     66.85994
    # ENER CROSS>      -176.10529      0.00000      0.00000      0.00000
    # ENER EXTERN>     -468.39449   -451.99472      0.00000      0.00000      0.00000
    #  ----------       ---------    ---------    ---------    ---------    ---------
