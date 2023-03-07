# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 10:54:41 2011

@author: tmeyer
"""

import numpy as np
import os
import re
from collections import defaultdict

import subprocess
import tempfile

# ToDo: first psf and then pdb is not checking for same number of atoms.
# ToDo: Detect models and warn or store them.
# ToDo: Warning for unresolved iCodes.



class Atom(object):
    def __init__(self, atom=None):
        self.data = {}
        if atom is None:
            self.data["hetatm"] = False
            self.data["name"] = None
            self.data["type"] = None
            self.data["resname"] = None
            self.data["segname"] = None
            self.data["resid"] = None
            self.data["icode"] = ' '
            self.data["coord"] = None
            self.data["index"] = None
            self.data["mass"] = None
            self.data["vdw"] = None
            self.data["charge"] = None
            self.data["chain"] = None
            self.data["occupancy"] = 1.0
            self.data["tempFactor"] = 0.0
            self.data["element"] = ' '
            self.data["altloc"] = ' '

            self.__list_position_in_host = None
        else:
            self.data = dict( atom.data )

    def __setitem__(self, index, value):
        self.data[index] = value
        return 1

    def __getitem__(self, index):
        return self.data[index]

    def keys(self):
        return self.data.keys()

    def has_key(self, key):
        return self.data.has_key(key)

    def __contains__(self, key):
        return key in self.data

    def copy(self):
        return Atom(self)



    # 1 -  6        Record name   "ATOM  "
    # 7 - 11        Integer       serial       Atom  serial number.
    #13 - 16        Atom          name         Atom name.
    #17             Character     altLoc       Alternate location indicator.
    #18 - 20        Residue name  resName      Residue name.
    #22             Character     chainID      Chain identifier.
    #23 - 26        Integer       resSeq       Residue sequence number.
    #27             AChar         iCode        Code for insertion of residues.
    #31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    #39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    #47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    #55 - 60        Real(6.2)     occupancy    Occupancy.
    #61 - 66        Real(6.2)     tempFactor   Temperature  factor.
    #77 - 78        LString(2)    element      Element symbol, right-justified.
    #79 - 80        LString(2)    charge       Charge  on the atom.
    # A PDB line written by Charmm:
    #ATOM      1  N   ASP     7      12.571  16.413  13.281  1.00 22.54      2OXP
    def write_pdb_line(self, set_charge=False):
        pdb_line = ""
        if not self['hetatm']:
            pdb_line +=  'ATOM'.ljust(6)
        else:
            pdb_line +=  'HETATM'.ljust(6)
        if self['index'] > 99999:
            pdb_line +=  '*****'
        else:
            pdb_line +=  str(self['index']).rjust(5)
        pdb_line +=  ' '
        name_offset = 2 - len(self['element'])
        # In case atom['element'] is not set.
        if name_offset == 2:
            name_offset = 1
        if len(self['name']) + name_offset > 4:
            name_offset = 0
        pdb_line +=  ''.rjust(name_offset)
        pdb_line +=  self['name'].ljust(4 - name_offset)
        pdb_line +=  self['altloc']
        pdb_line +=  self['resname'].ljust(4)
        #pdb_line +=  ' '
        #pdb_line +=  self['segname'][0].rjust(1)
        pdb_line +=  self['chain'].rjust(1)
        pdb_line +=  str(self['resid']).rjust(4)
        if self.has_key('icode'):
            pdb_line +=  self['icode'].rjust(1)
        else:
            pdb_line +=  ' '
        pdb_line +=  '   '
        pdb_line +=  '{0:.3f}'.format(self['coord'][0]).rjust(8)
        pdb_line +=  '{0:.3f}'.format(self['coord'][1]).rjust(8)
        pdb_line +=  '{0:.3f}'.format(self['coord'][2]).rjust(8)
        #pdb_line +=  ''.rjust(6)
        pdb_line +=  '{0:.2f}'.format(self['occupancy']).rjust(6)
        #pdb_line +=  ''.rjust(6)
        pdb_line +=  '{0:.2f}'.format(self['tempFactor']).rjust(6)

        # Fill the gap with segnames like Charmm is doing it.
        # pdb_line +=  ''.rjust(10)
        pdb_line +=  ''.rjust(6)
        pdb_line +=  self['segname'].ljust(4)

        #pdb_line +=  self['name'][0].rjust(2)
        pdb_line +=  self['element'].rjust(2)
        if not set_charge:
            #pdb_line +=  '{0:.1f}'.format(1.1).rjust(2)
            pdb_line +=  ''.rjust(2)
        else:
            #print self['charge']
            pdb_line +=  '{0:.1f}'.format( round(self['charge'], 1) ).rjust(2)
        return pdb_line

    def __str__(self):
        return self.write_pdb_line()


class SimpleAtomAccess(dict):

    """
     An super class for simple_atom_access_... classes.
    """

    def info(self):
        print "I'm:  " + str(type(self))

    def copy(self):
        # Todo: Macht das wirklich Sinn?
        # Find out what class I am.
        self_class = type(self)
        # Return copy of myself.
        return self_class(self)


class SimpleAtomAccessStructure(SimpleAtomAccess):

    """
    An object that allows an intuitive way to accessing chains/segments.
    """

    def __init__(self, *args, **kwargs):
        """
        Extends the base class function.
        """

        super(SimpleAtomAccessStructure, self).__init__(*args, **kwargs)

        if len(args) > 0 and isinstance(args[0], type(self)):
            self.chain_seg_dict = args[0].chain_seg_dict
        else:
            self.chain_seg_dict = {}

    def __getitem__(self, value):
        """
        Extends the base class function to support chain/segname interpretation.
        """
        # translate chain to segname

        if len(value) == 1:
            value = self.chain_seg_dict[value]

        return super(SimpleAtomAccessStructure, self).__getitem__(value)

    def add_chain_to_dict(self, chain, segname):
        self.chain_seg_dict[chain] = segname

    def has_key(self, key):
        # Look in segment names
        if super(SimpleAtomAccessStructure, self).has_key(key):
            return True
        # Look in chain names
        elif self.chain_seg_dict.has_key(key):
            return True
        else:
            return False

    def __contains__(self, key):
        return self.has_key(key)

    def keys(self):
        """
        Return list of segments.
        """
        return sorted( super(SimpleAtomAccessStructure, self).keys() )

    def segnames(self):
        """
        Same as keys()
        """
        return self.keys()

    def chainids(self):
        """
        Return list of chains.
        """
        return sorted( self.chain_seg_dict.keys() )

    def iter_segments(self):
        for segname in self.segnames():
            yield self[segname]

    def iter_chains(self):
        for chainid in self.chainids():
            yield self[chainid]

    def iter_residues(self):
        for seg in self.iter_segments():
            for res in seg.iter_residues():
                yield res

    def iter_atoms(self):
        for seg in self.iter_segments():
            for res in seg.iter_residues():
                for atm in res.iter_atoms():
                    yield atm


class SimpleAtomAccessSegment(SimpleAtomAccess):

    """
    An object that allows an intuitive way to accessing resids.
    """

    def __init__(self, *args, **kwargs):
        """
        Extends the base class function.
        """
        if len(args) > 0 and isinstance(args[0], type(self)):
            self.segname = args[0].segname
            self.chainid = args[0].chainid
        else:
            self.segname = None
            self.chainid = None
        super(SimpleAtomAccessSegment, self).__init__(*args, **kwargs)

    def __getitem__(self, key):
        """
        Extends the base class function.
        """
        if not isinstance(key, tuple):
            key = (key, ' ')
        return super(SimpleAtomAccessSegment, self).__getitem__(key)

    def __setitem__(self, key, value):
        """
        Extends the base class function.
        """
        if not isinstance(key, tuple):
            key = (key, ' ')
        super(SimpleAtomAccessSegment, self).__setitem__(key, value)

    def has_key(self, key):
        if not isinstance(key, tuple):
            key = (key, ' ')
        return super(SimpleAtomAccessSegment, self).has_key(key)

    def __contains__(self, key):
        return self.has_key(key)

    def keys(self):
        """
        Return list of (resid, icode) tuples.
        """
        return sorted( super(SimpleAtomAccessSegment, self).keys() )

    def resids(self):
        """
        Return list of resids. icodes are ignored, so a resid may occur more than once.
        """
        keys = super(SimpleAtomAccessSegment, self).keys()
        keys.sort()
        return [resid for resid, icode in keys]

    def get_segname(self):
        """
        Return the segname extracted from the first atom in the internal atom list.
        """
        return self.values()[0].values()[0]['segname']

    def get_chainid(self):
        """
        Return the chainid extracted from the first atom in the internal atom list.
        """
        return self.values()[0].values()[0]['chain']

    def get_first_atom(self):
        return self.iter_atoms().next()

    def iter_residues(self):
        for key in self.keys():
            yield self[key]

    def iter_atoms(self):
        for res in self.iter_residues():
            for atm in res.iter_atoms():
                yield atm


class SimpleAtomAccessResidue(SimpleAtomAccess):

    """
    An object that allows an intuitive way to accessing icode.
    """

    def __init__(self, *args, **kwargs):
        """
        Extends the base class function.
        """
        if len(args) > 0 and isinstance(args[0], type(self)):
            self.segname   = args[0].segname
            self.chainid   = args[0].chainid
            self.icode     = args[0].icode
            self.resid     = args[0].resid
            self.resname   = args[0].resname
            self.hetatm    = args[0].hetatm
        else:
            self.segname   = None
            self.chainid   = None
            self.icode     = None
            self.resid     = None
            self.resname   = None
            self.hetatm    = None
        super(SimpleAtomAccessResidue, self).__init__(*args, **kwargs)

    def has_key(self, key):
        if super(SimpleAtomAccessResidue, self).has_key(key):
            return True

    def __contains__(self, key):
        return self.has_key(key)

    def keys(self):
        """
        Overrides the base class function. Assumes that the atom name is consistent with the corresponding key.
        """
        atoms = sorted(self.values(), key=lambda a:a['index'])
        return [a['name'] for a in atoms]

    def atoms(self):
        """
        Same as keys()
        """
        return self.keys()

    def __getitem__(self, key):
        """
        Extends the base class function.
        """
        for atom in super(SimpleAtomAccessResidue, self).itervalues():
            if atom['name'] == key:
                return atom

        raise KeyError("Atom name '%s' is unknown in %s." % (key, str(self)))

    def __unicode__(self):
        return "%s-%i_%s" % (self.resname, self.resid, self.segname)
    def __str__(self):
        return self.__unicode__()
    def __repr__(self):
        return '<simple_atom_access_residue object for ' + self.__unicode__()  + '>'

    def __eq__(self, other):
        # if not isinstance(other, type(self)):
        #     print type(other)
        if other is None:
            return False
        atms_self  = [id(atm) for atm in  self.iter_atoms()]
        atms_other = [id(atm) for atm in other.iter_atoms()]
        return atms_self == atms_other

    def get_first_atom(self):
        return self.iter_atoms().next()

    def iter_atoms(self):
#                for atom in super(simple_atom_access_residue, self).itervalues():
#                    yield atom
        for atm in self.atoms():
            yield self[atm]

    def rename(self, resname):
        for atm in self.atoms():
            atom = self[atm]
            atom['resname'] = resname
        self.resname = resname

    def get_size_no_radii(self):
        min = list(self.get_first_atom()['coord'])
        max = list(self.get_first_atom()['coord'])
        for atom in self.iter_atoms():
            coord = atom['coord']
            for dim in range(3):
                if coord[dim] < min[dim]:
                    min[dim] = coord[dim]
                elif coord[dim] > max[dim]:
                    max[dim] = coord[dim]
        return (min, max)

    def get_size(self):
        min = list(self.get_first_atom()['coord'])
        max = list(self.get_first_atom()['coord'])
        for atom in self.iter_atoms():
            coord = atom['coord']
            vdw = atom['vdw']
            for dim in range(3):
                if (coord[dim] - vdw) < min[dim]:
                    min[dim] = coord[dim] - vdw
                elif (coord[dim] + vdw) > max[dim]:
                    max[dim] = coord[dim] + vdw
        return (min, max)


#noinspection PySimplifyBooleanCheck,PySimplifyBooleanCheck
class Simple_struct_parser(object):
    # version 0.1 (31.10.2011)

    def __init__(self):
        self.atoms = []
        self.psf_read = 0
        self.crd_read = 0
        self.par_read = 0
        
        self.pqr_read = 0

        self.radii_set = 0
        
        # For psf entries.
        self.vdw_table = {}

        # Futher information from par files.
        self.eps_table = {}

        # Information from top file.
        self.top_content = {}

        # For easy access to chains and resids.
        self.struct = {}

        # Stores all parsed files.
        self.file_log = []

        # Stores the PDB header. Only available if PDB has been read previously.
        # Structure: {<REMARK NUMBER> : [<COMMENT TEXT>]}
        self.header_dict = {}

        # Informations extracted from the PDB header. Only available if PDB has been read previously.
        self.header_info = {}
        # Information about missing entries obtained from REMARK 465.
        # Structure: [] with one entry for every missing segment: [(resname, chain, resid)]
        self.header_info['missing'] = []

        # Segments that had to be renumbered are stored here. The reason can be the use of iCodes or invalid resids.
        self.renames_segments = []

        # Connect entries at the end of a PDB file. List of touples. [(1,2), (1,3)]
        self.connections = []

        # A file parser should increase this number for every atom that is using a icode.
        self.icodes = 0
        

    def __getitem__(self, key):
        return self.struct[key]

    def keys(self):
        return self.struct.keys()
        
    def has_key(self, key):
        return self.struct.has_key(key)

    def __contains__(self, key):
        if key in self.struct:
            return True

        
    def copy(self, segname = None, residuelist = None, exclude = False):
        """
        Function to copy entries from a pdb

        Function will copy and return either everything (no arguments) or a list of segnames (segname = ...) or a list of residues (residuelist = ...)
        segnames or residues given can be either included (default, exclude = False) or excluded (set exclude=True) from copying

        args:
        [segname] = list of one or more segnames, e.g. ('ACHAIN', 'BCHAIN') or string with segname
        [residuelist] = array of residues with first giving the seg- or chain name and second the resid for each residue, e.g. [('ACHAIN', 1), ('ACHAIN', 3), ('BCHAIN', 7)]
        [exclude] = if exclude is set True, all segnames or residuelists are excluded, only rest is returned

        attention: [segname] and [residuelist] are exclusive, only one can be given, ([segname] has higher priority)
        """
        if type(segname) != list and segname is not None:
            segname = [segname]

        new = Simple_struct_parser()

        for atom_ref in self.atoms:
            atom = Atom(atom_ref)

            if segname is not None and not exclude:
                for sn in segname:
                    if atom['segname'] == sn:
                        new.new_atom( atom )
            elif segname is not None and exclude:
                keep=True
                for sn in segname:
                    if atom['segname'] == sn:
                        keep = False
                if keep == True:
                    new.new_atom( atom )

            elif residuelist is not None and not exclude:
                for residue in residuelist:
                    (segment, resid) = residue
                    if (atom['segname'] == segment or atom['chain'] == segment) and atom['resid'] == resid:
                        new.new_atom( atom )
            elif residuelist is not None and exclude:
                keep=True
                for residue in residuelist:
                    (segment, resid) = residue
                    if (atom['segname'] == segment or atom['chain'] == segment) and atom['resid'] == resid:
                        keep = False
                if keep == True:
                    new.new_atom( atom )
            else:
                new.new_atom( atom )

        new.psf_read  = self.psf_read
        new.crd_read  = self.crd_read
        new.par_read  = self.par_read
        new.pqr_read  = self.pqr_read

        new.vdw_table = dict(self.vdw_table)
        new.eps_table = dict(self.eps_table)

        new.top_content = {}
        if self.top_content:
            new.top_content['residues'] = {}
            for name, properties in self.top_content['residues'].iteritems():
                new.top_content['residues'][name] = {}
                new.top_content['residues'][name]['name'] = properties['name']
                new.top_content['residues'][name]['charge'] = properties['charge']
                if 'atoms' in properties.keys():
                    new.top_content['residues'][name]['atoms'] = {}
                    for atom_name, atom_properties in properties['atoms'].iteritems():
                        new.top_content['residues'][name]['atoms'][atom_name] = dict(atom_properties)
        else:
            new.top_content = self.top_content


        new.file_log = list(self.file_log)

        new.create_struct()

        return new

    def new_atom(self, atom):
        atom.__list_position_in_host = len(self.atoms)
        self.atoms.append(atom)
        return len(self.atoms)-1
        
#    def set_atom(self, index, atom):
#        #self.atoms[index] = 
#        pass

    def del_atom(self, atom, no_restruct=False):
#        for i, ref_atom in enumerate(self.atoms):
#            if ref_atom is atom:
#                self.atoms.pop(i)
#                return
        if type(atom) == list:
            del_indices = []
            for atm in atom:
                del_indices.append(atm.__list_position_in_host)

            for i in sorted(del_indices, reverse=True):
                self.atoms.pop(i)
            for i, atom in enumerate(self.atoms):
                atom.__list_position_in_host = i
        else:
            del_index = atom.__list_position_in_host
            self.atoms.pop(del_index)
            # Done in reassign_indices()
            # for i, atom in enumerate(self.atoms[del_index:]):
            #     atom.__list_position_in_host = i + del_index

        self.reassign_indices()
        if not no_restruct:
            self.create_struct()

        return 1

    def get_atom_ref(self, index):
        return self.atoms[index]

    def get_atom(self, index):
        return self.atoms[index].copy()
        
        
    def read_par(self, filename):
        self.file_log.append(filename)

        f = open(filename, 'r')
        
        for line in f:
            if line[0:9] == 'NONBONDED':
                break
        else:
             # print 'Error in "read_par": "NONBONDED" keyword not found!'
             return 0
            
        local_vdw_table = {}
        local_eps_table = {}
        sections = ['BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPER', 'CMAP']
        
        re_section = re.compile(r'(\w+).*')
        #!atom  ignored    epsilon      Rmin/2   ignored   eps,1-4       Rmin/2,1-4
        re_vdw_entry = re.compile(r'(\w+)\s+([\d\.-]+)\s+([\d\.-]+)\s+([\d\.-]+)')
        for line in f:
            mo = re_vdw_entry.match(line)
            if mo is not None:
                entries = mo.groups()
                atom_type  = entries[0]
                epsilon    = float(entries[2])
                vdw_radius = float(entries[3])
                
                if local_vdw_table.has_key(atom_type):
                    print 'Warning in "read_par": NONBONDED entry for '\
                            + atom_type + ' exists more than once in '\
                            + filename
                
                local_vdw_table[atom_type] = vdw_radius
                local_eps_table[atom_type] = epsilon

            
            # stop at the beginning of a new section
            mo = re_section.match(line)
            if mo is not None:
                if sections.count( mo.groups()[0] ) == 1:                    
                    break
                
        f.close()

        for atom_type in local_vdw_table.keys():
            vdw_radius = local_vdw_table[atom_type]
            epsilon = local_eps_table[atom_type]

            if self.vdw_table.has_key(atom_type):
                print 'Warning in "read_par": NONBONDED entry for '\
                        + atom_type + ' found in file ' + filename\
                        + ' exists and will be overwritten.'
            
            self.vdw_table[atom_type] = float(vdw_radius)
            self.eps_table[atom_type] = float(epsilon)

        self.par_read = 1

    def read_top(self, filename):
        """
        Parses the topology files.
        So far only names and charges are read.
        """

        if not self.top_content.has_key('residues'):
            self.top_content['residues'] = {}

        with open(filename, 'r') as f:
            last_residue = None
            for line in f:
#                    reg_m = reg.match(line)
                entries = line.split()
                if len(entries) == 0:
                    continue
                if entries[0] in ['RESI', 'PRES', 'END']:
                    if last_residue is not None:
                        # End of residue entry found. Store previous residue.
                        name = last_residue['name']
                        self.top_content['residues'][name] = last_residue
                        last_residue = None

                if entries[0] in ['RESI', 'PRES']:
                    # New residue found.
                    last_residue = {}
                    last_residue['name']   =           entries[1]
                    if (len(entries) > 2 and entries[2] != '!'):
                        last_residue['charge'] = int(float(entries[2]))
                elif entries[0] == 'ATOM':
                    # New ATOM entry found.
                    if last_residue is not None:
                        # The new ATOM entry belongs to a RESI entry. Adding it to the last found residue.
                        if not last_residue.has_key('atoms'):
                            last_residue['atoms'] = {}
                        name = entries[1]
                        last_residue['atoms'][name] = {}
                        last_residue['atoms'][name]['type']   =       entries[2]
                        last_residue['atoms'][name]['charge'] = float(entries[3])


    def read_xplor_psf(self, filename, allow_missmatch=False):
        # Todo: Take care of extended file format? Proabaly it does not matter.
        #  ATOMNO RESNO   RES  TYPE  X     Y     Z   SEGID RESID Weighting
        #    I5    I5  1X A4 1X A4 F10.5 F10.5 F10.5 1X A4 1X A4 F10.5
        #  ATOMNO RESNO   RES  TYPE  X     Y     Z   SEGID RESID Weighting
        #    I10   I10 2X A8 2X A8       3F20.10     2X A8 2X A8 F20.10

        self.file_log.append(filename)

        f = open(filename, 'r')
        
        num_of_atoms = 0
        for line in f:
            #   154896 !NATOM
            re_link = re.compile(r'\s+(\d+)\s+!NATOM')
            mo = re_link.match(line)
            if mo is not None:
                num_of_atoms = int( mo.groups()[0] )
                break
        if num_of_atoms == 0:
            print "Error in \"read_xplor_psf\": Could not determine the number of "\
                    + "atoms in psf file " + filename

        # xplor.psf                 
        #         9 A-PRO    1        GLN      HB2      HA       0.900000E-01   1.00800           1   0.00000     -0.301140E-02        
        # psf:
        #         1         1  GLN       N              79.5464806942       -7.7045483073       41.7368136925  A-PRO     1              48.3200000000
        re_link = re.compile(r"\s*(\d+)\s+([-_\w]+)\s+(\d+)"\
                + r"\s+(\w+)\s+(['\w]+)\s+(\w+)"\
                + r"\s+([\dE.-]+)\s+([\w.-]+)"\
                + r"\s+([\d.E-]+)\s+([\d.E-]+)\s+([\dE.-]+)\s+")
        atoms_found = 0        
        for line in f:
            mo = re_link.match(line)
            if mo is not None:
                entries = mo.groups()
                
                a = Atom()
                
                a["index"]   = int( entries[0] )
                a["segname"] = entries[1]
                a["resid"]   = int( entries[2] )
                a["resname"] = entries[3]
                a["name"]    = entries[4]
                a["type"]    = entries[5]
                a["charge"]  = float( entries[6] )
                a["mass"]    = float( entries[7] )
                a["coord"]   = None
                a["vdw"]     = None
                
                if self.crd_read == 0:
                    self.new_atom(a)
                else:
                    if allow_missmatch:
                        if self.struct.has_key( a["segname"] ):
                            if self.struct[ a["segname"] ].has_key( a["resid"] ):
                                a_stored = self.struct[ a["segname"] ][ a["resid"] ][ a["name"] ]
                            else:
                                continue
                        else:
                            continue
                    else:
                        a_stored = self.get_atom_ref(atoms_found)

                    # check that the correct atom is modified
                    for matching_key in ['index', 'resid', 'resname',
                                         'name', 'segname']:
                        if a[matching_key] != a_stored[matching_key]:
                            print 'Error in "read_xplor_psf": Entry "' + matching_key\
                                    + '" of atom ' + str(atoms_found)\
                                    + ' does not match with previous loaded coordinates ' \
                                    + 'content in ' + filename
                            f.close()
                            return 0

                    # add type, charge and mass
                    a_stored['type'] = a['type']
                    a_stored['charge'] = a['charge']
                    a_stored['mass'] = a['mass']
                    
                atoms_found += 1
            # if
        # for

        if num_of_atoms != atoms_found and not allow_missmatch:
            print 'Error in \"read_xplor_psf\": Could not find all expected atoms ('\
                    + str(atoms_found) + ' of ' + str(num_of_atoms) + ' in crd \
                    file ' + filename
            return 0
        
        self.psf_read = 1
        
        f.close()

        self.create_struct()
    
    
    def read_crd(self, filename):

        self.file_log.append(filename)

        f = open(filename, 'r')
        
        num_of_atoms = 0
        for line in f:
            #   154896 !NATOM
            # re_link = re.compile(r'\s+(\d+)\s+EXT')
            # mo = re_link.match(line)
            # if mo is not None:
            #     num_of_atoms = int( mo.groups()[0] )
            #     break
            if (len(line.split()) == 1 or len(line.split()) == 2) and line.split()[0] != '*':
                num_of_atoms = int( line.split()[0] )
                break

        if num_of_atoms == 0:
            print "Error in \"read_crd\": Could not determine the number of atoms in crd file " + filename


        #        17         1  GLN       HE22           75.6614477764       -8.3427226708       38.2630239971  A-PRO     1               0.0000000000
        # index, resid (continous), resname
        # name, x, y, z
        # segname, resid, ?
        re_link = re.compile(r"\s*(\d+)\s+(\d+)\s+(\w+)"\
                + r"\s+(['\w]+)\s+([\d\.-]+)\s+([\d\.-]+)\s+([\d\.-]+)"\
                + r"\s+([-_\w]+)\s+([\d]+)\s+([\d\.-E]+)\s+")
        atoms_found = 0        
        for line in f:
            mo = re_link.match(line)

            if mo is None:
                print '*****', line

            if mo is not None:
                entries = mo.groups()

                a = Atom()
                
                a["index"]   = int( entries[0] )
                a["resid"]   = int( entries[8] )
                a["resname"] = entries[2]
                a["name"]    = entries[3]
                x            = float( entries[4] )
                y            = float( entries[5] )
                z            = float( entries[6] )
                a["coord"]   = np.array( [x, y, z] , dtype=np.float)
                a["segname"] = entries[7]
                a["charge"]  = None
                a["mass"]    = None                    
                a["type"]    = None
                a["vdw"]     = None
                a["chain"]   = a["segname"][0]
                
                if self.psf_read == 0:
                    self.new_atom(a)
                else:
                    a_stored = self.get_atom_ref(atoms_found)
             
                    # check that the correct atom is modified
                    for matching_key in ['index', 'resid', 'resname',
                                         'name', 'segname']:
                        if a[matching_key] != a_stored[matching_key]:
                            print 'Error in "read_crd": Entry "' + matching_key\
                                    + '" of atom ' + str(atoms_found)\
                                    + ' does not match with previous loaded psf\
                                    content in ' + filename
                            f.close()
                            return 0
                    
                    # add coordiates
                    a_stored['coord'] = a['coord']
                    # add chain
                    a_stored['chain'] = a['chain']
                    

                atoms_found += 1
            # if
        # for
        
        f.close()

        # print num_of_atoms, atoms_found
        if num_of_atoms != atoms_found:
            print 'Error in "read_crd": Could not find all expected atoms ('\
                    + str(atoms_found) + ' of ' + str(num_of_atoms) + ' in crd \
                    file ' + filename
            
            return 0
        
        self.crd_read = 1

        self.create_struct()


    def read_pdb(self, filename, is_pqr=False):
        """Read a PDB file.

        is_pqr: If True the B-factor and occupancy entries are interpreted as charge and vdw radius.
        """

        self.file_log.append(filename)

        f = open(filename, 'r')
        
        atoms_found = 0

        segments_to_renumber = {}
        last_resid = defaultdict(int)
        resid_list = defaultdict(list)
        last_resid_entry = ''
        new_res = False

        header = []
        
        for line in f:
            if re.search(r"^REMARK", line):
                header.append(line)
                continue
            reg_m = re.search(r"(^ATOM)|(^HETATM)", line)
            if reg_m is not None:
            #if (line.find('ATOM') > -1) or (line.find('HETATM') > -1):
                if reg_m.groups()[0] is None:
                    is_hetatm = True
                else:
                    is_hetatm = False

                a = Atom()
                a['hetatm'] = is_hetatm

# 1 -  6        Record name   "ATOM  "
# 7 - 11        Integer       serial       Atom  serial number.
#13 - 16        Atom          name         Atom name.
#17             Character     altLoc       Alternate location indicator.
#18 - 20        Residue name  resName      Residue name.
#22             Character     chainID      Chain identifier.
#23 - 26        Integer       resSeq       Residue sequence number.
#27             AChar         iCode        Code for insertion of residues.
#31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
#39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
#47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
#55 - 60        Real(6.2)     occupancy    Occupancy.
#61 - 66        Real(6.2)     tempFactor   Temperature  factor.
#77 - 78        LString(2)    element      Element symbol, right-justified.
#79 - 80        LString(2)    charge       Charge  on the atom.

                entries = line

                if entries[6:11] != '*****':
                    a["index"]   = int( entries[6:11] )
                else:
                    a["index"]   = entries[6:11]
                a["name"]    = str( entries[12:16] ).strip()
                a["altloc"]  = str( entries[16] )
                a["resname"] = str( entries[17:21] ).strip()
                a["chain"]   = str( entries[21] )

                # Officially resid is only [22:26].
                icode_overwritten = False
                try:
                    int( entries[26] )
                    resid_str = entries[22:27]
                    icode_overwritten = True
                except:
                    resid_str  = entries[22:26]
                    a["icode"] = str( entries[26] )
                    if a["icode"].strip():
                        self.icodes += 1

                x            = float( entries[30:38] )
                y            = float( entries[38:46] )
                z            = float( entries[46:54] )
                a["coord"]   = np.array( [x, y, z], dtype=np.float )
                a["occupancy"]  = float( entries[54:60] )
                a["tempFactor"] = float( entries[60:66] )
                a["segname"]    = str( entries[72:76] ).strip(' \n')
                a["element"]    = str( entries[76:78].strip() )

                a["charge"]  = None
                a["mass"]    = None
                a["type"]    = None
                a["vdw"]     = None


                if a["chain"] == ' ' and not a["segname"]:
                    print 'Error in "re._crd": Nether Chain nor resname entry'\
                          + ' of atom ' + str(atoms_found)\
                          + ' exists in file ' + filename
                if not a["segname"]:
                    a["segname"] = a["chain"]
                if a["chain"] == ' ':
                    a["chain"] = a["segname"][0]

                # # print 'petra'
                # if a['segname'] in ['WC','WV','WP']:
                #     a['chain'] = a['segname'][1]
                # else:
                #     a['chain'] = a['segname'][0]

            #if icode_overwritten:
                #    if not segments_to_renumber.has_key(segname):
                #        print "Residues in segment %s will be renumbered, since icode was used." % segname
                #        segments_to_renumber[ segname ] = True

                # if a['resname'] == 'TIP3':
                #     segments_to_renumber[ segname ] = True

                segname = a["segname"]
                resid_ok = True
                try:
                    a["resid"] = int( resid_str )
                except ValueError:
                    a["resid"] = resid_str
                    resid_ok = False

                if resid_ok:
                    # Is this a new residue?
                    new_res = True
                    if (last_resid.has_key(segname)) and (last_resid[segname] == resid_str):
                        new_res = False

                    # If this is a new residue check if its resid has been used in this segment before.
                    if new_res:
                        if (resid_list.has_key(segname)) and (resid_str in resid_list[segname]):
                            resid_ok = False
                            # print segname
                            # print resid_str

                # Mark the segment,so it will be checked later. And set the entry for resid to a string, so it will be
                # changed later.
                if not resid_ok:
                    if not segments_to_renumber.has_key(segname):
                        print "Residues in segment %s will be renumbered, due to invalid resid: %s" % (segname, resid_str)
                    segments_to_renumber[ segname ] = True

                last_resid[segname] = resid_str
                resid_list[segname].append(resid_str)


                if not a["element"] or a["element"] == ' ':
                    # If the length of an atom name is 4 characters long the column that is normally reserved for
                    # 2-character elements has to be used. So it is more save to just consider the first letter as
                    # the element.
                    if len(a["name"]) < 4:
                        a["element"] = str( entries[12:14] ).strip()
                    else:
                        a["element"] = str( entries[12:13] ).strip()

                if is_pqr:
                    a["charge"]  = float( entries[54:60] )
                    a["vdw"]     = float( entries[60:66] )



                if self.psf_read == 0:
                    self.new_atom(a)
                else:
                    a_stored = self.get_atom_ref(atoms_found)
             
                    # check that the correct atom is modified
                    for matching_key in ['index', 'resid', 'resname',
                                         'name', 'segname']:
                        if a[matching_key] != a_stored[matching_key]:
                            print 'Error in "read_crd": Entry "' + matching_key\
                                    + '" of atom ' + str(atoms_found)\
                                    + ' does not match with previous loaded psf\
                                    content in ' + filename
                            f.close()
                            return 0
                    
                    # add coordiates
                    a_stored['coord'] = a['coord']
                    # add chain
                    a_stored['chain'] = a['chain']
                    # Todo: Add all pdb attributes like occupancy and B-Factor.
                
                atoms_found += 1
            # if
        # for
        
        f.close()
        
        self.crd_read = 1
        
        if is_pqr:
            self.pqr_read = 1



        # Repair chains with double or invalid resids.
        last_resid_entry = ''
        last_resid_set   = ''
        last_segname     = ''
        for atom in self.atoms:
            resid   = atom['resid']
            segname = atom['segname']
            # Only check segments stored in segments_to_renumber.
            if not segname in segments_to_renumber.keys():
                last_resid_entry = resid
                last_segname = segname
                continue

            new_res = True
            if (last_resid_entry == resid) and (last_segname == segname):
                new_res = False
            # If this is a new residue use the last resid increased by one. Or set it to 1 if it is a new segment.
            if new_res:
                if last_segname == segname:
                    new_resid = last_resid_set + 1
                else:
                    new_resid = 1

                atom['resid']  = new_resid
                last_resid_set = new_resid
            else:
                atom['resid'] = last_resid_set

            last_resid_entry = resid
            last_segname = segname



        self.create_struct()

        self.__analyse_pdb_header(header)

        self.renames_segments = segments_to_renumber

        return 1

    def __analyse_pdb_header(self, header):
#    REMARK 465
#    REMARK 465 MISSING RESIDUES
#    REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE
#    REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN
#    REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)
#    REMARK 465
#    REMARK 465   M RES C SSSEQI
#    REMARK 465     ALA A     1
#    REMARK 465     THR A     2
#    REMARK 465     SER A     3
#    REMARK 465     THR A     4
#    REMARK 465     LYS A     5
#    REMARK 465     LYS A     6
#    REMARK 465     LYS A    45
#    REMARK 465     HIS A    46
#    REMARK 465     PRO A    47
#    REMARK 465     LYS A    48
#    REMARK 465     LYS A    49
#    REMARK 465     GLY A    50
#    REMARK 465     GLU A   142
#    REMARK 465     ASP A   143
#    REMARK 465     ASN A   144
#    REMARK 465     ALA A   145
#    REMARK 465     ASP A   146
#    REMARK 465     SER A   147
#    REMARK 465     GLY A   148
#    REMARK 465     GLN A   149

        header_dict = defaultdict(list)

        reg = re.compile(r'^REMARK\s+(\d+)\s*(.*)$')
        for line in header:
            line = line.strip()
            reg_m = reg.match(line)
            if reg_m is not None:
                id      = reg_m.groups()[0]
                comment = reg_m.groups()[1]
                header_dict[id].append(comment)
        self.header_dict = header_dict

        # Check for missing sequences in REMARK 465
        missing_seg = []
        sequence = []
        reg = re.compile(r'\s*(\w+)\s+(\w)\s+(\d+)\s*')
        last_res_id = ('', '' , -999)
        for line in header_dict['465']:
            reg_m = reg.match(line)
            if reg_m is not None:
                resname = reg_m.groups()[0]
                chain   = reg_m.groups()[1]
                resid   = int(reg_m.groups()[2])
                res_id = (resname, chain, resid)

                if (last_res_id[1] == res_id[1]) and (last_res_id[2] == res_id[2] - 1):
                    # Add residue to last sequence.
                    sequence.append(res_id)
                    last_res_id = res_id
                else:
                    # End of sequence found. Start a new one.
                    if sequence:
                        missing_seg.append(sequence)
                    sequence = []
                    sequence.append(res_id)
                    last_res_id = res_id
        # Add final sequence.
        if sequence:
            missing_seg.append(sequence)

        self.header_info['missing'] = missing_seg

        return 1


    def read_pqr(self, filename):
        """
        Reads an APBS pqr file.
        """
        self.file_log.append(filename)

        f = open(filename, 'r')

        #APBS reads very loosely-formatted PQR files: all fields are whitespace-delimited rather than the strict
        #column formatting mandated by the PDB format. This more liberal formatting allows coordinates which are larger/smaller than ± 999 Å.
        #APBS reads data on a per-line basis from PQR files using the following format:
        #     0        1             2          3          X         4           5 6 7    8     9
        #Field_name Atom_number Atom_name Residue_name (Chain_ID) Residue_number X Y Z Charge Radius
        #ATOM   1172  HB1 LYS    84      22.091   9.803  -5.826  0.0000 0.0000
        atom_index = 0
        for line in f:
            entries = line.split()
            if entries[0] in ['ATOM', 'HETATM']:
                a = Atom()

                a["index"]   = int( entries[1] )
                a["name"]    = entries[2]
                a["resname"] = entries[3]

                resid_segname = entries[4]
                try:
                    a["resid"] = int(resid_segname)
                    shift = 0
                except ValueError:
                    a["segname"] = resid_segname
                    a["chain"]   = a["segname"][0]
                    a["resid"] = int(entries[5])
                    shift = 1

                x            = float( entries[5 + shift] )
                y            = float( entries[6 + shift] )
                z            = float( entries[7 + shift] )
                a["coord"]   = np.array( [x, y, z], dtype=np.float )

                a["charge"]  = float(entries[8 + shift])
                a["mass"]    = None
                a["type"]    = None
                a["vdw"]     = float(entries[9 + shift])


                if self.psf_read == 0:
                    self.new_atom(a)
                else:
                    a_stored = self.get_atom_ref(atom_index)

                    # check that the correct atom is modified
                    for matching_key in ['index', 'resid', 'resname',
                                         'name', 'segname']:
                        if a[matching_key] != a_stored[matching_key]:
                            print 'Error in "read_pqr": Entry "' + matching_key\
                                    + '" of atom ' + str(atom_index)\
                                    + ' does not match with previous loaded psf\
                                    content in ' + filename
                            f.close()
                            return 0

                    # add coordiates
                    a_stored['coord'] = a['coord']
                    # add chain
                    a_stored['chain'] = a['chain']
                    # add charge
                    a_stored['charge'] = a['charge']
                    # add vdw radius
                    a_stored['vdw'] = a['vdw']

                atom_index += 1

            # if
        # for

        f.close()

        self.pqr_read = 1

        self.create_struct()
    
    def write_pqr(self, filename, kb_style=False, no_water=False):
        """
            Writes a pqr file. if filename is 'None' no file is written but the content of the file
            is instead returned as string.

            Parameters:
            fileneme    String containing the filename or 'None'
            kb_style    If True the file is written in the PDB format as far sas possible.
                        If False the APBS definition of a PQR file is used..
        """

        # check if all required files have been read previously
        if not self.pqr_read == 1:
            # if self.psf_read == 0:
            #     print 'Error in "write_pqr": Atom charges are missing, please read'\
            #             + ' psf file first.'
            #     return -1
            if self.crd_read == 0:
                msg =  'Error in "write_pqr": Coordinates are missing, please read'\
                        + ' crd file first.'
                raise AssertionError(msg)
            if self.par_read == 0:
                msg = 'Error in "write_pqr": Atom radii are missing, please read'\
                        + ' parameter file first.'
                raise AssertionError(msg)
        
        # check that all charges and radii are defined
        for a in self.atoms:
            if a['charge'] is None:
                if self.top_content['residues'][a['resname']]['atoms'].has_key(a['name']):
                    charge = self.top_content['residues'][a['resname']]['atoms'][a['name']]['charge']
                else:
                    if self.top_content['residues']['NTER']['atoms'].has_key(a['name']):
                        charge = self.top_content['residues']['NTER']['atoms'][a['name']]['charge']
                    elif self.top_content['residues']['CTER']['atoms'].has_key(a['name']):
                        charge = self.top_content['residues']['CTER']['atoms'][a['name']]['charge']
                    else:
                        msg = 'Error in "write_pqr": Charge of atom could not be determined: '\
                            + str( a['index'] ) + ':' + str( a['name'] )
                        raise AssertionError(msg)

                a['charge'] = charge


            if a['charge'] is None:
                msg = 'Error in "write_pqr": Charge entry is missing in atom '\
                        + str( a['index'] ) + ':' + str( a['name'] )
                raise AssertionError(msg)

            # ToDo: Did not work with 2lzt
            # if a['type'] is None:
            #     if self.top_content['residues'][a['resname']]['atoms'].has_key(a['name']):
            #         atom_type = self.top_content['residues'][a['resname']]['atoms'][a['name']]['type']
            #     else:
            #         if self.top_content['residues']['NTER']['atoms'].has_key(a['name']):
            #             atom_type = self.top_content['residues']['NTER']['atoms'][a['name']]['type']
            #         elif self.top_content['residues']['CTER']['atoms'].has_key(a['name']):
            #             atom_type = self.top_content['residues']['CTER']['atoms'][a['name']]['type']
            #         else:
            #             print 'Error in "write_pqr": Charge of atom could not be determined: '\
            #                 + str( a['index'] ) + ':' + str( a['name'] )
            #
            #     a['type'] = atom_type


            if a['vdw'] is None and self.vdw_table[a['type']] is None:
                msg = 'Error in "write_pqr": Radius entry could bot be \
                        determined from parameter file in atom '\
                        + str( a['index'] ) + ':' + str( a['name'] )
                raise AssertionError(msg)


        lines = []

        if kb_style:
            # Karlsberg+ pqr are PDBs with radii and charge in the B-Factor and occupancy column.
            for atom in self.atoms:
                if no_water and atom['resname'] in ['HOH', 'TIP3']:
                    continue

                if atom['vdw'] is None:
                    vdw = self.vdw_table[atom['type']]
                else:
                    vdw = atom['vdw']

                new_atom = Atom(atom)
                new_atom['occupancy'] = new_atom['charge']
                new_atom['tempFactor'] = vdw
                lines.append(str(new_atom) + '\n')

    #             if not atom['hetatm']:
    #                lines.append( 'ATOM'.ljust(6) )
    #             else:
    #                 lines.append( 'HETATM'.ljust(6) )
    #             lines.append( str(atom['index']).rjust(5) )
    #             name_offset = 2 - len(atom['element'])
    #             # In case atom['element'] is not set.
    #             if name_offset == 2:
    #                 name_offset = 1
    #
    #             if len(atom['name']) + name_offset > 5:
    #                 name_offset = 0
    #             lines.append( ''.rjust(name_offset) )
    #             lines.append( atom['name'].ljust(5 - name_offset) )
    #
    # #            lines.append( atom['name'].rjust(5) )
    #             lines.append( ' ' )
    #             lines.append( atom['resname'].ljust(4) )
    #             lines.append( atom['segname'][0].rjust(1) )
    #             lines.append( str(atom['resid']).rjust(4) )
    #
    #             lines.append( '{0:.3f}'.format(atom['coord'][0]).rjust(12) )
    #             lines.append( '{0:.3f}'.format(atom['coord'][1]).rjust(8) )
    #             lines.append( '{0:.3f}'.format(atom['coord'][2]).rjust(8) )
    #
    # #            if kb_style:
    #             lines.append( '{0:.3f}'.format(atom['charge']).rjust(6) )
    #
    #             if atom['vdw'] is None:
    #                 lines.append( '{0:.3f}'.format( self.vdw_table[atom['type']] ).rjust(6) )
    #             else:
    #                 lines.append( '{0:.3f}'.format( atom['vdw'] ).rjust(6) )
    # #            else:
    # #                lines.append( '{0:.4f}'.format(atom['charge']).rjust(8) )
    # #
    # #                if atom['vdw'] is None:
    # #                    lines.append( '{0:.4f}'.format( self.vdw_table[atom['type']] ).rjust(7) )
    # #                else:
    # #                    lines.append( '{0:.4f}'.format( atom['vdw'] ).rjust(7) )
    #             lines.append('\n')
        else:
            # Use APBS style.
            #APBS reads very loosely-formatted PQR files: all fields are whitespace-delimited rather than the strict
            #column formatting mandated by the PDB format. This more liberal formatting allows coordinates which are larger/smaller than ± 999 Å.
            #APBS reads data on a per-line basis from PQR files using the following format:
            #Field_name Atom_number Atom_name Residue_name Chain_ID Residue_number X Y Z Charge Radius
            #ATOM   1172  HB1 LYS    84      22.091   9.803  -5.826  0.0000 0.0000
            lines.append("CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1\n")
            for atom in self.atoms:
                if no_water and atom['resname'] in ['HOH', 'TIP3']:
                    continue

                if not atom['hetatm']:
                   lines.append( 'ATOM'.ljust(6) )
                else:
                    lines.append( 'HETATM'.ljust(6) )
                if atom['index'] > 99999:
                    lines.append( "*****" )
                else:
                    lines.append( str(atom['index']).rjust(5) )

                lines.append(' ')

                name_offset = 2 - len(atom['element'])
                # In case atom['element'] is not set.
                if name_offset == 2:
                    name_offset = 1
                if len(atom['name']) + name_offset > 4:
                    name_offset = 0
                # lines.append( ''.rjust(name_offset) )
                lines.append( ' ' * name_offset )
                lines.append( atom['name'].ljust(4 - name_offset) )

    #            lines.append( atom['name'].rjust(5) )
                lines.append( ' ' )
                lines.append( atom['resname'].ljust(4) )

                # Chain IDs or segment names seem to be incomatible with at least the sasa binary.
                # lines.append( ' ' )
                # if atom['segname'] not in ['', ' ']:
                #    lines.append( atom['segname'] )
                # else:
                #    lines.append( atom['chain'] )

                lines.append( ' ' )
                l = len(str( atom['resid'] )) + 1
                if l <= 4 :
                    lines.append( str(atom['resid']).rjust(4) )
                else:
                    lines.append( str(atom['resid']).rjust( l ) )

                lines.append( ' '*3 )

                for i in [0,1,2]:
                    coor_str = '{0:.3f}'.format(atom['coord'][i])
                    l = len(coor_str) + 1
                    if l <= 9:
                        lines.append( coor_str.rjust(9) )
                    else:
                        lines.append( coor_str.rjust(l) )

                charge_str = '{0:.3f}'.format(atom['charge'])
                l = len( charge_str ) + 1
                if l <= 6:
                    lines.append( charge_str.rjust(6) )
                else:
                    lines.append( charge_str.rjust(l) )

                if atom['vdw'] is None:
                    vdw = self.vdw_table[atom['type']]
                else:
                    vdw = atom['vdw']
                vdw = '{0:.3f}'.format( vdw )
                l = len( vdw ) + 1
                if l <= 6:
                    lines.append( vdw.rjust(6) )
                else:
                    lines.append( vdw.rjust(l) )

                lines.append('\n')

        lines.append('END')

        if filename is not None:
            f = open(filename, 'w')
            for line in lines:
                f.write(line)
            f.close()
        else:
            s = ""
            for line in lines:
                s += line
            return s
        
        

    def write_pdb(self, filename, set_charge=False):
        # check if all required files have been read previously
        if self.crd_read == 0:
            print 'Error in "write_pdb": Coordinates are missing, please read crd or pdb file first.'
            return -1
        
        f = open(filename, 'w')

        if set_charge:
            if self.par_read == 0:
                print 'Error in "write_pdb": No psf was read. Charges will not be assigned.'
                set_charge = False
                
            # check that all charges and radii are defined
            for a in self.atoms:
                if a['charge'] is None:
                    print 'Error in "write_pdb": Charge entry is missing in atom '\
                            + str( a['index'] ) + ':' + str( a['name'] )
                    print '                      Charges will not be assigned.'
                    break

        
        #Field_name Atom_number Atom_name Residue_name Chain_ID Residue_number X Y Z Charge Radius
        #ATOM   1172  HB1 LYS    84      22.091   9.803  -5.826  0.0000 0.0000
        for atom in self.atoms:
#            if not atom['hetatm']:
#                f.write( 'ATOM'.ljust(6) )
#            else:
#                f.write( 'HETATM'.ljust(6) )
#            f.write( str(atom['index']).rjust(5) )
#            f.write( ' ' )
#            f.write( atom['name'].rjust(4) )
#            f.write( ' ' )
#            f.write( atom['resname'].ljust(3) )
#            f.write( ' ' )
#            #f.write( atom['segname'][0].rjust(1) )
#            f.write( atom['chain'].rjust(1) )
#            f.write( str(atom['resid']).rjust(4) )
#            if atom.has_key('icode'):
#                f.write( atom['icode'].rjust(1) )
#            else:
#                f.write( ' ' )
#            f.write( '   ' )
#            f.write( '{0:.3f}'.format(atom['coord'][0]).rjust(8) )
#            f.write( '{0:.3f}'.format(atom['coord'][1]).rjust(8) )
#            f.write( '{0:.3f}'.format(atom['coord'][2]).rjust(8) )
#
#            f.write( ''.rjust(6) )
#            #f.write( '{0:.2f}'.format(123.45).rjust(6) )
#            f.write( ''.rjust(6) )
#            #f.write( '{0:.2f}'.format(123.45).rjust(6) )
#            f.write( ''.rjust(10) )
#
#            f.write( atom['name'][0].rjust(2) )
#            if not set_charge:
#                #f.write( '{0:.1f}'.format(1.1).rjust(2) )
#                f.write( ''.rjust(2) )
#            else:
#                #print atom['charge']
#                f.write( '{0:.1f}'.format( round(atom['charge'], 1) ).rjust(2) )
            f.write( atom.write_pdb_line(set_charge) )
            f.write('\n')

        for connection in self.connections:
            line = "CONECT"
            for partner in connection:
                line += str(partner).rjust(5)
            line += '\n'
            f.write(line)
        f.write('END\n')
        f.close()

    def write_crd(self, filename, segname=[]):
        if type(segname) == str:
            segname = [segname]

        # Todo: use extended!

        # check if all required files have been read previously
        if self.crd_read == 0:
            print 'Error in "write_crd": Coordinates are missing, please read crd or pdb file first.'
            return -1

        #* Title
        # 1614
        #    1    1 ASP  N     12.57100  16.41300  13.28100 2OXP 7     22.54000
        #        17         1  GLN       HE22           75.6614477764       -8.3427226708       38.2630239971  A-PRO     1               0.0000000000
        # index, resid (continous), resname
        # name, x, y, z
        # segname, resid, ?
        # From http://scv.bu.edu/documentation/software-help/scientific-engineering/quantadocs/basic_ops/C_ff.html
        #1.   TITLE lines (character*80) begin with a *. Last title line is * followed by at least seven blanks. Natom --- defined as i5
        #2.   ATOM lines: atom# resid1 resnam atnam X Y Z segid resid2 bvalue
        #3.   format: I5,I5,1x,a4,1x,a4,3f10.5,1x,a4,1x,a4,f10.5
        #
        text = ''
        index = 0
        resid_continous = 0
        last_resid = None
        last_segname = None
        for atom in self.atoms:
            if segname:
                if not atom['segname'] in segname:
                    continue

            index += 1
            if atom['resid'] != last_resid or atom['segname'] != last_segname:
                resid_continous += 1
                last_resid = atom['resid']
                last_segname = atom['segname']

# #            f.write( str(atom['index']).rjust(7) )
#             text += ( str(index).rjust(5) )
#             text += ( str(resid_continous).rjust(5) )
#             text += (' ')
#             text += ( atom['resname'].ljust(4) )
#             text += (' ')
#             text += ( atom['name'].ljust(4) )
#             text += ( '{0:.5f}'.format(atom['coord'][0]).rjust(10) )
#             text += ( '{0:.5f}'.format(atom['coord'][1]).rjust(10) )
#             text += ( '{0:.5f}'.format(atom['coord'][2]).rjust(10) )
#             text += (' ')
#             text += ( atom['segname'][:4].ljust(4) )
#             text += (' ')
#             text += ( str(atom['resid']).ljust(4) )
#             text += ( '{0:.5f}'.format( 0.0 ).rjust(10) )
#             text += ('\n')

            # Use the extended file format as default.
            #      21119       868  POPC      H11S           49.6268386841        3.0000092983       -4.0521950722  MEMBRANE  91              0.0000000000
            # (2I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A8,F20.10)
            text += ( str(index).rjust(10) )
            text += ( str(resid_continous).rjust(10) )
            text += ('  ')
            text += ( atom['resname'].ljust(8) )
            text += ('  ')
            text += ( atom['name'].ljust(8) )
            text += ( '{0:.5f}'.format(atom['coord'][0]).rjust(20) )
            text += ( '{0:.5f}'.format(atom['coord'][1]).rjust(20) )
            text += ( '{0:.5f}'.format(atom['coord'][2]).rjust(20) )
            text += ('  ')
            text += ( atom['segname'][:4].ljust(8) )
            text += ('  ')
            text += ( str(atom['resid']).ljust(8) )
            text += ( '{0:.5f}'.format( 0.0 ).rjust(20) )
            text += ('\n')

        # text += ('\n')

        # Add number of atoms in first non-comment line
        # text = str( index ).rjust(5) + '\n' + text
        # For extendend file format.
        text = str( index ).rjust(10) + '  EXT\n' + text

        if filename is not None:
            f = open(filename, 'w')
            f.write(text)
            f.close()

        return text

    def write_apbs_par(self, filename):
        # check if all required files have been read previously
        if not self.pqr_read == 1:
            # if self.psf_read == 0:
            #     print 'Error in "write_pqr": Atom charges are missing, please read'\
            #             + ' psf file first.'
            #     return -1
            if self.crd_read == 0:
                print 'Error in "write_pqr": Coordinates are missing, please read'\
                        + ' crd file first.'
                return -1
            # if self.par_read == 0:
            #     print 'Error in "write_pqr": Atom radii are missing, please read'\
            #             + 'parameter file first.'
            #     return -1

        # check that all charges and radii are defined and set 'type field if necessary.'
        for a in self.atoms:
            if a['type'] is None:
                if self.top_content['residues'][a['resname']]['atoms'].has_key(a['name']):
                    type = self.top_content['residues'][a['resname']]['atoms'][a['name']]['type']
                else:
                    if self.top_content['residues']['NTER']['atoms'].has_key(a['name']):
                        type = self.top_content['residues']['NTER']['atoms'][a['name']]['type']
                    elif self.top_content['residues']['CTER']['atoms'].has_key(a['name']):
                        type = self.top_content['residues']['CTER']['atoms'][a['name']]['type']
                    else:
                        print 'Error in "write_pqr": Type of atom could not be determined: '\
                            + str( a['index'] ) + ':' + str( a['name'] )

                a['type'] = type

            if a['charge'] is None:
                if self.top_content['residues'][a['resname']]['atoms'].has_key(a['name']):
                    charge = self.top_content['residues'][a['resname']]['atoms'][a['name']]['charge']
                else:
                    if self.top_content['residues']['NTER']['atoms'].has_key(a['name']):
                        charge = self.top_content['residues']['NTER']['atoms'][a['name']]['charge']
                    elif self.top_content['residues']['CTER']['atoms'].has_key(a['name']):
                        charge = self.top_content['residues']['CTER']['atoms'][a['name']]['charge']
                    else:
                        print 'Error in "write_pqr": Charge of atom could not be determined: '\
                            + str( a['index'] ) + ':' + str( a['name'] )

                a['charge'] = charge



            if a['charge'] is None:
                print 'Error in "write_pqr": Charge entry is missing in atom '\
                        + str( a['index'] ) + ':' + str( a['name'] )

            if a['vdw'] is None and self.vdw_table[a['type']] is None:
                print 'Error in "write_pqr": Radius entry could bot be \
                        determined from parameter file in atom '\
                        + str( a['index'] ) + ':' + str( a['name'] )

            if self.eps_table[a['type']] is None:
                print 'Error in "write_pqr": Epsilon entry could bot be \
                        determined from parameter file in atom '\
                        + str( a['index'] ) + ':' + str( a['name'] )



        f = open(filename, 'w')
        f.write("WAT	OW	0.000000	1.7683	0.6364\n")

        # Format:
        #  <RESIDUE> <ATOM> <CHARGE> <RADIUS> <EPSILON>
        # where <RESIDUE> is the residue name, <ATOM> is the atom name, <CHARGE>
        # is the atomic charge in e, <RADIUS> is the van der Waals radius in
        # A, and <EPSILON> is the van der Waals well depth in kJ/mol.
        ##########################################################################
        #ALK	C	0.597200	1.9080	0.4577

        # Collect atom types and their properties.
        # Structure: atom_types['resname']['atom name'][
        #                                              'charge'
        #                                              'radius'
        #                                              'epsilon'
        #                                               ]
        atom_types = {}
        for atom in self.atoms:
            resname = atom['resname']
            name    = atom['name']
            if not atom_types.has_key(resname):
                atom_types[resname] = {}
            if not atom_types[resname].has_key(name):
                atom_types[resname][name] = {}
                atom_types[resname][name]['charge'] = atom['charge']

                if atom['vdw'] is None:
                    vdw = self.vdw_table[atom['type']]
                else:
                    vdw = atom['vdw']
                atom_types[resname][name]['radius'] = vdw
                atom_types[resname][name]['epsilon'] = self.eps_table[atom['type']]

        for resname in atom_types.keys():
            for name in atom_types[resname].keys():
                line = "%s   %5s   %7s   %7s   %7s" % (resname, name,\
                                                     "%2.3f" % atom_types[resname][name]['charge'],\
                                                     "%2.3f" % atom_types[resname][name]['radius'],\
                                                     "%2.3f" % atom_types[resname][name]['epsilon'],\
                                                        )
                f.write(line + '\n')

        f.close()

    def shift_segment(self, segname, value):
        """
        Increase the resid of all residues in segment 'segname' by 'value'
        """
        #WARUNGEN FALLS icodes existieren!!
        for atm in self.struct[segname].iter_atoms():
            atm['resid'] += value

    def shift_chain(self, chainid, value):
        """
        Increase the resid of all resdiues in chain 'chainid' by 'value'
        """
        self.shift_segment(chainid, value)

#        for atom in self.atoms:
#            if atom['chain'] == chain_id:
#                atom['resid'] += value

    def reassign_indices(self):
        for i, atm in enumerate(self.atoms):
            atm['index'] = i + 1
            atm.__list_position_in_host = i



    def create_struct(self):
        """
        Create and fill the "simple_atom_access" object "struct" of this structure.
        """

        self.struct = SimpleAtomAccessStructure()
        
        for atom in self.atoms:
            # Segment / Chain
            segname = atom['segname']
            chain = atom['chain']
            if not self.struct.has_key(segname):
                self.struct[segname] = SimpleAtomAccessSegment()
            # This command is only allowed the corresponding segname has been added to struct.
            self.struct.add_chain_to_dict(chain, segname)

            # Residue
            resid = atom['resid']
            if atom.has_key('icode'):
                icode = atom['icode']
            else:
                icode = ' '
            if not self.struct[segname].has_key( (resid, icode) ):
                self.struct[segname][resid, icode] = SimpleAtomAccessResidue()

            # Atom
#            if not self.struct[segname][resid].has_key(icode):
#                self.struct[segname][resid][icode] = {}
            name = atom['name']
            self.struct[segname][resid, icode][name] = atom

        for seg in self.struct.iter_segments():
            first_atom = seg.iter_atoms().next()
            seg.segname = first_atom['segname']
            seg.chainid   = first_atom['chain']
            for res in seg.iter_residues():
                first_atom = res.iter_atoms().next()
                res.segname   = first_atom['segname']
                res.chainid   = first_atom['chain']
                res.icode     = first_atom['icode']
                res.resid     = first_atom['resid']
                res.resname   = first_atom['resname']
                res.hetatm    = first_atom['hetatm']

                if 'sasa' in first_atom:
                    backbone_atoms = ['N', 'HN', 'CA', 'HA', 'C', 'O',
                                      'HT1', 'HT2', 'HT3',
                                      'OT1', 'OT2',
                                      'CT', 'NT', 'CAT', 'HNT']
                    sasa = 0
                    side_chain_sasa = 0
                    for atm in res.iter_atoms():
                        sasa += atm['sasa']

                        if atm['name'] not in backbone_atoms:
                            side_chain_sasa += atm['sasa']

                    res.sasa = sasa
                    res.side_chain_sasa = side_chain_sasa
                else:
                    res.sasa = None
                    res.side_chain_sasa = None

        return 1

    def resolve_icodes(self):
        # self.create_struct()
        
        # store the original resid
        for atom in self.atoms:
            atom['resid_org'] = atom['resid']
            atom['icode_org'] = atom['icode']

        restart = True
        while restart:
            for atom in self.struct.iter_atoms():
                if atom.has_key('icode') and atom['icode'] != ' ':
                    resid_ref   = atom['resid']
                    segname_ref = atom['segname']
                    icode_ref   = atom['icode']

                    # Increase resid of this residue and all following residues in
                    # the same chain by one.
                    for res in self.struct[segname_ref].iter_residues():
                        if (res.resid == resid_ref and res.icode >= icode_ref) or res.resid > resid_ref:
                            for atm in res.iter_atoms():
                                atm['resid'] += 1

                    # Remove icode entry. Old resids are still valid.
                    for atm in self[segname_ref][resid_ref, icode_ref].iter_atoms():
                        atm['icode'] = ' '

                    # recreate struct dictionary.
                    self.create_struct()
                    restart = True
                    break
            else:
                restart = False

#    def resolve_icodes(self):
#        self.create_struct()
#
#        # store the original resid
#        for atom in self.atoms:
#            atom['resid_old'] = atom['resid']
#            atom['icode_old'] = atom['icode']
#
#        icode_removed_dict           = {}
#        increase_resid_in_chain_dict = {}
#
#        for atom in self.atoms:
#            if atom.has_key('icode') and atom['icode'] != ' ':
#                resid   = atom['resid']
#                chainid = atom['chain']
#                icode   = atom['icode']
#
#                # Increase resid of this residue and all following residues in
#                # the same chain by one.
#                for (rid, residue) in self[chainid].iteritems():
#                    if rid < resid:
#                        continue
#
#                    for (ic, ic_a) in residue.iteritems():
#                        if rid == resid and ic < icode:
#                            continue
#
#                        for a in ic_a.itervalues():
#                            a['resid'] += 1
#
#                # Remove icode entry. Old resids are still valid.
#                for a in self[chainid][resid][icode].itervalues():
#                    a['icode'] = ' '
#
#                # recreate struct dictionary.
#                self.create_struct()

    def remove_altloc(self):
        del_list = []
        for atm in self.atoms:
            if not atm['altloc'] in [' ', 'A']:
                del_list.append(atm)
        if del_list:
            self.del_atom(del_list)

    def rename(self, old_resname, new_resname):
        for seg in self.struct.iter_segments():
            for res in seg.iter_residues():
                if res.resname == old_resname:
                    for atm in res.iter_atoms():
                        atm['resname'] = new_resname
        self.create_struct()
        return 1

    def get_size(self):
        min = list(self.atoms[0]['coord'])
        max = list(self.atoms[0]['coord'])
        for atom in self.atoms:
            coord = atom['coord']
            for dim in range(3):
                if coord[dim] < min[dim]:
                    min[dim] = coord[dim]
                elif coord[dim] > max[dim]:
                    max[dim] = coord[dim]
        return (min, max)


    # def sasa(self):
    #     import sys
    #     sys.path.append('/user/tmeyer/workspace/script/protein_toolbox/asa')
    #
    #     import molecule
    #     import asa
    #
    #     mol = molecule.Molecule(self)
    #     atoms = mol.atoms()
    #     molecule.add_radii(atoms)
    #
    #     # n_sphere = 960
    #     n_sphere = 200
    #     asas = asa.calculate_asa(atoms, 1.4, n_sphere)
    #     print "%.1f angstrom squared." % sum(asas)
    #
    #     for asa, atom in zip(asas, atoms):
    #         atom.bfactor = asa
    #     mol.write_pdb(args[1])

    def sasa(self, estimate=True, tmp_folder=None, sasa_bin=None, no_water=True):
        if sasa_bin is None:
            sasa_bin = '/scratch/scratch/tmeyer/kbplus2/sasa_static.bin'
        # if tmp_folder is None:
        #     tmp_folder = tempfile.mkdtemp() + '/'
            # tmp_folder = '/tmp/'

        if estimate:
            default_3 = {'CAL' : 1.3}
            default_2 = {'CL' : 1.9}
            default_1 = {'C' : 2.0,\
	                     'O' : 1.7,\
                         'N' : 1.8,\
	                     'H' : 1.2,\
	                     'S' : 1.2,\
	                     'Z' : 1.1,\
	                     'P' : 2.15\
                        }
            for atm in self.atoms:
                if atm['vdw'] is None:
                    name = atm['name']
                    if name[0:3] in default_3.keys():
                        atm['vdw'] = default_3[name[0:3]]
                    elif name[0:2] in default_2.keys():
                        atm['vdw'] = default_2[name[0:2]]
                    elif name[0] in default_1.keys():
                        atm['vdw'] = default_1[name[0]]
                    else:
                        print("No default radius defined for %s." % name)
                atm['charge'] = 0.0
            self.par_read = 1

        if tmp_folder is None:
            f, filename = tempfile.mkstemp()
            # File handle is not needed
        else:
            filename = tmp_folder + 'sasa_no_water.pqr'
        self.write_pqr(filename, no_water=no_water)

        process = subprocess.Popen("tcsh",\
                shell=True,\
                stdin=subprocess.PIPE,\
                stdout=subprocess.PIPE,\
                stderr=subprocess.PIPE,\
                )
        commands = []
        commands.append("%s %s %f\n" % (sasa_bin, filename, 1.4))
        commands.append("exit\n")

        for command in commands:
            process.stdin.write(command)

        output = []
        while True:
            next_line = process.stdout.readline()

            # if shell terminates
            if not next_line:
                process = None
                break

            output.append(next_line)

        # for line in process.stderr:
        #     print line

        # Atom 5 = 8.95333 A^2
        atm_sasa = {}
        reg = re.compile(r'Atom (\d+) = ([.\d]+) A\^2')
        for line in output:
            reg_m = reg.match(line)
            if reg_m is not None:
                index = int(reg_m.groups()[0])
                sasa  = float(reg_m.groups()[1])
                atm_sasa[index] = sasa

        c = 0
        for atm in self.atoms:
            if atm['resname'] not in ['HOH', 'TIP3']:
                atm['sasa'] = atm_sasa[c]
                c += 1

        self.create_struct()

        if tmp_folder is None:
            os.remove(filename)

    def create_interaction_graph(self, filename, inet, change_occupancy=False):
        if inet.residue_names is None:
            str = "inet.residue_names must be set to create interaction graph!"
            raise AssertionError(str)

        # Locate all residues in this structure.
        residues = []
        residues_center = []
        for residue_name in inet.residue_names:
            (resname, resid, segname) = re.split(r'[-_]' ,residue_name)
            resid = int(resid)
            residue = self.struct[segname][resid]
            residues.append(residue)
            # Find geometric center.
            center = np.zeros(3)
            nr_of_atoms = 0
            for atom in residue.iter_atoms():
                center += atom['coord']
                nr_of_atoms += 1
            center = np.dot(center, 1.0/nr_of_atoms)
            # print center
            residues_center.append(center)
        
        template_atom = Atom()        
        template_atom["name"] = 'O'
        template_atom["type"] = 'O'
        template_atom["resname"] = 'INT'
        template_atom["segname"] = 'A'
        template_atom["chain"] = 'A'
        template_atom["hetatm"] = True

        if change_occupancy:
            for atom in self.struct.iter_atoms():
                atom["occupancy"] = 0.0
        self.create_struct()

        ipdb = Simple_struct_parser()
        c = 0
        new_resid = 1
        for cluster in inet.last_clusters:
            for i, res1 in enumerate(cluster):
                residue1 = residues[res1]
                center1 = residues_center[res1]


                # The next step is only required, if there are interaction partners.
                if inet.last_interactions.has_key(res1):
                    if change_occupancy:
                        for atom in residue1.iter_atoms():
                            atom['occupancy'] = new_resid / 100.0
                        new_resid += 1
                else:
                    continue

                for res2 in cluster[i+1:]:
                    if not inet.last_interactions[res1].has_key(res2):
                        continue

                    residue2 = residues[res2]
                    center2 = residues_center[res2]
                    
                    new_atom1 = Atom(template_atom)
                    new_atom1["resid"] = c
                    new_atom1["index"] = c
                    c1 = c
                    c += 1
                    new_atom1["coord"] = center1
                    # Store the strongest interaction in units of 100 kJ/mol
                    new_atom1["occupancy"] = inet.last_interactions[res1][res2] / 100.0
                    ipdb.new_atom(new_atom1)

                    new_atom2 = Atom(template_atom)
                    new_atom2["resid"] = c
                    new_atom2["index"] = c
                    c2 = c
                    c += 1
                    new_atom2["coord"] = center2
                    # Store the strongest interaction in units of 100 kJ/mol
                    new_atom2["occupancy"] = inet.last_interactions[res1][res2] / 100.0
                    ipdb.new_atom(new_atom2)

                    ipdb.connections.append((c1, c2))
                    ipdb.connections.append((c2, c1))
        ipdb.crd_read = True
        ipdb.create_struct()
        ipdb.write_pdb(filename)

class MDAnalysis_ssp(object):
    def __init__(self, ssp, tmp_folder="/tmp/mdanalysis/"):
        import MDAnalysis
        import os

        if not os.path.exists(tmp_folder):
            os.mkdir(tmp_folder)
        if tmp_folder[-1] != -1:
            tmp_folder += '/'

        if type(ssp) != str:
            self.ssp = ssp

            tmp_filename = tmp_folder + 'tmp_structure_file.pdb'
            ssp.write_pdb(tmp_filename)
        else:
            tmp_filename = ssp

        self.protein = MDAnalysis.Universe(tmp_filename, format='PDB')

    def res_type(self, residue):
        st = residue.find('-')
        resname = residue[0:st]

        acid = ['ASP', 'GLU', 'EPP', 'DPP', 'CTER', 'TYR' ]
        base = ['ARG', 'LYS', 'HIS', 'HSP', 'HSD', 'HSE', 'TYR']

        if resname in acid:
            return 'acid'
        elif resname in base:
            return 'base'
        else:
            assert("Unknown reidue name: %s" % resname)

    def get_residues_in_sb(self, residue_list=[], kb_resnames=False, tyr_acid = False, tyr_base = False, cutoff = 4):
        """
        Returns a dictionary that contains information about residues in salt bridges, residues not in salt bridges and salt bridge pairs

        Lists have the following structure:
        residues_in_salt_bridges = ['LYS-36_A', 'ASP-45_A', ...]
        residues_not_in_salt_bridges = ['GLU-80_A',...]
        residues_in_salt_bridges = [ [(<Residue_LYS, 36>, 'LYS-36_A'), (<Residue_ASP, 45>,'ASP-45_A')], [...], ...]
        where <Residue_RESNAME, RESID> is a residue object of MDAnalysis.Universe class

        Special:
        Tyrosine residue can be considered as and acid (deprotonated) or a "base" when it forms a hydrogen bond
        @tyr_acid: True/False
        @tyr_base:  True/False


        # @return: list, list, list
        """

        #todo: Possibly replace MDanalysis with CHARMM or redefine search

        sel_str_arg = "(resname ARG and (name NH1 or name NH2))"
        sel_str_lys = "(resname LYS and name NZ)"
        sel_str_his = "((resname HIS or resname HSE or resname HSP or resname HSD) and (name ND1 or name NE2))"
        sel_str_base = sel_str_arg + " or " + sel_str_lys + " or " + sel_str_his
        sel_str_asp = "((resname ASP or resname DPP) and (name OD1 or name OD2))"
        sel_str_glu = "((resname GLU or resname EPP) and (name OE1 or name OE2))"
        sel_str_acid =  sel_str_asp + " or " + sel_str_glu

        if tyr_acid:
            sel_str_tyr = "(resname TYR and name OH)"
            sel_str_acid +=  " or " + sel_str_tyr

        if tyr_base:
            sel_str_tyr = "(resname TYR and name OH)"
            sel_str_base +=  " or " + sel_str_tyr

        chainid_used = False
        for residue in residue_list:
            st = residue.find('-')
            en = residue.find('_')
            resname = residue[0:st]
            resid = int(residue[st + 1: en])
            chain_seg = residue[en+1:]

            chain_seg_ws = chain_seg
            if len(chain_seg) == 1:
                chain_seg_ws += '*'
                chainid_used = True

            if resname == 'CTE':
                sel_str = "((resid %i and segid  %s) and (name OT1 or name OT2))" % (resid, chain_seg_ws)
                sel_str_acid += " or " + sel_str
            if resname == 'NTE':
                sel_str = "((resid %i and segid %s) and name N )" % (resid, chain_seg_ws)
                sel_str_base += " or " + sel_str

        sel_acids = self.protein.select_atoms("(%s) and around %i (%s)" % (sel_str_acid,  cutoff, sel_str_base))
        sel_base = self.protein.select_atoms("(%s) and around %i (%s)" % (sel_str_base, cutoff, sel_str_acid))

        # for res in sel_base.residues:
        #     print res
        # for res in sel_acids.residues:
        #     print res

        residues_in_salt_bridges = []
        for res in (sel_acids + sel_base).residues:
            if kb_resnames:
                resname = res.name
                resname = resname.replace('ASP', 'DPP')
                resname = resname.replace('GLU', 'EPP')
                resname = resname.replace('HIS', 'HSP')
                resname = resname.replace('HSE', 'HSP')
                resname = resname.replace('HSD', 'HSP')
            else:
                resname = res.name
            if chainid_used:
                segname = res.segments[0].name[0]
            else:
                segname = res.segments[0].name

            residue = "%s-%i_%s" % (resname, res.resnum, segname)

            if len(residue_list) != 0:
                if residue in residue_list:
                    residues_in_salt_bridges.append( (res, residue) )
                    # residues_in_salt_bridges.append( residue )
            else:
                residues_in_salt_bridges.append( (res, residue) )
                # residues_in_salt_bridges.append( residue )


        residues_not_in_salt_bridges = list(residue_list)
        for res, residue in residues_in_salt_bridges:
            if residue in residues_not_in_salt_bridges:
                residues_not_in_salt_bridges.remove(residue)


        # Group residues by salt bridges.
        salt_bridges = []

        acid_base_sel = "(%s or %s)" % (sel_str_base, sel_str_acid)
        for (res, residue) in residues_in_salt_bridges:
            # if not salt_bridges:
            #     salt_bridges.append( [(res, residue)] )
            #     print(res, residue)
            # else:
            sel_str = "(segid %s and resid %i and %s)" % (res.segments[0].name, res.resnum, acid_base_sel)

            match_found = -1
            to_merge = []
            for i, sb in enumerate(salt_bridges):
                for (sb_res, sb_residue) in sb:
                    sel_str_comp = "(segid %s and resid %i and %s)" % (sb_res.segments[0].name, sb_res.resnum, acid_base_sel)
                    atms = self.protein.select_atoms("%s and around %i %s" % (sel_str, cutoff, sel_str_comp))
                    if len(atms.atoms) > 0 and (self.res_type(residue) != self.res_type(sb_residue)):
                        if match_found == -1:
                            sb.append( (res, residue) )
                            match_found = i
                        else:
                            if match_found != i:
                                # This is the second alt bridge with a hit match. The two salt bridges are merged.
                                to_merge_entry = (min(match_found, i), max(match_found, i))
                                if not to_merge_entry in to_merge:
                                    to_merge.append(to_merge_entry)

            if match_found == -1:
                salt_bridges.append( [(res, residue)] )

            to_delete = []
            for (i, j) in to_merge:
                for entry in salt_bridges[j]:
                    salt_bridges[i].append( tuple(entry) )
                to_delete.append(j)

            to_delete.sort(reverse=True)
            for i in to_delete:
                salt_bridges.pop(i)


        # Convert residues_in_salt_bridges into list of strings with residue descriptions.
        residues_in_salt_bridges = [x[1] for x in residues_in_salt_bridges]

        return (residues_in_salt_bridges, residues_not_in_salt_bridges, salt_bridges)


# 1 -  6        Record name   "ATOM  "
# 7 - 11        Integer       serial       Atom  serial number.
#13 - 16        Atom          name         Atom name.
#17             Character     altLoc       Alternate location indicator.
#18 - 20        Residue name  resName      Residue name.
#22             Character     chainID      Chain identifier.
#23 - 26        Integer       resSeq       Residue sequence number.
#27             AChar         iCode        Code for insertion of residues.
#31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
#39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
#47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
#55 - 60        Real(6.2)     occupancy    Occupancy.
#61 - 66        Real(6.2)     tempFactor   Temperature  factor.
#77 - 78        LString(2)    element      Element symbol, right-justified.
#79 - 80        LString(2)    charge       Charge  on the atom.





if __name__ == '__main__':

    # s = Simple_struct_parser()
    # s.read_pdb('/scratch/scratch/tmeyer/md_pka/md_pka_manager/sb_runs_2/2lzt/kbp/done/frame1/c_pH7_frame1.pqr', is_pqr=True)
    # s.read_pdb('/scratch/scratch/tmeyer/md_pka/md_pka_manager/sb_runs_2/2lzt/kbp/done/frame1/c_pH7_frame1.pqr', is_pqr=False)
    # s.read_par('/scratch/scratch/tmeyer/CHARMM_NAMD/toppar_36/par_all22_prot.inp')
    # s.read_top('/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/top.inp')
    # s.read_top('/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/patches.rtf')

    # top = []
    # top.append("/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/top.inp")
    # top.append("/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/patches.rtf")
    # par = []
    # par.append("/scratch/scratch/tmeyer/CHARMM_NAMD/toppar_36/gbsw/par_all22_prot_gbsw.inp")
    # par.append("/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/patches.prm")


    # s.read_pdb("/scratch/scratch/pdb/pdb_bio_merged/v5/1v54.pdb1")
    # for t in top:
    #     s.read_top(t)
    # for p in par:
    #     print p
    #     s.read_par(p)



    # import charmm
    #
    # filename = "/scratch/scratch/pdb/pdb_bio_merged/v5/1v54.pdb1"
    #
    # workdir = "/tmp/charmm_tmp/"
    # c = charmm.Charmm_manager(workdir=workdir, top=top, par=par)
    # c.add_structure(filename)
    #
    # c.add_decision('rename__HOH_TIP3', 'keep')
    # c.add_decision('gap__model_gaps', 'keep')
    # c.check_structures(quiet=False)


    # s = c.structure
    # for t in top:
    #     s.read_top(t)
    # for p in par:
    #     s.read_par(p)
    # s.write_pqr('/user/tmeyer/temp/structure1.pdb')

    # mda = MDAnalysis_ssp(s)






    # asas = calculate_asa(atoms, 1.4, n_sphere)
    # print "%.1f angstrom squared." % sum(asas)

    # l = generate_sphere_points(960)
    #
    #
    # protein = Simple_struct_parser()
    # atm = Atom()
    # atm.data["hetatm"] = False
    # atm.data["name"] = 'C'
    # atm.data["type"] = 'C'
    # atm.data["resname"] = 'ALA'
    # atm.data["segname"] = 'A'
    # atm.data["resid"] = 1
    # atm.data["icode"] = ' '
    # atm.data["coord"] = None
    # atm.data["index"] = None
    # # atm.data["mass"] = None
    # # atm.data["vdw"] = None
    # # atm.data["charge"] = None
    # atm.data["chain"] = 'A'
    # atm.data["occupancy"] = 1.0
    # atm.data["tempFactor"] = 0.0
    # atm.data["element"] = 'C'
    # atm.data["altloc"] = ' '
    # for a_coor in l:
    #     # print a_coor
    #     atm_new = Atom(atm)
    #     atm_new.data["coord"] = np.array(a_coor)
    #     protein.new_atom(atm_new)
    # protein.crd_read = True
    # protein.write_pdb("/user/tmeyer/temp/structure.pdb")

        




    # s.write_apbs_par('/scratch/scratch/tmeyer/md_pka/md_pka_manager/sb_runs_2/2lzt/kbp/confE/run_apolar/apbs.par')
    # s.write_pqr('/scratch/scratch/tmeyer/md_pka/md_pka_manager/sb_runs_2/2lzt/kbp/confE/run_apolar/apbs.pqr')

    # a = [1,2,2,4,5,6,7,8,9]
    # del_list = []
    # for i,x in enumerate(a):
    #     if x in [2,6]:
    #         del_list.append(i)
    #
    # print a
    # for i in sorted(del_list, reverse=True):
    #     a.pop(i)
    #
    # print a

    s = Simple_struct_parser()
    

    # s.read_pdb('/scratch/scratch/tmeyer/projects/hemaglutinin/H5/md/lig/c_pH7_frame105.reference.pqr', is_pqr=True)
#    s.write_pqr('/scratch/scratch/tmeyer/projects/hemaglutinin/H5/md/lig/input.pqr')
    #s.create_struct()
    #s.resolve_icodes()
    #s.write_pdb('/scratch/scratch/tmeyer/projects/hemaglutinin/H5/kb_runs/3S11/3s11_noi.pdb')
    
    # s.read_par(f_par)
    # s.read_crd(f_crd)
    # s.read_xplor_psf(f_psf_xplor)

    s.read_pdb('/scratch/scratch/tmeyer/md_benchmark/cco_2gsm_in_water.pdb')
    s.write_crd('/scratch/scratch/tmeyer/md_benchmark/cco_2gsm_in_water.crd')


    ## filename = '/scratch/scratch/pdb/pdb_bio_merged/lz/2lzt.pdb1'
    #filename = '/scratch/scratch/tmeyer/md_pka/md_pka_manager/sb_runs_2/1a2p_c_n_acids_c_his/md/protein_in_water.pdb'
    #
    #top = []
    #top.append("/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/top.inp")
    #top.append("/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/patches.rtf")
    #par = []
    #par.append("/scratch/scratch/tmeyer/CHARMM_NAMD/toppar_36/gbsw/par_all22_prot_gbsw.inp")
    # par.append("/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/patches.prm")
    #
    # import charmm
    #
    # # filename = "/scratch/scratch/pdb/pdb_bio_merged/v5/1v54.pdb1"
    #
    # workdir = "/tmp/charmm_tmp/"
    # c = charmm.Charmm_manager(workdir=workdir, top=top, par=par)
    # c.add_structure(filename)
    #
    # c.add_decision('rename__HOH_TIP3', 'keep')
    # c.add_decision('gap__model_gaps', 'keep')
    # c.check_structures(quiet=False)


    # s = c.structure
    # for t in top:
    #     s.read_top(t)
    # for p in par:
    #     s.read_par(p)
    # s.write_apbs_par('/user/tmeyer/temp/out.pqr')

    # s.write_apbs_par('/user/tmeyer/temp/out.pqr')


    #s.read_pdb(filename)
    #
    #s.sasa(estimate=True)
    #
    #print s.struct['A'][86].sasa
    #print s.struct['A'][75].sasa
    #print s.struct['A'][59].sasa

# test for simple_struct_parser:
#f_crd = '/user/tmeyer/workspace/projects/hämaglutinin/H3_1HGD_with_sb/1hgd-11-1/new_1hgd-11-1-final.crd'
#f_psf_xplor = '/user/tmeyer/workspace/projects/hämaglutinin/H3_1HGD_with_sb/1hgd-11-1/new_1hgd-11-1-final.xplor.psf'
#f_par = '/user/tmeyer/workspace/projects/hämaglutinin/H3_1HGD_with_sb/1hgd-11-1/par_all27_prot_na.prm'
#struct = simple_struct_parser()
#
#struct.read_par(f_par)
#struct.read_crd(f_crd)
#struct.read_xplor_psf(f_psf_xplor)
#
#struct.write_pqr('/user/tmeyer/out.pqr')


#print struct.atoms[5037]['resid']
#print struct.atoms[5037]['name']
#print struct.atoms[5037]['coord']
#print struct.atoms[5037]['mass']
        
