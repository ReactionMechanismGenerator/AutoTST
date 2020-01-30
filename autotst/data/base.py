#!/usr/bin/python
# -*- coding: utf-8 -*-

##########################################################################
#
#   AutoTST - Automated Transition State Theory
#
#   Copyright (c) 2015-2020 Richard H. West (r.west@northeastern.edu)
#   and the AutoTST Team
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
##########################################################################

"""
This module contains functionality for working with transition state family functional
groups, including support for using group additivity to estimate TS geometries.
"""

import os
import logging
import codecs
import numpy as np
import cclib.io
from copy import deepcopy
import rdkit
import rmgpy.data.base
import rmgpy.reaction
import rmgpy.molecule
import rmgpy.data.kinetics.groups 
import rmgpy.species
import rmgpy.reaction 

##########################################################################


class QMData():
    """
    A class that acts as a container for .ts objects
    """

    def __init__(self,
                 ground_state_degeneracy=0,
                 number_of_atoms=0,
                 steric_energy=None,
                 molecular_mass=(0, "amu"),
                 energy=(0, 'eV/molecule'),
                 atomic_numbers=np.array([]),
                 rotational_constants=([], "cm^-1"),
                 atom_coords=([[]], "angstrom"),
                 frequencies=([], "cm^-1"),
                 source=None,
                 method=None,
                 label=""):

        self.ground_state_degeneracy = ground_state_degeneracy
        self.number_of_atoms = number_of_atoms
        self.steric_energy = steric_energy
        self.molecular_mass = molecular_mass
        self.energy = energy
        self.atomic_numbers = atomic_numbers
        self.rotational_constants = rotational_constants
        self.atom_coords = atom_coords
        self.frequencies = frequencies
        self.source = source
        self.method = method
        self.label = label

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        string = 'QMData('
        string += "ground_state_degeneracy={0!r}, ".format(
            self.ground_state_degeneracy)
        string += "number_of_atoms={0!r}, ".format(self.number_of_atoms)
        string += "steric_energy={0!r}, ".format(self.steric_energy)
        string += "molecular_mass={0!r}, ".format(self.molecular_mass)
        string += "energy={0!r}, ".format(self.energy)
        string += "atomic_numbers={0}, ".format(
            "{0!r}".format(self.atomic_numbers).replace(" ", ""))
        string += "rotational_constants={0}, ".format("{0}".format(
            self.rotational_constants).replace("\n", "").replace(" ", ""))
        string += "atom_coords={0}, ".format("{0}".format(
            self.atom_coords).replace("\n", "").replace(" ", ""))
        string += "frequencies={0}, ".format("{0}".format(
            self.frequencies).replace("\n", "").replace(" ", ""))
        string += "source={0!r}, ".format(self.source)
        string = string[:-2] + ')'
        return string

    def get_qmdata(self, file_path=None):
        "A helper function to fill in the qmdata using CCLib"

        parser = cclib.io.ccread(file_path, loglevel=logging.ERROR)

        self.ground_state_degeneracy = parser.mult
        self.atom_numbers = parser.atomnos
        self.number_of_atoms = len(parser.atomnos)
        self.atom_coords = (parser.atomcoords[-1], "angstrom")
        self.steric_energy = None  # Need to fix this
        self.molecular_mass = (parser.atommasses.sum(), "amu")
        self.energy = (parser.scfenergies[-1], "eV/molecule")
        self.atomic_numbers = parser.atomnos
        self.rotational_constants = ([], "cm^-1")  # Need to fix this
        self.frequencies = (parser.vibfreqs, "cm^-1")
        self.source = "AutoTST"
        self.method = parser.metadata["functional"]


class DistanceData():
    """
    A class for storing distance matrix data, for geometry estimation
    """

    def __init__(self, distances={}, uncertainties=None, method=None):
        self.distances = distances
        self.uncertainties = uncertainties
        self.method = method
        self.comment = ''
        assert isinstance(distances, dict), "distances should be a dict"
        if method:
            assert isinstance(method, str), "method should be a string"

    def __repr__(self):
        strings = ["DistanceData("]

        strings.append("distances={")
        for key in sorted(self.distances.keys()):
            strings.append("{0!r}: {1:.6f},".format(key, self.distances[key]))
        strings.append("}")

        if self.uncertainties is not None:
            strings.append(", uncertainties={")
            for key in sorted(self.uncertainties.keys()):
                strings.append("{0!r}: {1:.6f},".format(
                    key, self.uncertainties[key]))
            strings.append("}")

        if self.method:
            strings.append(", method={0!r}".format(self.method))
        if self.comment:
            strings.append(", comment={0!r}".format(self.comment))
        strings.append(")")
        return ''.join(strings)

    def add(self, other):
        """Adds the `other` distances to these."""
        assert len(
            self.distances) == len(
            other.distances), "self and other must have the same size dictionary of distances, but self={0!r} and other={1!r}".format(
            self, other)
        for key, value in other.distances.items():
            self.distances[key] += value
        if self.uncertainties and other.uncertainties:
            for key, value in other.uncertainties.items():
                self.uncertainties[key] += value
        else:
            self.uncertainties = None

    def __copy__(self):
        return DistanceData(
            distances=self.distances.copy(),
            uncertainties=self.uncertainties.copy(),
            method=self.method)


class TransitionStates(rmgpy.data.base.Database):
    """
    loads, and contains, both Depository and Groups
    """

    def __init__(self):
        self.groups = None
        self.depository = None
        self.family = None

    def load(self, path, local_context, global_context):
        """
        Load the TS database
        """
        if local_context is None:
            local_context = {}
            local_context['DistanceData'] = DistanceData

        fpath = os.path.join(path, 'TS_training', 'reactions.py')

        logging.debug(
            "Loading transitions state family training set from {0}".format(fpath))

        depository = TransitionStateDepository(
            label='{0}/TS_training'.format(path.split('/')[-1]))  # 'intra_H_migration/TS_training')
        depository.load(fpath, local_context, global_context)
        self.depository = depository

        fpath = os.path.join(path, 'TS_groups.py')
        logging.info(
            "Loading transitions state family groups from {0}".format(fpath))
        # 'intra_H_migration/TS_groups')
        groups = TSGroups(label='{0}/TS_groups'.format(path.split('/')[-1]))
        groups.load(fpath, local_context, global_context)

        self.family.forward_template.reactants = [
            groups.entries[entry.label] for entry in self.family.forward_template.reactants]
        # self.family.forward_template.products = [groups.entries[entry.label] for entry in self.family.forward_template.products]
        self.family.entries = groups.entries
        self.family.groups = groups
        groups.num_reactants = len(self.family.forward_template.reactants)
        self.groups = groups

    def estimate_distances(self, reaction):
        """
        Return estimated DistanceData for the given reaction
        """
        # Should check depository first, but for now just go straight to group
        # additive estimate:
        return self.groups.estimate_distances_using_group_additivity(reaction)

    def save_transition_state_groups(self, path, entry_name='entry'):
        """
        Save the current database to the file at location `path` on disk. The
        optional `entryName` parameter specifies the identifier used for each
        data entry.
        """
        entries = self.groups.get_entries_to_save()

        # Write the header
        f = codecs.open(path, 'w', 'utf-8')
        f.write('#!/usr/bin/env python\n')
        f.write('# encoding: utf-8\n\n')
        f.write('name = "{0}"\n'.format(self.groups.name))
        f.write('short_desc = u"{0}"\n'.format(self.groups.short_desc))
        f.write('long_desc = u"""\n')
        f.write(self.groups.long_desc)
        f.write('\n"""\n\n')

        # Save the entries
        for entry in entries:
            self.save_entry(f, entry)

        # Write the tree
        if len(self.groups.top) > 0:
            f.write('tree(\n')
            f.write('"""\n')
            f.write(self.generate_old_tree(self.groups.top, 1))
            f.write('"""\n')
            f.write(')\n\n')

        f.close()
        return

    def save_entry(self, f, entry):
        """
        Save an `entry` in the kinetics database by writing a string to
        the given file object `f`.
        """
        import arkane.output

        def sort_efficiencies(efficiencies0):
            efficiencies = {}
            for mol, eff in efficiencies0.items():
                if isinstance(mol, str):
                    # already in SMILES string format
                    smiles = mol
                else:
                    smiles = mol.to_smiles()

                efficiencies[smiles] = eff
            keys = list(efficiencies.keys())
            keys.sort()
            return [(key, efficiencies[key]) for key in keys]

        f.write('entry(\n')
        f.write('    index = {0:d},\n'.format(entry.index))
        if entry.label != '':
            f.write('    label = "{0}",\n'.format(entry.label))

        # Entries for kinetic rules, libraries, training reactions
        # and depositories will have a rmgpy.reaction.Reaction object for its item
        if isinstance(entry.item, rmgpy.reaction.Reaction):
            # Write out additional data if depository or library
            # kinetic rules would have a Group object for its reactants instead of Species
            if isinstance(entry.item.reactants[0], rmgpy.species.Species):
                # Add degeneracy if the reaction is coming from a depository or kinetics library
                f.write('    degeneracy = {0:.1f},\n'.format(
                    entry.item.degeneracy))
                if entry.item.duplicate:
                    f.write('    duplicate = {0!r},\n'.format(
                        entry.item.duplicate))
                if not entry.item.reversible:
                    f.write('    reversible = {0!r},\n'.format(
                        entry.item.reversible))
                if entry.item.allow_pdep_route:
                    f.write('    allow_pdep_route = {0!r},\n'.format(
                        entry.item.allow_pdep_route))
                if entry.item.elementary_high_p:
                    f.write('    elementary_high_p = {0!r},\n'.format(
                        entry.item.elementary_high_p))
                if entry.item.allow_max_rate_violation:
                    f.write('    allow_max_rate_violation = {0!r},\n'.format(
                        entry.item.allow_max_rate_violation))
            # Entries for groups with have a group or logicNode for its item
        elif isinstance(entry.item, rmgpy.molecule.Group):
            f.write('    group = \n')
            f.write('"""\n')
            f.write(entry.item.to_adjacency_list())
            f.write('""",\n')
        elif isinstance(entry.item, rmgpy.data.base.LogicNode):
            f.write('    group = "{0}",\n'.format(entry.item))
        else:
            raise rmgpy.data.base.DatabaseError(
                "Encountered unexpected item of type {0} while saving database.".format(entry.item.__class__))

        # Write distances
        if isinstance(entry.data, DistanceData):
            data_str = arkane.output.prettify(repr(entry.data))
            data_str = data_str.replace('\n', '\n    ')
            f.write('    distances = {},\n'.format(data_str))
        else:
            assert False

        # Write reference
        if entry.reference is not None:
            reference = entry.reference.to_pretty_repr()
            lines = reference.splitlines()
            f.write('    reference = {0}\n'.format(lines[0]))
            for line in lines[1:-1]:
                f.write('    {0}\n'.format(line))
            f.write('    ),\n'.format(lines[0]))

        if entry.reference_type != "":
            f.write('    reference_type = "{0}",\n'.format(entry.reference_type))
        if entry.rank is not None:
            f.write('    rank = {0},\n'.format(entry.rank))

        if entry.short_desc.strip() != '':
            f.write('    short_desc = u"""')
            try:
                f.write(entry.short_desc.encode('utf-8'))
            except:
                f.write(entry.short_desc.strip().encode(
                    'ascii', 'ignore') + "\n")
            f.write('""",\n')

        if entry.long_desc.strip() != '':
            f.write('    long_desc = \n')
            f.write('u"""\n')
            try:
                f.write(entry.long_desc.strip().encode('utf-8') + "\n")
            except:
                f.write(entry.long_desc.strip().encode(
                    'ascii', 'ignore') + "\n")
            f.write('""",\n')

        f.write(')\n\n')

        return


##########################################################################


class TransitionStateDepository(rmgpy.data.base.Database):
    """
    A class for working with an RMG transition state depository. Each depository
    corresponds to a reaction family (a :class:`KineticsFamily` object). Each
    entry in a transition state depository involves a reaction defined either by a
    real reactant and product species.
    """

    def __init__(self, label='', name='', short_desc='', long_desc=''):
        rmgpy.data.base.Database.__init__(self, label=label, name=name,
                          short_desc=short_desc, long_desc=long_desc)

    def __repr__(self):
        return '<TransitionStateDepository "{0}">'.format(self.label)

    def load(self, path, local_context=None, global_context=None):

        rmgpy.data.base.Database.load(self, path, local_context, global_context)

        # Load the species in the kinetics library
        species_dict = self.get_species(os.path.join(
            os.path.dirname(path), 'dictionary.txt'))
        # Make sure all of the reactions draw from only this set
        entries = list(self.entries.values())
        for entry in entries:
            # Create a new reaction per entry
            rxn = entry.item
            rxn_string = entry.label
            # Convert the reactants and products to Species objects using the
            # speciesDict
            reactants, products = rxn_string.split('=')
            reversible = True
            if '<=>' in rxn_string:
                reactants = reactants[:-1]
                products = products[1:]
            elif '=>' in rxn_string:
                products = products[1:]
                reversible = False
            assert reversible == rxn.reversible
            for reactant in reactants.split('+'):
                reactant = reactant.strip()
                if reactant not in species_dict:
                    raise rmgpy.data.base.DatabaseError(
                        'RMGSpecies {0} in kinetics depository {1} is missing from its dictionary.'.format(
                            reactant, self.label))
                # For some reason we need molecule objects in the depository
                # rather than species objects
                rxn.reactants.append(species_dict[reactant])
            for product in products.split('+'):
                product = product.strip()
                if product not in species_dict:
                    raise rmgpy.data.base.DatabaseError(
                        'RMGSpecies {0} in kinetics depository {1} is missing from its dictionary.'.format(
                            product, self.label))
                # For some reason we need molecule objects in the depository
                # rather than species objects
                rxn.products.append(species_dict[product])

            if not rxn.is_balanced():
                raise rmgpy.data.base.DatabaseError(
                    'Reaction {0} in kinetics depository {1} was not balanced! Please reformulate.'.format(
                        rxn, self.label))

    def load_entry(self,
                  index,
                  reactant1=None,
                  reactant2=None,
                  reactant3=None,
                  product1=None,
                  product2=None,
                  product3=None,
                  distances=None,
                  degeneracy=1,
                  label='',
                  duplicate=False,
                  reversible=True,
                  reference=None,
                  reference_type='',
                  short_desc='',
                  long_desc='',
                  rank=None,
                  ):

        reaction = rmgpy.reaction.Reaction(reactants=[], products=[
        ], degeneracy=degeneracy, duplicate=duplicate, reversible=reversible)

        entry = rmgpy.data.base.Entry(
            index=index,
            label=label,
            item=reaction,
            data=distances,
            reference=reference,
            reference_type=reference_type,
            short_desc=short_desc,
            long_desc=long_desc.strip(),
            rank=rank,
        )
        self.entries['{0:d}:{1}'.format(index, label)] = entry
        return entry

    def save_entry(self, f, entry):
        """
        Save an `entry` in the kinetics database by writing a string to
        the given file object `f`.
        """
        import arkane.output

        def sort_efficiencies(efficiencies0):
            efficiencies = {}
            for mol, eff in efficiencies0.items():
                if isinstance(mol, str):
                    # already in SMILES string format
                    smiles = mol
                else:
                    smiles = mol.to_smiles()

                efficiencies[smiles] = eff
            keys = list(efficiencies.keys())
            keys.sort()
            return [(key, efficiencies[key]) for key in keys]

        f.write('entry(\n')
        f.write('    index = {0:d},\n'.format(entry.index))
        if entry.label != '':
            f.write('    label = "{0}",\n'.format(entry.label))

        # Entries for kinetic rules, libraries, training reactions
        # and depositories will have a Reaction object for its item
        if isinstance(entry.item, rmgpy.reaction.Reaction):
            # Write out additional data if depository or library
            # kinetic rules would have a Group object for its reactants instead of Species
            if isinstance(entry.item.reactants[0], rmgpy.species.Species):
                # Add degeneracy if the reaction is coming from a depository or kinetics library
                f.write('    degeneracy = {0:.1f},\n'.format(
                    entry.item.degeneracy))
                if entry.item.duplicate:
                    f.write('    duplicate = {0!r},\n'.format(
                        entry.item.duplicate))
                if not entry.item.reversible:
                    f.write('    reversible = {0!r},\n'.format(
                        entry.item.reversible))
                if entry.item.allow_pdep_route:
                    f.write('    allow_pdep_route = {0!r},\n'.format(
                        entry.item.allow_pdep_route))
                if entry.item.elementary_high_p:
                    f.write('    elementary_high_p = {0!r},\n'.format(
                        entry.item.elementary_high_p))
                if entry.item.allow_max_rate_violation:
                    f.write('    allow_max_rate_violation = {0!r},\n'.format(
                        entry.item.allow_max_rate_violation))
            # Entries for groups with have a group or logicNode for its item
        elif isinstance(entry.item, rmgpy.molecule.Group):
            f.write('    group = \n')
            f.write('"""\n')
            f.write(entry.item.to_adjacency_list())
            f.write('""",\n')
        elif isinstance(entry.item, rmgpy.data.base.LogicNode):
            f.write('    group = "{0}",\n'.format(entry.item))
        else:
            raise rmgpy.data.base.DatabaseError(
                "Encountered unexpected item of type {0} while saving database.".format(entry.item.__class__))

        # Write distances
        if isinstance(entry.data, DistanceData):
            data_str = arkane.output.prettify(repr(entry.data))

            data_str = data_str.replace('\n', '\n   ')

            f.write('    distances = {},\n'.format(data_str))
        else:
            assert False

        # Write reference
        if entry.reference is not None:
            reference = entry.reference.to_pretty_repr()
            lines = reference.splitlines()
            f.write('    reference = {0}\n'.format(lines[0]))
            for line in lines[1:-1]:
                f.write('    {0}\n'.format(line))
            f.write('    ),\n'.format(lines[0]))

        if entry.reference_type != "":
            f.write('    reference_type = "{0}",\n'.format(entry.reference_type))
        if entry.rank is not None:
            f.write('    rank = {0},\n'.format(entry.rank))

        if entry.short_desc.strip() != '':
            f.write('    short_desc = u"""')
            try:
                f.write(entry.short_desc.encode('utf-8'))
            except:
                f.write(entry.short_desc.strip().encode(
                    'ascii', 'ignore') + "\n")
            f.write('""",\n')

        if entry.long_desc.strip() != '':
            f.write('    long_desc = \n')
            f.write('u"""\n')
            try:
                f.write(entry.long_desc.strip().encode('utf-8') + "\n")
            except:
                f.write(entry.long_desc.strip().encode(
                    'ascii', 'ignore') + "\n")
            f.write('""",\n')

        f.write(')\n\n')

        return


##########################################################################

class TSGroups(rmgpy.data.base.Database):
    """
    A class for working with group additivity values for transition state distances.
    """

    def __init__(self,
                 entries=None,
                 top=None,
                 label='',
                 name='',
                 short_desc='',
                 long_desc='',
                 forward_template=None,
                 forward_recipe=None,
                 reverse_template=None,
                 reverse_recipe=None,
                 forbidden=None
                 ):
        rmgpy.data.base.Database.__init__(self, entries, top, label, name, short_desc, long_desc)
        self.num_reactants = 0

    def __repr__(self):
        return '<TSGroups "{0}">'.format(self.label)

    def load_entry(
            self,
            index,
            label,
            group,
            distances,
            reference=None,
            reference_type='',
            short_desc='',
            long_desc=''):
        if group[0:3].upper() == 'OR{' or group[0:4].upper(
        ) == 'AND{' or group[0:7].upper() == 'NOT OR{' or group[0:8].upper() == 'NOT AND{':
            item = rmgpy.data.base.make_logic_node(group)
        else:
            item = rmgpy.molecule.Group().from_adjacency_list(group)
        self.entries[label] = rmgpy.data.base.Entry(
            index=index,
            label=label,
            item=item,
            data=distances,
            reference=reference,
            reference_type=reference_type,
            short_desc=short_desc,
            long_desc=long_desc.strip(),
        )

    def get_reaction_template(self, reaction):
        """
        For a given `reaction` with properly-labeled :class:`Molecule` objects
        as the reactants, determine the most specific nodes in the tree that
        describe the reaction.
        """
        # from .family import TemplateReaction
        #assert isinstance(reaction, TemplateReaction), "Can only match TemplateReactions"
        # Get forward reaction template and remove any duplicates
        forward_template = self.top[:]
        temporary = []
        symmetric_tree = False
        for entry in forward_template:
            if entry not in temporary:
                temporary.append(entry)
            else:
                # duplicate node found at top of tree
                # eg. R_recombination: ['Y_rad', 'Y_rad']
                assert len(
                    forward_template) == 2, 'Can currently only do symmetric trees with nothing else in them'
                symmetric_tree = True
        forward_template = temporary

        # Descend reactant trees as far as possible
        template = []
        for entry in forward_template:
            # entry is a top-level node that should be matched
            group = entry.item

            # To sort out "union" groups, descend to the first child that's not a logical node
            # ...but this child may not match the structure.
            # eg. an R3 ring node will not match an R4 ring structure.
            # (but at least the first such child will contain fewest labels - we hope)
            if isinstance(entry.item, rmgpy.data.base.LogicNode):
                group = entry.item.get_possible_structures(self.entries)[0]

            # list of atom labels in highest non-union node
            atom_list = group.get_all_labeled_atoms()

            for reactant in reaction.reactants:
                if isinstance(reactant, rmgpy.species.Species):
                    reactant = reactant.molecule[0]
                # Match labeled atoms
                # Check this reactant has each of the atom labels in this group
                if not all([reactant.contains_labeled_atom(label)
                            for label in atom_list]):
                    continue  # don't try to match this structure - the atoms aren't there!
                # Match structures
                atoms = reactant.get_all_labeled_atoms()

                matched_node = self.descend_tree(reactant, atoms, root=entry)
                if matched_node is not None:
                    template.append(matched_node)
                # else:
                #    logging.warning("Couldn't find match for {0} in {1}".format(entry,atomList))
                #    logging.warning(reactant.to_adjacency_list())

        # Get fresh templates (with duplicate nodes back in)
        forward_template = self.top[:]
        if self.label.lower().startswith('r_recombination'):
            forward_template.append(forward_template[0])

        # Check that we were able to match the template.
        # template is a list of the actual matched nodes
        # forwardTemplate is a list of the top level nodes that should be
        # matched
        if len(template) != len(forward_template):
            logging.warning(
                'Unable to find matching template for reaction {0} in reaction family {1}'.format(
                    str(reaction), str(self)))
            logging.warning(" Trying to match " + str(forward_template))
            logging.warning(" Matched " + str(template))
            print(str(self), template, forward_template)
            for n, reactant in enumerate(reaction.reactants):
                print("Reactant", n)
                print(reactant.to_adjacency_list() + '\n')
            for n, product in enumerate(reaction.products):
                print("Product", n)
                print(product.to_adjacency_list() + '\n')
            raise rmgpy.data.kinetics.groups.KineticsError(reaction)

        for reactant in reaction.reactants:
            if isinstance(reactant, rmgpy.species.Species):
                reactant = reactant.molecule[0]
            # reactant.clear_labeled_atoms()

        return template

    def estimate_distances_using_group_additivity(self, reaction):
        """
        Determine the appropriate transition state distances for a reaction
        with the given `template` using group additivity.
        """

        template = self.get_reaction_template(reaction)
        reference_distances = self.top[0].data  # or something like that

        # Start with the generic distances of the top-level nodes
        # Make a copy so we don't modify the original
        ts_distances = deepcopy(reference_distances)

        # Now add in more specific corrections if possible
        for entry in template:
            comment_line = "Matched node "
            while not entry.data.distances and entry not in self.top:
                # Keep climbing tree until you find a (non-top) node with
                # distances.
                comment_line += "{0} >> ".format(entry.label)
                entry = entry.parent
            if entry.data.distances and entry not in self.top:
                ts_distances.add(entry.data)
                comment_line += "{0} ({1})".format(entry.label,
                                                   entry.long_desc.split('\n')[0])
            elif entry in self.top:
                comment_line += "{0} (Top node)".format(entry.label)
            ts_distances.comment += comment_line + '\n'

        return ts_distances

    def generate_group_additivity_values(self, training_set):
        """
        Generate the group additivity values using the given `training_set`,
        a list of 2-tuples of the form ``(template, kinetics)``. You must also
        specify the `kunits` for the family and the `method` to use when
        generating the group values. Returns ``True`` if the group values have
        changed significantly since the last time they were fitted, or ``False``
        otherwise.
        """
        # keep track of previous values so we can detect if they change
        old_entries = dict()
        for label, entry in list(self.entries.items()):
            if entry.data is not None:
                old_entries[label] = entry.data

        # Determine a complete list of the entries in the database, sorted as
        # in the tree
        group_entries = self.top[:]

        for entry in self.top:
            # Entries in the TS_group.py tree
            group_entries.extend(self.descendants(entry))

        # Determine a unique list of the groups we will be able to fit
        # parameters for
        group_list = []
        for template, distances in training_set:
            for group in template:
                if isinstance(group, str):
                    group = self.entries[group]
                if group not in self.top:
                    group_list.append(group)
                    group_list.extend(self.ancestors(group)[:-1])
        group_list = list(set(group_list))
        group_list.sort(key=lambda x: x.index)

        if True:  # should remove this IF block, as we only have one method.
            # Initialize dictionaries of fitted group values and uncertainties
            group_values = {}
            group_uncertainties = {}
            group_counts = {}
            group_comments = {}
            for entry in group_entries:
                group_values[entry] = []
                group_uncertainties[entry] = []
                group_counts[entry] = []
                group_comments[entry] = set()

            # Generate least-squares matrix and vector
            A = []
            b = []

            # ['d12', 'd13', 'd23']
            distance_keys = sorted(training_set[0][1].distances.keys())
            distance_data = []
            for template, distance_data in training_set:
                d = [distance_data.distances[key] for key in distance_keys]
                distance_data.append(d)

                # Create every combination of each group and its ancestors with
                # each other
                combinations = []
                for group in template:
                    groups = [group]
                    # Groups from the group.py tree
                    groups.extend(self.ancestors(group))
                    combinations.append(groups)
                combinations = get_all_combinations(combinations)
                # Add a row to the matrix for each combination
                for groups in combinations:
                    Arow = [1 if group in groups else 0 for group in group_list]
                    Arow.append(1)
                    brow = d
                    A.append(Arow)
                    b.append(brow)

                    for group in groups:
                        if isinstance(group, str):
                            group = self.entries[group]
                        group_comments[group].add("{0!s}".format(template))

            if len(A) == 0:
                logging.warning(
                    'Unable to fit kinetics groups for family "{0}"; no valid data found.'.format(
                        self.label))
                return
            A = np.array(A)
            b = np.array(b)
            distance_data = np.array(distance_data)

            x, residues, rank, s = np.linalg.lstsq(A, b)

            for t, distance_key in enumerate(distance_keys):

                # Determine error in each group
                stdev = np.zeros(len(group_list) + 1, np.float64)
                count = np.zeros(len(group_list) + 1, np.int)

                for index in range(len(training_set)):
                    template, distances = training_set[index]
                    d = np.float64(distance_data[index, t])
                    dm = x[-1, t] + sum([x[group_list.index(group), t]
                                         for group in template if group in group_list])
                    variance = (dm - d)**2
                    for group in template:
                        groups = [group]
                        groups.extend(self.ancestors(group))
                        for g in groups:
                            if g.label not in [top.label for top in self.top]:
                                ind = group_list.index(g)
                                stdev[ind] += variance
                                count[ind] += 1
                    stdev[-1] += variance
                    count[-1] += 1

                import scipy.stats
                ci = np.zeros(len(count))
                for i in range(len(count)):
                    if count[i] > 1:
                        stdev[i] = np.sqrt(stdev[i] / (count[i] - 1))
                        ci[i] = scipy.stats.t.ppf(
                            0.975, count[i] - 1) * stdev[i]
                    else:
                        stdev[i] = None
                        ci[i] = None
                # Update dictionaries of fitted group values and uncertainties
                for entry in group_entries:
                    if entry == self.top[0]:
                        group_values[entry].append(x[-1, t])
                        group_uncertainties[entry].append(ci[-1])
                        group_counts[entry].append(count[-1])
                    elif entry.label in [group.label for group in group_list]:
                        index = [group.label for group in group_list].index(
                            entry.label)
                        group_values[entry].append(x[index, t])
                        group_uncertainties[entry].append(ci[index])
                        group_counts[entry].append(count[index])
                    else:
                        group_values[entry] = None
                        group_uncertainties[entry] = None
                        group_counts[entry] = None

            # Store the fitted group values and uncertainties on the associated
            # entries
            for entry in group_entries:
                if group_values[entry] is not None:
                    if not any(
                        np.isnan(
                            np.array(
                                group_uncertainties[entry]))):
                        # should be entry.data.* (e.g.
                        # entry.data.uncertainties)
                        uncertainties = np.array(group_uncertainties[entry])
                        uncertainty_type = '+|-'
                    else:
                        uncertainties = {}
                    # should be entry.*
                    short_desc = "Fitted to {0} distances.\n".format(
                        group_counts[entry][0])
                    long_desc = "\n".join(group_comments[entry.label])
                    distances_dict = {key: distance for key, distance in zip(
                        distance_keys, group_values[entry])}
                    uncertainties_dict = {
                        key: distance for key, distance in zip(
                            distance_keys, uncertainties)}
                    entry.data = DistanceData(
                        distances=distances_dict,
                        uncertainties=uncertainties_dict)
                    entry.short_desc = short_desc
                    entry.long_desc = long_desc
                else:
                    entry.data = DistanceData()

        changed = False
        for label, entry in list(self.entries.items()):
            if entry.data is not None:
                continue  # because this is broken:
                if label in old_entries:
                    old_entry = old_entries[label][label][0]
                    for key, distance in entry.data.items():
                        diff = 0
                        for k in range(3):
                            diff += abs(distance[0][k] / old_entry[k] - 1)
                        if diff > 0.01:
                            changed = True
                            entry.history.append(event)
            else:
                changed = True
                entry.history.append(event)
        return True  # because the thing above is broken
        return changed
