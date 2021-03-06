#!/usr/bin/env python3
# In this module there are stored the customized classes for the Model and Chain objects that are included in the Bio.PDB module from Biopython.
#Bio.PDB is a Biopython module that focuses on working with crystal structures of biological macromolecules.

from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio import pairwise2
from comparch.utilities import protein_dict, dna_dict, rna_dict

class UpdModel(Model):
    """Custom  model class that inherits attributes and methods from the biopython model class.
    Considerations: it allows to have chains with repeated IDs."""
    def add(self, entity):
        """Adds a child to the Entity."""
        entity.set_parent(self)
        self.child_list.append(entity)


class UpdChain(Chain):
    """Custom chain class that inherits attributes and methods from the biopython chain class."""
    protein_dict = protein_dict
    dna_dict = dna_dict
    rna_dict = rna_dict

    def __init__(self, chainObject):
        """Method to initialize the object after its creation. The attributes are needed to correctly set the relationships between parent and child classes.
        Considerations: in order to set CustomChain as child of CustomModel, CustomChain will have to remove
        its original biopython chain model as parent."""
        self.child_list = chainObject.child_list
        self._id = chainObject._id
        self.parent = chainObject.parent
        self.xtra = chainObject.xtra
        self.level = chainObject.level


    def __hash__(self):
        """Returns the hash value of the object. Used to avoid the existence of more than one instance of the object with the same ID in dictionaries or sets."""
        return self._id.__hash__()


    def get_sequence(self):
        """Checks what type of sequence has the chain according to the nature of its first residue.
        Reads the sequence in the format present in the object and returns it in a more conventional format."""
        sequence = ""
        first_residue = self.child_list[0].resname.strip()

        if first_residue not in self.protein_dict:
            if "D" in first_residue:
                for res in self:
                    if res.id[0] == " ":
                        sequence += self.dna_dict[res.resname.strip()]
            else:
                for res in self:
                    if res.id[0] == " ":
                        sequence += self.rna_dict[res.resname.strip()]
        else:
            for res in self:
                if res.id[0] == " ":
                    sequence += self.protein_dict[res.resname]
        return sequence


    def has_equivalent(self, other_sequence):
        """Checks if the sequence of the current chain object is equivalent (very similar or identical) to another sequence.
        Considerations: pairwise alignment is performed and the cutoff is set to 0.95. If the alignment score is above the cutoff,
        the function returns True. Otherwise it returns False."""
        cutoff = 0.95
        sequence1 = self.get_sequence()
        sequence2 = other_sequence

        alignment = pairwise2.align.globalxx(sequence1, sequence2)[0]
        align_score = alignment[2]
        align_length = (len(alignment[0]))

        cutoff_chain = align_score / align_length
        if cutoff_chain > cutoff:
            return True
        else:
            return False

    def get_common_atoms(self, other):
        """Compares the list of atoms of two chains and returns an even tuple of atoms."""
        self_atoms = sorted(self.get_atoms())
        other_atoms = sorted(other.get_atoms())
        len_self = len(self_atoms)
        len_other = len(other_atoms)

        if len_self > len_other:
            return self_atoms[:len_other], other_atoms
        elif len_other > len_self:
            return self_atoms, other_atoms[:len_self]
        else:
            return self_atoms, other_atoms
