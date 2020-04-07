#!/usr/bin/env python3
# In this module there are stored the custom exception classes that are used in the application.


class Directory_Not_Found(OSError):
    """Custom OSError exception that is raised when the input directory with PDB files does not exist."""
    def __init__(self, directory):
        """Method to initialize the object after its creation.

        Keyword arguments:
        directory -- name of the directory that does not exist."""
        self.directory = directory

    def __str__(self):
        """Returns the string representation of the object with the error message."""
        return "Directory %s does not exist. Please select a valid directory." %(self.directory)


class FASTA_Not_Found(OSError):
    """Custom OSError exception that is raised when the input FASTA file, if provided, does not exist."""
    def __init__(self, fasta_file):
        """Method to initialize the object after its creation.

        Keyword arguments:
        fasta_file -- name of the FASTA file that does not exist."""
        self.fasta_file = fasta_file

    def __str__(self):
        """Returns the string representation of the object with the error message."""
        return "FASTA file %s does not exist. Please select a valid FASTA file." %(self.fasta_file)


class No_Input_Directory_Provided(ValueError):
    """Custom ValueError exception that is raised when a directory is not provided as input."""
    pass

class No_PDB_Found(ValueError):
    """Custom ValueError exception that is raised when the input directory does not have PDB files."""
    pass

class Incorrect_Number_Chains(ValueError):
    """Custom ValueError exception that is raised when the PDB model object contains more than 2 chains."""
    pass
