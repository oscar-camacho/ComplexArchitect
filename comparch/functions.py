#!/usr/bin/env python3
#In this module there are stored the functions required to run the application.

from Bio import SeqIO, PDB, pairwise2
from comparch.pdb_classes import UpdChain, UpdModel
from comparch.exception_classes import *
from comparch.optimize_functions import *
from comparch.utilities import ASCII_characters
import os, sys, re, random, copy


def data_extraction(directory, fasta_file, verbose=False):
    """Takes a directory with PDB files, generates PDB objects of the models for each file and returns a list of the PDB models.
    It also checks if the FASTA file, if provided, exists.

    Keyword arguments:
    directory -- path of the directory where the PDB files are
    fasta_file -- FASTA file with the non-repeated sequences of the chains that conform the complex.
    verbose -- boolean, prints to stderr the progress of the program

    Exceptions raised:
    Directory_Not_Found -- raised if name of the directory does not exist
    No_Input_Directory_Provided -- raised if directory is no provided
    No_PDB_Found -- raised if directory does not contain PDB files
    Incorrect_Number_Chains -- raised if a PDB object contains more than 2 chains
    FASTA_Not_Found -- raised if name of the FASTA file does not exist"""
    if directory is not None:
        if verbose:
            sys.stderr.write("Looking for PDB input files from %s\n" %directory)
        if os.path.isdir(directory):
            path = directory + '/'
            pdb_files = [path + pdb_file for pdb_file in os.listdir(path) if pdb_file[-4:] == '.pdb']
        else:
            raise Directory_Not_Found(directory)
    else:
        raise No_Input_Directory_Provided("""No input directory provided. You should introduce the name of the directory containing the PDB files for each interacting pair of the complex.""")

    num_files = len(pdb_files)   #Printing number of processed filenames
    if verbose:
        sys.stderr.write('%s PDB files found.\n' %(num_files))

    parser = PDB.PDBParser(PERMISSIVE=1, QUIET=True)   #Storing PDB information
    if verbose:
        sys.stderr.write("Storing PDB information...\n")
    if num_files != 0:
        pdb_models = []
        for pdb_file in pdb_files:
            pdb_models.append(parser.get_structure('pdb_model', pdb_file))
            if verbose:
                sys.stderr.write("  %s finished\n" %(pdb_file))
    else:
        raise No_PDB_Found("No PDB files found. Please make sure the given directory contains PDB files.")

    for pdb_model in pdb_models:    #Ensuring all pdb files contain 1 or 2 chains
        if len(pdb_model.child_list) > 2:
            raise Incorrect_Number_Chains("A PDB input file does not contain two chains. Please, all PDB files must only contain two chains.")
    if verbose:
        sys.stderr.write("PDB information stored.\n\n")

    if fasta_file:
        if not os.path.isfile(fasta_file):
            raise FASTA_Not_Found(fasta_file)
    return pdb_models


def new_id(id_list):
    """Checks for characters in an ASCII characters set that are not in a given list of IDs. Returns a new character that will serve as chain ID.

    Keyword arguments:
    id_list -- list of chain IDs that is progressively being apdated"""
    characters = ASCII_characters
    for character in characters:
        if character not in id_list:
            return character


def find_equivalent_chains(pdb_models, fasta_file, verbose=False):
    """Updates original PDB models to customized PDB models and return the list of those. Unifies the ID chains in a way that those chains with sequences
    that are very similar or identical (equivalent chains) receive the same ID.

    Keyword arguments:
    pdb_models -- list of objects of the PDB models with 2 interacting chains, extracted from the PDB files
    verbose -- boolean, prints to stderr the progress of the program

    Considerations: the function works in two different scenarios.
    FASTA file is provided -- pairwise alignments between chain sequences and FASTA sequences are performed. If they are equivalent, the chain receives the ID from the FASTA sequence.
    FASTA file is not provided -- pairwise alignments are performed between the chain sequences. Unique chains are identified and those chains that are equivalent share the same ID."""
    if verbose:
        sys.stderr.write("Identifying equivalent chains and unifying IDs.\n")
        sys.stderr.write("Performing pairwise sequence alignments...\n")
    unique_chains = []
    id_list = []
    for i in range(len(pdb_models)):
        pdb = pdb_models[i]
        pdb_model = UpdModel(str(i))
        for chain in pdb.get_chains():
            chain = UpdChain(chain)
            chain.parent = None
            if fasta_file:
                for record in SeqIO.parse(fasta_file, "fasta"):
                    equivalence = chain.has_equivalent(record.seq)
                    if equivalence:
                        m = re.search("(?<=\:)(.).*?", record.id)
                        id = m.group()
                        chain.id = id
            else:
                if not unique_chains:
                    chain.id = new_id(id_list)
                    id_list.append(chain.id)
                    unique_chains.append(chain)
                else:
                    for unique_chain in unique_chains:
                        equivalence = chain.has_equivalent(unique_chain.get_sequence())
                        if equivalence:
                            chain.id = unique_chain.get_id()
                            break
                        if unique_chain == unique_chains[-1]:
                            chain.id = new_id(id_list)
                            id_list.append(chain.id)
                            unique_chains.append(chain)
            pdb_model.add(chain)
        pdb_models[i] = pdb_model
    if verbose:
        sys.stderr.write("Equivalent chains correctly identified and IDs unified.\n")
    return pdb_models


def get_unique_chains(pdb_models):
    """Collects all the chains with unified IDs from the models and returns a list of unique chains formed only by a signle representant of the chains that shares a specific ID.

    Keyword arguments
    pdb_models -- PDB models containing chains with unified IDs."""
    chain_list = []
    for pdb_model in pdb_models:
        for chain in pdb_model.get_chains():
            chain_list.append(chain)
    unique_chains = set(chain_list)
    sorted_unique_chains = sorted(unique_chains)
    return sorted_unique_chains


def get_stoichiometry(stoich_string):
    """Transforms a string containing the stoichiometry of the complex given by the user into a dictionary where keys are the chain IDs and values
    are the number of times the chain has to be present in the complex. Returns the stoichiometry dictionary.

    Keyword arguments
    stoich_string -- desired stoichiometry for the complex given by the user as a string
    Exceptions:
    ValueError: if the format of the stoichiometry string is incorrect, prints a message and exits the program."""
    stoich_dict = {}
    try:
        stoich_list = stoich_string.split(",")
        for stoich in stoich_list:
            chain, number = stoich.split(":")
            stoich_dict[chain] = int(number)
        return stoich_dict
    except ValueError:
        sys.stderr.write("Stoichiometry format is wrong. Please follow this format: A:1,B:3,C:2, ...\n")
        sys.exit(1)


def get_stoich_without_fasta(pdb_models):
    """Alternative function that asks for the stoichiometry of the complex and return a stoichiometry dictionary.

    Keyword arguments
    pdb_models -- PDB models containing chains with unified IDs according to sequence similarity
    Considerations: If the user does not provide a FASTA file, probably can not identify the working sequences to add a stoichiometry. Then the function shows
    the sequences of the unique chains with its corresponding generated IDs and asks for the stoichiometry taking into account the given IDs."""
    sys.stderr.write("\nATENTION! You did not provide a FASTA file. However, you may want to specify the stoichiometry of the complex.\n")
    sys.stderr.write("In order to correctly identify the sequences and introduce the stoichiometry, the sequences of the unique chains are provided below:\n")
    unique_chains = get_unique_chains(pdb_models)
    for unique_chain in unique_chains:
        sys.stderr.write("\n" + unique_chain.id + ":\n")
        sys.stderr.write(unique_chain.get_sequence() + "\n")
    answer = str(input("\nDo you want to specify the stoichiometry? [yes/no]\n"))
    if answer == "yes":
        stoich_string = str(input("Please enter the desired stoichimetry:"))
        stoich_dict = get_stoichiometry(stoich_string)
        return stoich_dict
    else:
        return None


def data_transformation(pdb_models, verbose):
    """Transforms and aggregates the information of the PDB models in order to be more easy to get access to it.
    Takes the list of PDB models with the unified IDs according to sequence cimilarity and returns a dictionary where keys are the unique IDs that share some chains after getting
    unified by them and values are a list of tuples with information related to the chains that have that ID. The tuple consists on the model object that contains the chain with the specified ID and
    its interacting chain, the chain object that has the specified ID and the other chain object from the one (the one that interacts with the chain with the specified ID).
    Example:
        A => [(<Model id=0>, <Chain id=A>, <Chain id=A>), (<Model id=1>, <Chain id=A>, <Chain id=A>), (<Model id=2>, <Chain id=A>, <Chain id=B>)]
        B => [(<Model id=0>, <Chain id=B>, <Chain id=A>), (<Model id=2>, <Chain id=B>, <Chain id=B>)]

    Keyword arguments:
    pdb_models -- PDB models containing chains with unified IDs according to sequence similarity
    verbose -- boolean, prints to stderr the progress of the program"""
    interaction_dict = {}
    count = 0
    for pdb_model in pdb_models:
        for chain in pdb_model.get_chains():
            count += 1
            if count % 2 == 0:
                chain_index = 0
            else:
                chain_index = 1
            interaction_dict.setdefault(chain.id, [])
            if list(pdb_model.id) in [list(a[0].id) for a in interaction_dict[chain.id]]:
                continue
            interacting_chain = [other_chain for other_chain in pdb_model.get_chains()][chain_index]
            interaction_dict[chain.id].append((pdb_model, chain, interacting_chain))
    return interaction_dict



def starting_model(pdb_models, verbose):
    """Counts the number of times a specific unique chain is present in the set of PDB models. According to that frequency, scores the different models inferating
    which would be the model with more possible interactions, becoming the starting model of the complex in construction.
    If there is more than one model with maximum possible interactions, the starting model is randomly selected among the best ones.

    Keyword arguments:
    pdb_models -- PDB models containing chains with unified IDs according to sequence similarity
    verbose -- boolean, prints to stderr the progress of the program"""
    chain_count = {}
    possible_models = []
    for pdb_model in pdb_models:
        for chain in pdb_model.get_chains():
            if chain.id not in chain_count:
                chain_count[chain.id] = 0
            chain_count[chain.id] += 1
    models_interactions = {}
    for pdb_model in pdb_models:
        model_interactions = sum(chain_count[chain.id] for chain in pdb_model.child_list)
        models_interactions.setdefault(pdb_model, model_interactions)
    max_model_interaction = max(models_interactions.values())
    for model, model_interactions in models_interactions.items():
        if model_interactions == max_model_interaction:
            possible_models.append(model)
    if verbose:
        sys.stderr.write("The maximum number of possible interactions belongs to the model(s) {}\n".format([model.id for model in possible_models]))
    starting_model = random.choice(possible_models)
    return starting_model



def get_model_stoichiometry(model):
    """Generates a dictionary with the id chain as key and the number of repetitions of this chain as values. Returns the stoichiometry dictionary for the
    macrocomplex mmodel in construction.

    Keyword arguments
    model -- macrocomplex model that is in construction process"""
    stoichiometry = {}
    for chain in model:
        stoichiometry.setdefault(chain.id, 0)
        stoichiometry[chain.id] += 1
    return stoichiometry


def has_clashes(move_atoms, model, clash_dist):
    """Compares the backbone atoms of the moving chain with the backbone atoms of the model. Returns a boolean value based on the presence and abundancy of clashes.

    Keyword arguments:
    move_atoms -- set of atoms from the chain that is potentially going to be added to the macrocomplex
    model -- current macrocomplex model in construction
    clas_dist -- minimum clash distance between 2 atoms. The default minimum is 2 A

    Considerations:
    Alpha carbons are used as backbone atoms.
    Threshold: if 3% of atoms of the chain that is been comparing to the model shows clashes, <True> is returned. Hence, that chain won't be added to the model."""
    backbone = {"CA", "C1\'"}
    chain_atoms = [atom for atom in move_atoms if atom.id in backbone]
    model_atoms = [atom for atom in model.get_atoms() if atom.id in backbone]
    ns = PDB.NeighborSearch(model_atoms)
    clashes = 0
    for atom in chain_atoms:
        clashes += bool(ns.search(atom.coord, float(clash_dist)))
    if clashes/len(chain_atoms) >= 0.03:
        return True
    else:
        return False



def complex_builder(interaction_dict, pdb_models, num_models, max_chains, stoich_dict, clash_dist, verbose):
    """Function that iteratively builds a macrocomplex model from the pairs of interactions.
    First, it selects the starting pair of interacting chains of the macrocomplex. From that, the chains of the macrocomplex are being iterated looking for possible
    interactions according to the interaction dictionary. If it does not encounter clashes, the two chains that are equivalent superimpose, resulting in the addition of the
    other interacting chain to the complex. That specific model of interacting chains is removed from a list of pending interactions, since it is already added. If clashes are found,
    the model is not added, waiting. Returns a list of resulting macrocomplex models.

    Keyword arguments:
    interaction_dict -- dictionary of tuples with the possible interactions a given chain can perform.
    pdb_models -- PDB models containing chains with unified IDs according to sequence similarity
    num_models -- maximum number of models the function is going to produce
    max_chains -- maximum number of chains allowed in the models
    stoich_dict -- dictionary with the stoichiometry specified by the user
    clash_dist -- minimum clash distance between 2 atoms. The default minimum is 2 A
    verbose -- boolean, prints to stderr the progress of the program

    Considerations:
    Through all the process, the current stoichiometry of the building complex is analyzed. If a specific stoichiometry is provided, that is compared with the one of the macrocomplex.
    If the stoichiometry is fulfilled, the process stops.
    If there are no possible interactions left, the process stops.
    If the maximum number of chains is overpassed, the process stops.
    If no chains are added to complex during three iterations through the chains, the process stops"""
    output_objects = []
    for i in range(1, int(num_models) + 1):
        if verbose:
            sys.stderr.write("Building Macrocomplex " + str(i) + " ...\n")
        macrocomplex = starting_model(pdb_models, verbose)
        for key, tuple in interaction_dict.items():
            for tuple in interaction_dict[key]:
                if tuple[0].id == str(macrocomplex.id):
                    index = interaction_dict[key].index(tuple)
                    del interaction_dict[key][index]
        if verbose:
            sys.stderr.write("Building it from model {}, which contains chains {} and {}.\n\n".format(macrocomplex.id, [chain.id for chain in macrocomplex.get_chains()][0], [chain.id for chain in macrocomplex.get_chains()][1]))
        model_stoich = get_model_stoichiometry(macrocomplex)
        macrocomplex.id = "Model_" + str(i)
        run = True  # While this variable is true, the program will keep trying to add chains to the macrocomplex
        u = 1
        num_of_chains = 2
        num_empty_chains = 0
        while run:
            for chain in macrocomplex:    # Iterating through each chain of the macrocomplex
                if num_of_chains < int(max_chains):    # If the number of chains still hasn't reached the maximum allowed
                    interaction_copy = interaction_dict.copy()
                    if len(interaction_dict[chain.id]) != 0:   # Check if the chain has possible interactions left to make
                        random.shuffle(interaction_dict[chain.id])   # Shuffle the interactions list (avoiding repetitive behaviour)
                        for tuple in interaction_dict[chain.id]:
                            fix = tuple[1]    # Chain instance that corresponds to the same chain in macrocomplex
                            move = tuple[2]   # Chain instance that interacts with the other
                            if stoich_dict:
                                move_chain_id = move.id
                                model_stoich.setdefault(move_chain_id, 0)
                                model_number_chain = model_stoich[move_chain_id]
                                stoich_dict.setdefault(move_chain_id, 0)
                                if stoich_dict[move_chain_id] <= model_number_chain:  # Don't add the chain if surpasses the stoichiometry given
                                    if verbose:
                                        sys.stderr.write("  The current chain already fulfills the stoichiometry of the complex.\n")
                                        sys.stderr.write("  Chain " + move.id + " from model " + tuple[0].id + " NOT ADDED.\n\n")
                                    continue #   go to the next interaction tuple
                            sup = PDB.Superimposer()
                            chain_atoms, fix_atoms = chain.get_common_atoms(fix)
                            sup.set_atoms(chain_atoms, fix_atoms)  # Generate the superposition
                            sup.apply(move)
                            move_atoms = sorted(move.get_atoms())
                            if not has_clashes(move_atoms, macrocomplex, clash_dist):
                                if verbose:
                                    sys.stderr.write("  Succesful superposition between " + str(chain.id) + " from macrocomplex and chain " + fix.id + " from model " + tuple[0].id + ".\n")
                                    sys.stderr.write("  Chain " + str(num_of_chains) + " added: Chain " + move.id + " from model " + tuple[0].id + ".\n\n")
                                move.parent = None
                                macrocomplex.add(move)
                                model_stoich.setdefault(move.id, 0)
                                model_stoich[move.id] += 1
                                num_of_chains += 1
                                index = interaction_dict[chain.id].index(tuple)
                                del interaction_dict[chain.id][index]   # Deleting the model that has just been added from the interaction_dict
                                for redundant_tuple in interaction_dict[move.id]:
                                    if redundant_tuple[0].id == tuple[0].id:
                                        index = interaction_dict[move.id].index(redundant_tuple)
                                        del interaction_dict[move.id][index]
                            else:
                                if verbose:
                                    sys.stderr.write("  Unsuccesful superposition between " + str(chain.id) + " from macrocomplex and chain " + fix.id + " from model " + tuple[0].id + ".\n")
                                    sys.stderr.write("  Chain " + move.id + " from model " + tuple[0].id + " NOT ADDED.\n\n")
                    else:
                        if verbose:
                            sys.stderr.write("No possible interactions with chain " + chain.id + " left.\n")
                        num_empty_chains += 1
                    if stoich_dict:
                        if stoich_dict == model_stoich:
                            run = False
                            break
                else:
                    run = False  # When the maximum chain treshold is reached, stop running
                    break
            if num_empty_chains >= len(macrocomplex):  # If all chains of the macrocomplex have no possible interactions to make, stop running
                run = False
            if interaction_dict == interaction_copy:  # After 3 iterations without being able to add any chain to the macrocomplex, stop running
                u += 1
                if u == 3:
                    run = False
                    break

        stoichiometry_string = ""
        for key in sorted(model_stoich.keys()):
            stoichiometry_string += key + ":" + str(model_stoich[key]) + ","
        stoichiometry_string = stoichiometry_string[:-1]
        sys.stderr.write("\nStoichiometry of macrocomplex " + str(i) + " is: " + stoichiometry_string + ".\n")
        sys.stderr.write("Macrocomplex " + str(i) + " finished.\n")
        output_objects.append(macrocomplex)
    return output_objects




def save_results(out_models, output, directory, verbose):
    """Saves the resulting models into PDB files. Creates a specific directory for the model if it does not exist.
    Additionally, each chain receives a new ID in order to distinguish those chains that were equivalent.

    Keyword arguments:
    out_models -- list of the resulting model objects created by the program
    output -- name of the output model/file given by the user
    verbose -- boolean, prints to stderr the progress of the program"""
    u = 1
    if verbose:
        sys.stderr.write("Saving models...\n")
    io = PDB.PDBIO()
    if not os.path.exists(directory):
        os.makedirs(directory)
    for i in range(len(out_models)):
        id_list = []
        final_model = UpdModel(str(i))
        old_model = out_models[i]
        for chain in old_model.get_chains():
            new_chain = chain.copy()
            new_chain.id = new_id(id_list)
            id_list.append(new_chain.id)
            final_model.add(new_chain)
        io.set_structure(final_model)
        io.save(directory + "/" + output + "_" + str(u) + ".pdb")
        if verbose:
            sys.stderr.write("  " + output + "_" + str(u) + ".pdb saved\n")
        u += 1



def main(input, output, fasta, stoichiometry, model_num, chains, clash_dist, optimize, verbose):
    """Main function of the program that includes all the essential processes. It extracts the information from the input files in a list of PDB model objects.
    If stoichiometry is provided in the command line, it gets it directly from there. PDB model objects are updated, being unified the IDs of the chains according to their
    similarity (equivalence). If no FASTA input nor stoichiometry provided, it asks again for the stoichiometry in case the user wants to specify it. Information from model objects is
    reestructured in an interaction dictionary taking into account the equivalences between chains. "complex_builder" function builds possible commplexes from the information
    of the different pairs of interactions, specifying different types of parameters. Each resulting complex model is saved in PDB format in a specific directory with output name given.
    Finally, if the user chooses to optimize the models, they are refined with MODELLER tools and a DOPE plot comparison between the unoptimized and optimized model is made.

    Keyword arguments (all correspond to command line arguments)
    input -- input directory where PDB files with the pair interactions are located.
    output -- name of the output file. No extension is needed.
    fasta -- FASTA file with the sequences of the chains that will conform the macrocomplex.
    stoichiometry -- desired stoichiometry of the resulting complex.
    model_num -- maximum number of models the program is going to build
    chains -- maximum number of chains that the resulting model complex is going to have
    clash_dist -- minimum clash distance between 2 atoms. The default minimum is 2 A
    optimize -- boolean, refines the resulting model complexes and creates a DOPE plot comparison between the unoptimized and the optimized model
    verbose -- boolean, prints to stderr the progress of the program
    """
    sys.stderr.write("Initializing the program...\n")
    try:
        pdb_models = data_extraction(input, fasta, verbose)
    except (Directory_Not_Found, No_Input_Directory_Provided, No_PDB_Found, Incorrect_Number_Chains, FASTA_Not_Found) as e:
        print(e, file=sys.stderr)
        sys.exit(1)
    if stoichiometry:
        stoich_dict = get_stoichiometry(stoichiometry)
    else:
        stoich_dict = None
    updated_pdb_models = find_equivalent_chains(pdb_models, fasta, verbose)
    if fasta is None and stoichiometry is None:
        stoich_dict = get_stoich_without_fasta(updated_pdb_models)
    interaction_dict = data_transformation(updated_pdb_models, verbose)
    output_models = complex_builder(interaction_dict, pdb_models, model_num, chains, stoich_dict, clash_dist, verbose)
    directory = output + "_results"
    save_results(output_models, output, directory, verbose)
    if optimize:
        if verbose:
            sys.stderr.write("\nRefining the resulting complexes. This process can take a while.\n")
        path_list = []
        for pdb_file in os.listdir(directory):
            path = os.path.join(directory, pdb_file)
            path_list.append(path)
        for pdb_file in path_list:
            optimization(pdb_file)
        if verbose:
            sys.stderr.write("Refinement process finished.")
    sys.stderr.write("\nProgram finished.\n")
