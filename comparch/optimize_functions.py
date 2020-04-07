#!/usr/bin/env python3
# In this module there are stored the functions necessary to opimize the resulting complex model and the energy plot.

from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients, molecular_dynamics, actions
import os,sys
import matplotlib.pyplot as plt
import pandas as pd


def blockPrint():
    """Function that blocks the printing to stdout in order to avoid all the additional data and processes of MODELLER."""
    sys.stdout = open(os.devnull, 'w')


def enablePrint():
    """Function that restores the printing to stdout."""
    sys.stdout = sys.__stdout__


def DOPE_comparison_plot(path, data_unopt, data_opt):
    """Creates a DOPE comparison plot between the unoptimized and the optimized model. Reads the energy profile data of each model from the corresponding files and saves
    the plot in <png> format in the current directory

    Keyword arguments
    path -- whole path of the resulting complex models with its proper name except the extension of the file
    data_unopt -- file with the energy profile data from the unoptimized model
    data_opt -- file with the energy profile data from the optimized model"""
    code = path.split('/')[-1]
    plt.clf()
    # Storing data from the unoptimized model
    data = pd.read_csv(data_unopt, sep="\s+", header=None)
    line_norm, = plt.plot(data[0], data[1], linewidth=1, label='Not Optimized')
    #Storing data from the optimized model
    data_o = pd.read_csv(data_opt, sep="\s+", header=None)
    line_opt, = plt.plot(data_o[0], data_o[1], linewidth=1, label='Optimized')
    plt.legend(handles=[line_norm, line_opt])

    plt.title(code + ' Energy Profile')
    plt.xlabel('Residue')
    plt.ylabel('DOPE_Energy')
    plt.axhline(linewidth=1, linestyle=':', color='r')
    #Creating the plot
    plt.savefig(path + '_EnergyProfile_plot.png', dpi=300)


def optimization(pdb_file):
    """Optimizes the stereochemistry of a given model saved in a pdb file, including non-bonded contacts. The function is derived from the MODELLER package.
    It also creates a plot with a comparison of the DOPE of the unoptimized model vs the optimized model.

    Keyword arguments
    pdb_file -- whole path of the resulting complex model in pdb format. It includes the directory where it is stored.

    Considerations:
    get_dope_profile() function: returns the individual residue components of the DOPE score, which can be used to detect poor regions of the model.
    get_smoothed() function: modifies the resulting energy profile. Each residue's energy is smoothed by weighted window averaging, using its own energy and the
                             energy of window residues either side of it; energies of residues are weighted by how close they are to the center residue in sequence.
    The model is refined using conjugate gradients."""
    blockPrint()
    env = environ()
    env.io.atom_files_directory = ['../atom_files'] # poner directorio donde iran las pdbs
    env.edat.dynamic_sphere = True
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')

    path, ext = pdb_file.split('.')
    dir, code = path.split('/')
    mdl = complete_pdb(env, pdb_file)

    atmsel = selection(mdl)   #Selecting all atoms

    # Calculating initial energy and energy profile for the unoptimized structure
    mpdf_ini = atmsel.energy()   # Initial energy
    dope_ini = mdl.assess_normalized_dope()  # Initial z-score
    mdl_ep_ini = atmsel.get_dope_profile()
    mdl_ep_ini_smoothed = mdl_ep_ini.get_smoothed(window=50)
    ep_txt_path = dir + '/' + code + '_DOPE_EnergyProfile.txt'
    mdl_ep_ini_smoothed.write_to_file(ep_txt_path)

    #Generate the restraints:
    mdl.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)
    # Create optimizer objects and set defaults for all further optimizations
    cg = conjugate_gradients(output='NO_REPORT')
    # Open a file to get basic stats on each optimization
    trcfil = open(dir + '/optimization_stats_' + code + '.txt', 'w')
    # Run CG on the all-atom selection
    cg.optimize(atmsel, max_iterations=20, actions=actions.trace(5, trcfil))

    # Finish off with some more CG
    cg.optimize(atmsel, max_iterations=20,actions=[actions.trace(5, trcfil)])

    # Calculate final energy and energy profile for the optimized structure
    mpdf_final = atmsel.energy()
    dope_final = mdl.assess_normalized_dope()

    enablePrint()
    print("\nUnoptimized model of " + code + ":")
    print("  Normalized DOPE z-score: " + str(dope_ini))

    print ("\nOptimized model of " + code + ":")
    print("  Normalized DOPE z-score: " + str(dope_final))
    #Writing out the optimized model in pdb format
    blockPrint()
    mdl.write(file = dir + '/' + code + '_optimized.pdb')

    mdl_ep_fin = atmsel.get_dope_profile()
    mdl_ep_fin_smoothed = mdl_ep_fin.get_smoothed(window=50)
    opt_ep_txt_path = path + '_optimized_DOPE_EnergyProfile.txt'
    mdl_ep_fin_smoothed.write_to_file(opt_ep_txt_path)
    DOPE_comparison_plot(path, ep_txt_path, opt_ep_txt_path)
