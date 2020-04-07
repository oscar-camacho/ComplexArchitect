#Idea of script for the GUI, being adapted from the example given by the Python teacher

import tkinter
import tkinter.filedialog
import tkinter.messagebox
import os,sys
from Bio import SeqIO
#from Bio.Alphabet import generic_dna, generic_protein
from Bio.Alphabet import IUPAC
import collections
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import re


#Initialize main window

class Application(tkinter.Frame):

    def select_PDB_directory(self):
        filetext="Select"
        directory = tkinter.filedialog.askdirectory(title = "Select the directory with PDB files")
        print (root.directory)
        bP["text"]= simpleName(root.directory) if root.directory else filetext
        return directory

    def select_fasta(self):
        fasta_filename = tkinter.filedialog.askopenfilename(title="Select a %s FASTA file",
                            filetypes=[("FASTA files","*.fasta"),("FASTA files","*.fa")] )
        return fasta_filename


    def quit(self):
        if tkinter.messagebox.askyesno("Quit","Are you sure you want to exit ComplexArchitect?"):
            tkinter.Frame.quit(self)

    def show_help(self):
        tkinter.messagebox.showinfo(title="Help", message="Here I should show the document with the help of the program")

    def show_about(self):
        tkinter.messagebox.showinfo(title="About", message="This program has been done in the PYT subject")


    def create_menu(self):
        self.menubar = tkinter.Menu(self)

        #CREATE THE FILEMENU
        filemenu = tkinter.Menu(self.menubar)
        filemenu.add_command(label="Select directory with PDB files", command=self.select_PDB_directory)
        filemenu.add_separator()
        filemenu.add_command(label="QUIT", command=self.quit)

        #CREATE THE HELP MENU
        helpmenu = tkinter.Menu(self.menubar)
        helpmenu.add_command(label="Help", command=None)
        helpmenu.add_command(label="About", command=None)
        self.menubar.add_cascade(label="File", menu=filemenu)
        self.menubar.add_cascade(label="Help", menu=helpmenu)

        self.master.config(menu=self.menubar)





    def open_FASTA(self):
        fasta_filename_path = tkinter.filedialog.askopenfilename(title="Select a FASTA file",
                            filetypes=[("FASTA files","*.fasta"),("FASTA files","*.fa")] )



    def create_sequence_textbox(self):
        text_frame = tkinter.LabelFrame(self, text="Sequences", width=400, padx=5, pady=5)
        self.sequence_text = tkinter.Text(text_frame)
        self.sequence_text.pack( )
        text_frame.grid( row=0, column=1 )













    def create_options_frame(self):

        self.options = tkinter.LabelFrame(self, text="Options")


        frame = tkinter.Frame(self.options)
        label_output = tkinter.Label(frame, text="Output name:")
        label_fasta = tkinter.Label(frame, text="FASTA file:")
        entry_fasta = tkinter.Button(frame, text="Select", command=self.open_FASTA)
        clear_fasta = tkinter.Button(frame, text="Clear", command=None)
        label_fasta_path = tkinter.Label(frame, textvariable=None)
        label_num_chains = tkinter.Label(frame, text="Max number of chains:")
        label_num_models = tkinter.Label(frame, text="Number of models:")
        label_current_dir = tkinter.Label(frame, text="Selected directory:")
        label_current_dir_path = tkinter.Label(frame, textvariable=None)
        self.label_optimize = tkinter.Checkbutton(frame, text="Optimize model", onvalue=True, offvalue=False, variable=None)
        entry_output = tkinter.Entry(frame, texstetvariable=None)
        self.entry_max_chains = tkinter.Spinbox(frame, from_=100, to=1000)
        self.entry_num_models = tkinter.Spinbox(frame, from_=1, to=100)
        self.entry_run = tkinter.Button(frame, text="Build model!", command=None)
        label_stoich = tkinter.Label(frame, text="Stoichiometry (optional):")
        entry_label_stoich = tkinter.Entry(frame, textvariable=None)

        #Option grids
        label_fasta.grid(row=1, column=0, sticky="w")
        entry_fasta.grid(row=1, column=1, sticky="w")
        label_fasta_path.grid(row=1, column=2)
        clear_fasta.grid(row=1, column=1, sticky="e")
        label_stoich.grid(row=5, column=0, sticky="w")
        entry_label_stoich.grid(row=5, column=1, columnspan=2, sticky="w")
        label_current_dir.grid(row=0, column=0, sticky="w")
        label_current_dir_path.grid(row=0, column=1, sticky="w", columnspan=3)
        label_output.grid(row=2, column=0, sticky="w")
        entry_output.grid(row=2, column=1, sticky="w")
        label_num_chains.grid(row=2, column=0, sticky="w")
        self.entry_max_chains.grid(row=2, column=1, sticky="w")
        label_num_models.grid(row=4, column=0, sticky="w")
        self.entry_num_models.grid(row=4, column=1, sticky="w")
        self.label_optimize.grid(row=1, column=4, sticky="w", columnspan=3)
        self.entry_run.grid(row=2, column=6, sticky="w")
        frame.grid(row=0)

        self.options.grid(row=0, column=0)



    def createWidgets(self):
        self.create_menu()
        self.create_options_frame()
        #self.create_sequence_textbox()
        #self.create_left_frame()
        #self.create_barplot_options_frame()

        #self.create_graphic_canvas()
        self.grid(row=0)

    def __init__(self, master=None, **kwargs):
        """Initalizates the app"""
        tkinter.Frame.__init__(self, master, **kwargs)

        self.master.wm_title("ComplexArchitect")
        self.master.resizable(width=True, height=True)
        logo = tkinter.Image("photo", file = os.path.dirname(os.path.realpath(__file__))+"/logo.png")
        self.master.iconphoto(False, logo)

        self.master.geometry("1000x500")

        #DEFINE ATTRIBUTES
        self.seq_list = []
        self.global_frequencies = {}
        self.seq_listbox = None

        self.menubar = None
        self.options = None
        self.sequence_text = None


        self.createWidgets()



root = tkinter.Tk()
app = Application(master=root, padx=10, pady=10)

app.mainloop()
