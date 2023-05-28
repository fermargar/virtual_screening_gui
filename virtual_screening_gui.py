import tkinter as tk
import customtkinter
from tkinter import *
from tkinter import filedialog
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import MolToFile
import ipywidgets as widgets
import subprocess
import os
import sys

#Definition of appearance mode with customtkinter
customtkinter.set_appearance_mode("System")  # Modes: "System" (standard), "Dark", "Light"
customtkinter.set_default_color_theme("blue")  # Themes: "blue" (standard), "green", "dark-blue"

#Initialization of App with customtkinter.
class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()
        
        #Definition of window size and name of the app
        self.geometry("1280x600")
        self.title("2D/3D Virtual Screening")
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure((2, 3), weight=0)
        self.grid_rowconfigure((0, 1, 2), weight=1)
        
        #Definition of default content for dropdown menus in PharmScreen configuration
        minimization_var = customtkinter.StringVar(value="none")
        conformer_var = customtkinter.StringVar(value="No")
        charges_var = customtkinter.StringVar(value="none")
        logp_var = customtkinter.StringVar(value="none")
        minimization_var_library = customtkinter.StringVar(value="none")
        conformer_var_library = customtkinter.StringVar(value="No")
        charges_var_library = customtkinter.StringVar(value="none")
        logp_var_library = customtkinter.StringVar(value="none")
        alignment_method = customtkinter.StringVar(value="field")
        similarity_method = customtkinter.StringVar(value="tanimoto")
        field_1 = customtkinter.StringVar(value="hydroele")
        field_2 = customtkinter.StringVar(value="hydrocav")
        field_3 = customtkinter.StringVar(value="hbond")
        
        #It's not taking the value 8 by default. Value empty and get error till enter value in GUI field.
        threads_number = "8"
             
        # Creation of tab view for PharmScreen and RDKit
        self.program_selection = customtkinter.CTkTabview(self, width=500)
        self.program_selection.grid(row=0, column=0, padx=(10, 10), pady=(20, 0), sticky="nsew")
        self.program_selection.add("3D VS (PharmScreen)")
        self.program_selection.add("2D VS (RDKit)")
        
        # Creation of tab view for PharmScreen 
        self.tabview_pharmscreen = customtkinter.CTkTabview(self.program_selection.tab("3D VS (PharmScreen)"), width=800)
        self.tabview_pharmscreen.grid(row=2, column=0, padx=(10, 10), pady=(20, 0), sticky="nsew", columnspan=2)
        self.tabview_pharmscreen.add("Prepare Reference")
        self.tabview_pharmscreen.add("Prepare Library")
        self.tabview_pharmscreen.add("Virtual Screening")
        
        # Creation of tab view for RDKit
        self.tabview_rdkit = customtkinter.CTkTabview(self.program_selection.tab("2D VS (RDKit)"), width=500)
        self.tabview_rdkit.grid(row=0, column=0, padx=(10, 10), pady=(20, 0), sticky="nsew")
        self.tabview_rdkit.add("Run Similarity")
        
        # Creation of section for options for the GUI 
        self.tabview2 = customtkinter.CTkTabview(self, width=150)
        self.tabview2.grid(row=0, column=1, padx=(10, 10), pady=(20, 0), sticky="nsew")
        self.tabview2.add("GUI Options")
        
        # PHARMSCREEN FUNCTIONS
        # Open file explorer for library, reference and their prepared versions
        def open_reference_file_explorer():
            global reference, reference_path
            reference_path = filedialog.askopenfilename()
            reference = reference_path.split("/")[-1]
            label_reference.config(text=reference.split("/")[-1])
        
        def open_library_file_explorer():
            global library, library_path
            library_path = filedialog.askopenfilename()
            library = library_path.split("/")[-1]
            label_library.config(text=library.split("/")[-1])
        
        def open_prepared_reference_file_explorer():
            global prepared_reference, prepared_reference_path
            prepared_reference_path = filedialog.askopenfilename()
            prepared_reference = prepared_reference_path.split("/")[-1]
            label_prepared_reference.config(text=prepared_reference.split("/")[-1])
        
        def open_prepared_library_file_explorer():
            global prepared_library, prepared_library_path
            prepared_library_path = filedialog.askopenfilename()
            prepared_library = prepared_library_path.split("/")[-1]
            label_prepared_library.config(text=prepared_library.split("/")[-1])
            
        def open_license_file_explorer():
            global license_phs, license_phs_path
            license_phs_path = filedialog.askopenfilename()
            license_phs = license_phs_path.split("/")[-1]
            label_license_phs.config(text=license_phs.split("/")[-1])
            
        # Command-line options to run PharmScreen for reference and library preparation and virtual screening. There are different versions depending if the command line includes additional parameters to the default options (i.e. --tversky for virtual screening)
        def run_pharmscreen_reference():
            if conformer_var.get() == "No":
                print(threads_number.get())
                os.system(f"pharmscreen --key {license_phs_path} --input {reference_path} --minimize {minimization_var.get()} --single --name reference_prep --charges {charges_var.get()} --logp {logp_var.get()} --nameconf --odir {dir_ref_out.get()} --threads {threads_number.get()} &")
            else:
                os.system(f"pharmscreen --key {license_phs_path} --input {reference_path} --minimize {minimization_var.get()} --single --name reference_prep --charges {charges_var.get()} --logp {logp_var.get()} --genconf --nameconf --odir {dir_ref_out.get()} --threads {threads_number.get()} &")
        def run_pharmscreen_library():
            if conformer_var_library.get() == "No":
                os.system(f"pharmscreen --key {license_phs_path} --input {library_path} --minimize {minimization_var_library.get()} --single --name library_prep --charges {charges_var_library.get()} --logp {logp_var_library.get()} --nameconf --odir {dir_lib_out.get()} --threads {threads_number.get()} &")
            else:
                os.system(f"pharmscreen --key {license_phs_path} --input {library_path} --minimize {minimization_var_library.get()} --single --name library_prep --charges {charges_var_library.get()} --logp {logp_var_library.get()} --genconf --nameconf --odir {dir_lib_out.get()} --threads {threads_number.get()} &")
        def run_pharmscreen_vs():
            if similarity_method.get() == "tanimoto":
                os.system(f"pharmscreen --key {license_phs_path} --ref {prepared_reference_path} --input {prepared_library_path} --logp userdefined --nameconf --odir {dir_vs_out.get()} --align {alignment_method.get()} --field1 {field_1.get()} --field2 {field_2.get()} --field3 {field_3.get()} --threads {threads_number.get()} &")
            else:
                os.system(f"pharmscreen --key {license_phs_path} --ref {prepared_reference_path} --input {prepared_library_path} --logp userdefined --nameconf --tversky --odir {dir_vs_out.get()} --align {alignment_method.get()} --field1 {field_1.get()} --field2 {field_2.get()} --field3 {field_3.get()} --threads {threads_number.get()} &")

        
        # RDKIT FUNCTIONS
        # Open file explorer and open smile files.
        def open_rdkit_reference_file_explorer():
            global rdkit_reference, rdkit_reference_path
            rdkit_reference_path = filedialog.askopenfilename()
            rdkit_reference = rdkit_reference_path.split("/")[-1]
            label_rdkit_reference.config(text=rdkit_reference.split("/")[-1])        

        def open_rdkit_library_file_explorer():
            global rdkit_library, rdkit_library_path
            rdkit_library_path = filedialog.askopenfilename()
            rdkit_library = rdkit_library_path.split("/")[-1]
            label_rdkit_library.config(text=rdkit_library.split("/")[-1])        

        # Read the smiles with RDKit, save them in a list and run morgan fingerprint similarity. Similarity printed in terminal.
        def read_lines_smiles(output_file):
            molecules1 = [(Chem.MolFromSmiles(line.split(' ')[0]), line.split(' ')[1])
                        for line in open(rdkit_reference, 'r') if line.strip()] 
            with open(rdkit_library, 'r') as f2, open(output_file, 'w') as output:
                mol_counter = 0  # Counter for the number of molecules processed
                Draw.DrawingOptions.bondLineWidth = 1.2  # Adjust line width of bonds in the molecular drawing

                for line in f2:
                    line = line.strip()
                    if line:
                        smiles, name = line.split(' ')
                        molecule = Chem.MolFromSmiles(smiles)
                        if molecule is not None:
                            canonical_smiles = Chem.MolToSmiles(molecule) # Calculate the canonical smiles
                            morgan_fp = AllChem.GetMorganFingerprintAsBitVect(molecule, 4, nBits=2048) # Calculates Morgan Fingerprint 4
                            similarity_scores = [(Chem.DataStructs.TanimotoSimilarity(
                                morgan_fp, AllChem.GetMorganFingerprintAsBitVect(m1, 4, nBits=2048)), m1_name)
                                for m1, m1_name in molecules1]
                            similarity_scores.sort(key=lambda x: x[0],reverse=True)
                            results_printed = False
                            for score, m1_name in similarity_scores:
                                formatted_score = "{:.5f}".format(round(score, 5))
                                output.write(f"{name};{m1_name};{formatted_score}")
                                results_printed = True
                            if results_printed:
                                output.write('\n')
                            
            # Sort the lines in the output file based on the similarity score
            lines = open(output_file, 'r').readlines()
            lines = [line for line in lines if ':' in line]  # Exclude lines without a similarity score
            lines.sort(key=lambda x: float(x.split(':')[-1].strip()) if ':' in x else -1, reverse=True)
            with open(output_file, 'w') as output:
                output.writelines(lines)
                                
                
        def output_rdkit():
            output_file = "2D_similarity_ranking.csv"  # Specify the desired output file path
            read_lines_smiles(output_file)
                
        # Close window
        def close_pharmscreen():
            self.destroy()
        
        # Change appearance of the window
        def change_appearance_mode_event(new_appearance_mode: str):
            customtkinter.set_appearance_mode(new_appearance_mode)
        
        
        ######################################
        ########## PHARMSCREEN GUI ###########
        ######################################
        
        
        #Load library file and define the number of threads to run the experiments (preparation and virtual screening)
        button_license_phs = customtkinter.CTkButton(self.program_selection.tab("3D VS (PharmScreen)"), text="Select License", command=open_license_file_explorer, corner_radius=5, fg_color="#bcbcbc", text_color="black")
        button_license_phs.grid(padx=5, pady=5, row=0, column=0, sticky="W")
        label_license_phs = tk.Label(self.program_selection.tab("3D VS (PharmScreen)"), text="", bg="white", width=20, height=1, anchor="w")
        label_license_phs.grid(padx=5, pady=5, row=1, column=0, sticky="W")        
        
        threads_label = customtkinter.CTkLabel(self.program_selection.tab("3D VS (PharmScreen)"),text="Number of Threads")
        threads_label.grid(row=0, column=1, padx=5, pady=5, sticky="W")
        threads_number = customtkinter.CTkEntry(self.program_selection.tab("3D VS (PharmScreen)"))
        threads_number.grid(row=1, column=1, padx=5, pady=5, sticky="W")
 
        
        #Load reference, define parameters and run PharmScreen
        button_reference = customtkinter.CTkButton(self.tabview_pharmscreen.tab("Prepare Reference"), text="Select reference", command=open_reference_file_explorer, corner_radius=5, fg_color="#bcbcbc", text_color="black")
        button_reference.grid(pady=5, padx=5, row=0, column=0, sticky="W")
        label_reference = tk.Label(self.tabview_pharmscreen.tab("Prepare Reference"), text="", bg="white", width=20, height=1, anchor="w")
        label_reference.grid(padx=5, pady=5, row=0, column=1)
        
        minimization_label = customtkinter.CTkLabel(self.tabview_pharmscreen.tab("Prepare Reference"), text="Minimization method:")
        minimization_label.grid(row=1, column=0, padx=5, pady=5, sticky="W")
        minimization_dropdown = customtkinter.CTkOptionMenu(self.tabview_pharmscreen.tab("Prepare Reference"), variable=minimization_var, values=["none", "uff", "mmff", "am1", "rm1"])
        minimization_dropdown.grid(row=1, column=1, padx=5, pady=5, sticky="W")
        
        conformer_label = customtkinter.CTkLabel(self.tabview_pharmscreen.tab("Prepare Reference"), text="Calculate conformers:")
        conformer_label.grid(row=2, column=0, padx=5, pady=5, sticky="W")
        conformer_dropdown = customtkinter.CTkOptionMenu(self.tabview_pharmscreen.tab("Prepare Reference"), variable=conformer_var, values=["No", "Yes"])
        conformer_dropdown.grid(row=2, column=1, padx=5, pady=5, sticky="W")
        
        charges_label = customtkinter.CTkLabel(self.tabview_pharmscreen.tab("Prepare Reference"), text="Charges method:")
        charges_label.grid(row=3, column=0, padx=5, pady=5, sticky="W")
        charges_dropdown = customtkinter.CTkOptionMenu(self.tabview_pharmscreen.tab("Prepare Reference"), variable=charges_var, values=["none", "userdefined", "gasteiger", "bcc", "esp"])
        charges_dropdown.grid(row=3, column=1, padx=5, pady=5, sticky="W")
        
        logp_label = customtkinter.CTkLabel(self.tabview_pharmscreen.tab("Prepare Reference"), text="LogP method:")
        logp_label.grid(row=4, column=0, padx=5, pady=5, sticky="W")
        logp_dropdown = customtkinter.CTkOptionMenu(self.tabview_pharmscreen.tab("Prepare Reference"), variable=logp_var, values=["none", "userdefined", "at", "mst"])
        logp_dropdown.grid(row=4, column=1, padx=5, pady=5, sticky="W")
        
        dir_ref_out_label_library = customtkinter.CTkLabel(self.tabview_pharmscreen.tab("Prepare Reference"), text="Reference output folder")
        dir_ref_out_label_library.grid(row=5, column=0, padx=5, pady=5, sticky="W")
        dir_ref_out = customtkinter.CTkEntry(master=self.tabview_pharmscreen.tab("Prepare Reference"), placeholder_text="")
        dir_ref_out.grid(row=5, column=1, padx=5, pady=5, sticky="W")
        
        button_run_pharmscreen_reference = customtkinter.CTkButton(self.tabview_pharmscreen.tab("Prepare Reference"), text="Prepare reference", command=run_pharmscreen_reference, corner_radius=5, fg_color="#e69138", text_color="black")
        button_run_pharmscreen_reference.grid(row=6, column=0, padx=5, pady=5, sticky="W")
        
        #Load library, define parameters and run PharmScreen
        button_library = customtkinter.CTkButton(self.tabview_pharmscreen.tab("Prepare Library"), text="Select library", command=open_library_file_explorer, corner_radius=5, fg_color="#bcbcbc", text_color="black")
        button_library.grid(padx=5, pady=5, row=0, column=0, sticky="W")
        label_library = tk.Label(self.tabview_pharmscreen.tab("Prepare Library"), text="", bg="white", width=20, height=1, anchor="w")
        label_library.grid(padx=5, pady=5, row=0, column=1)
        
        minimization_label_library = customtkinter.CTkLabel(self.tabview_pharmscreen.tab("Prepare Library"), text="Minimization method:")
        minimization_label_library.grid(row=1, column=0, padx=5, pady=5, sticky="W")
        minimization_dropdown_library = customtkinter.CTkOptionMenu(self.tabview_pharmscreen.tab("Prepare Library"), variable=minimization_var_library, values=["none", "uff", "mmff", "am1", "rm1"])
        minimization_dropdown_library.grid(row=1, column=1, padx=5, pady=5, sticky="W")
        
        conformer_label_library = customtkinter.CTkLabel(self.tabview_pharmscreen.tab("Prepare Library"), text="Calculate conformers:")
        conformer_label_library.grid(row=2, column=0, padx=5, pady=5, sticky="W")
        conformer_dropdown_library = customtkinter.CTkOptionMenu(self.tabview_pharmscreen.tab("Prepare Library"), variable=conformer_var_library, values=["No", "Yes"])
        conformer_dropdown_library.grid(row=2, column=1, padx=5, pady=5, sticky="W")
        
        charges_label_library = customtkinter.CTkLabel(self.tabview_pharmscreen.tab("Prepare Library"), text="Charges method:")
        charges_label_library.grid(row=3, column=0, padx=5, pady=5, sticky="W")
        charges_dropdown_library = customtkinter.CTkOptionMenu(self.tabview_pharmscreen.tab("Prepare Library"), variable=charges_var_library, values=["none", "userdefined", "gasteiger", "bcc", "esp"])
        charges_dropdown_library.grid(row=3, column=1, padx=5, pady=5, sticky="W")
        
        logp_label_library = customtkinter.CTkLabel(self.tabview_pharmscreen.tab("Prepare Library"), text="LogP method:")
        logp_label_library.grid(row=4, column=0, padx=5, pady=5, sticky="W")
        logp_dropdown_library= customtkinter.CTkOptionMenu(self.tabview_pharmscreen.tab("Prepare Library"), variable=logp_var_library, values=["none", "userdefined", "at", "mst"])
        logp_dropdown_library.grid(row=4, column=1, padx=5, pady=5, sticky="W")
        
        dir_lib_out_label_library = customtkinter.CTkLabel(self.tabview_pharmscreen.tab("Prepare Library"), text="Library output folder")
        dir_lib_out_label_library.grid(row=5, column=0, padx=5, pady=5, sticky="W")
        dir_lib_out = customtkinter.CTkEntry(master=self.tabview_pharmscreen.tab("Prepare Library"), placeholder_text="")
        dir_lib_out.grid(row=5, column=1, padx=5, pady=5, sticky="W")
        
        button_run_pharmscreen_library = customtkinter.CTkButton(self.tabview_pharmscreen.tab("Prepare Library"), text="Prepare library", command=run_pharmscreen_library, corner_radius=5, fg_color="#e69138", text_color="black")
        button_run_pharmscreen_library.grid(row=6, column=0, padx=5, pady=5, sticky="W")
        
        
        #Load prepared reference and library, define parameters and run PharmScreen
        button_prepared_reference = customtkinter.CTkButton(self.tabview_pharmscreen.tab("Virtual Screening"), text="Select prepared reference", command=open_prepared_reference_file_explorer, corner_radius=5, fg_color="#bcbcbc", text_color="black")
        button_prepared_reference.grid(row=0, column=0, padx=5, pady=5, sticky="W")
        label_prepared_reference = tk.Label(self.tabview_pharmscreen.tab("Virtual Screening"), text="", bg="white", width=20, height=1, anchor="w")
        label_prepared_reference.grid(row=0, column=1, padx=5, pady=5, sticky="W")
        
        button_prepared_library = customtkinter.CTkButton(self.tabview_pharmscreen.tab("Virtual Screening"), text="Select prepared library", command=open_prepared_library_file_explorer, corner_radius=5, fg_color="#bcbcbc", text_color="black")
        button_prepared_library.grid(row=1, column=0, padx=5, pady=5, sticky="W")
        label_prepared_library = tk.Label(self.tabview_pharmscreen.tab("Virtual Screening"), text="", bg="white", width=20, height=1, anchor="w")
        label_prepared_library.grid(row=1, column=1, padx=5, pady=5, sticky="W")
       
        alignment_label = customtkinter.CTkLabel(self.tabview_pharmscreen.tab("Virtual Screening"), text="Alignment method:")
        alignment_label.grid(row=2, column=0, padx=5, pady=5, sticky="W")
        alignment_label_dropdown= customtkinter.CTkOptionMenu(self.tabview_pharmscreen.tab("Virtual Screening"), variable=alignment_method, values=["skip", "structure", "field", "montecarlo", "hq-montecarlo"])
        alignment_label_dropdown.grid(row=2, column=1, padx=5, pady=5, sticky="W")

        similarity_label = customtkinter.CTkLabel(self.tabview_pharmscreen.tab("Virtual Screening"), text="Similarity method:")
        similarity_label.grid(row=3, column=0, padx=5, pady=5, sticky="W")
        similarity_label_dropdown= customtkinter.CTkOptionMenu(self.tabview_pharmscreen.tab("Virtual Screening"), variable=similarity_method, values=["tanimoto", "tversky"])
        similarity_label_dropdown.grid(row=3, column=1, padx=5, pady=5, sticky="W")
        
        field_1_label = customtkinter.CTkLabel(self.tabview_pharmscreen.tab("Virtual Screening"), text="Field 1:")
        field_1_label.grid(row=4, column=0, padx=5, pady=5, sticky="W")
        field_1_label_dropdown= customtkinter.CTkOptionMenu(self.tabview_pharmscreen.tab("Virtual Screening"), variable=field_1, values=["electro", "shape", "hydro", "hydroele", "hydrocav", "hydrovdw", "hbond"])
        field_1_label_dropdown.grid(row=4, column=1, padx=5, pady=5, sticky="W")

        field_2_label = customtkinter.CTkLabel(self.tabview_pharmscreen.tab("Virtual Screening"), text="Field 2:")
        field_2_label.grid(row=5, column=0, padx=5, pady=5, sticky="W")
        field_2_label_dropdown= customtkinter.CTkOptionMenu(self.tabview_pharmscreen.tab("Virtual Screening"), variable=field_2, values=["electro", "shape", "hydro", "hydroele", "hydrocav", "hydrovdw", "hbond"])
        field_2_label_dropdown.grid(row=5, column=1, padx=5, pady=5, sticky="W")
        
        field_3_label = customtkinter.CTkLabel(self.tabview_pharmscreen.tab("Virtual Screening"), text="Field 3:")
        field_3_label.grid(row=6, column=0, padx=5, pady=5, sticky="W")
        field_3_label_dropdown= customtkinter.CTkOptionMenu(self.tabview_pharmscreen.tab("Virtual Screening"), variable=field_3, values=["electro", "shape", "hydro", "hydroele", "hydrocav", "hydrovdw", "hbond"])
        field_3_label_dropdown.grid(row=6, column=1, padx=5, pady=5, sticky="W")

        dir_vs_out_label_library = customtkinter.CTkLabel(self.tabview_pharmscreen.tab("Virtual Screening"), text="VS output folder")
        dir_vs_out_label_library.grid(row=7, column=0, padx=5, pady=5, sticky="W")
        dir_vs_out = customtkinter.CTkEntry(self.tabview_pharmscreen.tab("Virtual Screening"), placeholder_text="")
        dir_vs_out.grid(row=7, column=1, padx=5, pady=5, sticky="W")
        
        button_run_pharmscreen_vs = customtkinter.CTkButton(self.tabview_pharmscreen.tab("Virtual Screening"), text="Run Virtual Screening", command=run_pharmscreen_vs, corner_radius=5, fg_color="#e69138", text_color="black")
        button_run_pharmscreen_vs.grid(row=8, column=0, columnspan=5, padx=5, pady=5, sticky="W")
        
        
        # Theme selection button
        theme_mode = customtkinter.CTkLabel(self.tabview2.tab("GUI Options"), text="Theme")
        theme_mode.grid(row=0, column=0, padx=5, pady=5, sticky="W")
        appearance_mode_optionemenu = customtkinter.CTkOptionMenu(self.tabview2.tab("GUI Options"), values=["Light", "Dark", "System"], command=change_appearance_mode_event)
        appearance_mode_optionemenu.grid(row=0, column=1, padx=5, pady=5, sticky="W")
        
        #slider = widgets.IntSlider()
        #slider.grid(row=1, column=0, padx=5, pady=5, sticky="W")
        
        ############################################
        ################# RDKITGUI #################
        ############################################
        
        # Load two smile files (one or more smiles per file) and run tanimoto similarity with morgan fingerprints
        rdkit_button_reference = customtkinter.CTkButton(self.tabview_rdkit.tab("Run Similarity"), text="Select reference", command=open_rdkit_reference_file_explorer, corner_radius=5, fg_color="#bcbcbc", text_color="black")
        rdkit_button_reference.grid(pady=5, padx=5, row=0, column=0, sticky="W")
        label_rdkit_reference = tk.Label(self.tabview_rdkit.tab("Run Similarity"), text="", bg="white", width=20, height=1, anchor="w")
        label_rdkit_reference.grid(padx=5, pady=5, row=0, column=1)

        rdkit_button_library = customtkinter.CTkButton(self.tabview_rdkit.tab("Run Similarity"), text="Select library", command=open_rdkit_library_file_explorer, corner_radius=5, fg_color="#bcbcbc", text_color="black")
        rdkit_button_library.grid(pady=5, padx=5, row=1, column=0, sticky="W")
        label_rdkit_library = tk.Label(self.tabview_rdkit.tab("Run Similarity"), text="", bg="white", width=20, height=1, anchor="w")
        label_rdkit_library.grid(padx=5, pady=5, row=1, column=1)
        
        button_rdkit_prepare_reference = customtkinter.CTkButton(self.tabview_rdkit.tab("Run Similarity"), text="Run similarity", command=output_rdkit, corner_radius=5, fg_color="#e69138", text_color="black")
        button_rdkit_prepare_reference.grid(row=2, column=0, columnspan=5, padx=5, pady=5, sticky="W")


if __name__ == "__main__":
    app = App()
    app.mainloop()
