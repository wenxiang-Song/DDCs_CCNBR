# -*- coding: utf-8 -*-
import os
from io import BytesIO
import tkinter as tk
from tkinter import *
from tkinter import messagebox
from tkinter import ttk
import numpy as np
import pandas as pd
import pickle
import itertools
import xlsxwriter
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from e3fp.pipeline import fprints_from_mol
from e3fp.fingerprint.fprint import mean
import shutil
from PIL import Image, ImageTk
class WelcomeWindow:
    def __init__(self, master):
        self.master = master
        master.title('CCNBR')
        master.config(background='#FEFEFE')
        master.resizable(False, False)

        image = Image.open("./Data/page1.png")
        img = image.resize((1020, 600))
        self.my_img = ImageTk.PhotoImage(img)
        self.background_label = Label(master, image=self.my_img)

        self.background_label.grid(row=0, column=0, rowspan=4)

        self.enter_button1 = tk.Button(master, text='Predict Drug-drug Cocrystal', command=self.open_ddc_wind, font=('Arial', 15),
                                      bg='#FFFFFF', relief=tk.RAISED, width=25, borderwidth=4)
        self.enter_button1.grid(row=1, column=0, pady=20, rowspan=2, sticky=S)
        self.enter_button2 = tk.Button(master, text='Predict Cocrystal Coformer', command=self.open_ccf_wind, font=('Arial', 15),
                                      bg='#FFFFFF', relief=tk.RAISED, width=25, borderwidth=4)
        self.enter_button2.grid(row=2, column=0, pady=20)

    def open_ccf_wind(self):
        self.master.withdraw()
        root = tk.Toplevel(self.master)
        app = CCFpredict(root)
    def open_ddc_wind(self):
        self.master.withdraw()
        root = tk.Toplevel(self.master)
        app = DDCpredict(root)
class CCFpredict:
    def __init__(self, master):
        self.master = master
        master.title('Predict Cocrystal Coformer')
        master.config(background='#FFFFFF')
        master.resizable(False, False)
        self.model_var = tk.StringVar(value="Please select the model SMINBR or 3D-SMINBR.")
        self.input_smi_var = tk.StringVar(value="For example，The SMILES string for aspirin: CC(=O)OC1=CC=CC=C1C(=O)O")
        self.hit_num_var = tk.IntVar(value=20)

        image = Image.open("./Data/page2.png")
        img = image.resize((1025, 100))
        self.my_img = ImageTk.PhotoImage(img)
        self.background_label = Label(master, image=self.my_img)
        self.background_label.grid(row=0, column=0, columnspan=10)

        self.canvas = tk.Canvas(master, width=1020, height=320, bg="#FFFFFF")
        self.canvas.grid(row=1, column=0, rowspan=8, columnspan=10)

        self.lab_line1 = tk.Label(master, text="The recommendation model you can chose to perform prediction:",
                                  font=('Arial', 12, "bold"), fg='#233C6F', bg='#FFFFFF')
        self.lab_line1.grid(row=1, column=0, pady=5, columnspan=10)

        self.model_combobox = ttk.Combobox(master, textvariable=self.model_var, values=["SMINBR", "3D-SMINBR"],
                                           state='readonly', font=('Arial', 12), background='#FFFFFF', foreground='#000000', width=63)
        self.model_combobox.grid(row=2, column=0, pady=5, columnspan=10)

        enter_smiles = '''Please input a SMILES string for a molecule to predict:'''
        self.lab_line2 = tk.Label(master, text=enter_smiles, font=('Arial', 12, "bold"), fg='#233C6F', bg='#FFFFFF')
        self.lab_line2.grid(row=3, column=0, pady=5, columnspan=10)

        self.input_smi_entry = tk.Entry(master, textvariable=self.input_smi_var, font=('Arial', 12), bg='#FFFFFF', fg="#000000", width=65)
        self.input_smi_entry.grid(row=4, column=0, pady=5, columnspan=10)

        self.lab_line3 = tk.Label(master, text="The number of recommendations you want to get:",
                                  font=('Arial', 12, "bold"), fg='#233C6F', bg='#FFFFFF')
        self.lab_line3.grid(row=5, column=0, padx=10, pady=5, columnspan=10)

        self.hit_num_entry = tk.Entry(master, textvariable=self.hit_num_var, font=('Arial', 12), bg='#FFFFFF', fg="#000000" , width=65)
        self.hit_num_entry.grid(row=6, column=0, padx=10, pady=5, columnspan=10)

        self.run_button = tk.Button(master, text="Submit", command=self.run_recommendation, font=('Arial', 13, "bold"),
                                    bg='#FFFFFF', relief=tk.RAISED, width=10)
        self.run_button.grid(row=7, column=0, padx=10, pady=5, columnspan=10)



        title_sminbr = """
        
        Substructure-Molecular-Interaction Network-Based Recommendation (SMINBR)"""
        self.lab_long = tk.Label(master, text=title_sminbr,
                                  font=('Arial', 11, "bold"), fg='#000000', bg='#FFFFFF')
        self.lab_long.grid(row=9, column=0, columnspan=10)

        SMINBR_text = """
    We introduced chemical substructures into the molecular-interaction network and firstly constructed a substructure-molecular-interaction network. In this integrated network, a simple and efficient three-length path strategy was designed and validated. The mechanism of SMINBR could capture the structural complementarities that drive the formation of intermolecular interactions, thus yielding good performance. SMINBR could simultaneously utilize molecular-interaction relationships and structural information to perform multi-component crystal prediction. Moreover, with substructures to bridge new chemical entities (NCEs) and the existing molecular-interaction network, SMINBR can directly predict multi-component crystals for the NCEs.
        
    Reference: Lulu Zheng,# Bin Zhu,# Zengrui Wu, Fang Liang, Minghuang Hong, Guixia Liu, Weihua Li, Guobin Ren,* and Yun Tang*. SMINBR: An Integrated Network and Chemoinformatics Tool Specialized for Prediction of Two-Component Crystal Formation. J. Chem. Inf. Model. 2021, 61, 4290-4302.
        """
        self.lab_long = tk.Label(master, text=SMINBR_text,
                                  font=('Arial', 10), fg='#1C1C1C', bg='#FFFFFF', width=130, justify='left', wraplength=1000)
        self.lab_long.grid(row=10, column=0, columnspan=10)



        self.lab_long = tk.Label(master, text="      3D-Substructure-Molecular-Interaction Network-Based Recommendation (3D-SMINBR)",
                                  font=('Arial', 11, "bold"), fg='#000000', bg='#FFFFFF')
        self.lab_long.grid(row=11, column=0, columnspan=10)
        SMINBR3d_text = """
    Given that the 3D conformations of molecules play a significant role in the crystal formation and growth, the 2D chemical substructures in SMINBR were replaced by 3D chemical substructures and a 3D-substructure-molecular-interaction network was constructed. Moreover, the recommendation algorithm was modified and a edge weight parameter was introduced in 3D-SMINBR, which improves the model performance. 3D-SMINBR will be an effective tool to guide experimental cocrystal design.
        
    Reference: Lulu Zheng,# Bin Zhu,# Zengrui Wu, Fang Liang, Minghuang Hong, Guixia Liu, Weihua Li, Guobin Ren,* and Yun Tang*. Pharmaceutical Cocrystal Discovery via 3D-SMINBR: A New Network Recommendation Tool Augmented by 3D Molecular Conformations. J. Chem. Inf. Model. 2023, 63, 4301−4311.
        """
        self.lab_long = tk.Label(master, text=SMINBR3d_text,
                                  font=('Arial', 10), fg='#1C1C1C', bg='#FFFFFF', width=130, justify='left', wraplength=1000)
        self.lab_long.grid(row=12, column=0, columnspan=10)


        self.return_button = Button(master, text='Go to DDC Model', command=self.return_to_ddc, font=('Arial', 12),
                                    bg='#FFFFFF', relief=tk.RAISED)
        self.return_button.grid(row=2, column=9, pady=20)
    def return_to_ddc(self):
        self.master.withdraw()
        root = tk.Toplevel(self.master)
        app = DDCpredict(root)




    def run_recommendation(self):
        try:
            Model = self.model_var.get()
            input_smi = self.input_smi_var.get()
            if not input_smi:
                messagebox.showinfo("warning", "Please enter the molecule SMILES correctly !")
                return
            hit_num = self.hit_num_var.get()
            mol_interact_relations = pickle.load(open('./Data/data1_molecular_link_index.pkl', 'rb'))
            mol_index = pickle.load(open('./Data/data2_smiles_index.pkl', 'rb'))
            mol_syno = pickle.load(open('./Data/data3_synonym_smiles.pkl', 'rb'))
            can_mol_index = pickle.load(open('./Data/data4_smiles_index.pkl', 'rb'))
            mol_fps_2d = pickle.load(open('./Data/data6_fingerprint_2d.pkl', 'rb'))
            mol_fps_3d = pickle.load(open('./Data/data5_fingerprint_3d.pkl', 'rb'))
            can_mol_dict = dict(can_mol_index)
            input_mol = Chem.MolFromSmiles(input_smi)
            can_input_smi = Chem.MolToSmiles(input_mol, isomericSmiles=False)
            input_smi_index = can_mol_dict.get(can_input_smi, len(can_mol_index))

            print(len(mol_interact_relations))
            print(len(mol_index))

            fprint_params = {'bits': 1024,
                             'level': -1,
                             'first': 1,
                             'radius_multiplier': 1.189,
                             'stereo': True,
                             'counts': False,
                             'include_disconnected': True,
                             'rdkit_invariants': False,
                             'remove_duplicate_substructs': True,
                             'exclude_floating': True}
            if Model == 'SMINBR':
                substr_weight = 1.5
                mol_input_fps = mol_fps_2d
                if input_smi_index >= len(mol_index):
                    input_mol = Chem.MolFromSmiles(can_input_smi)
                    input_fp = [x for x in
                                AllChem.GetMorganFingerprintAsBitVect(input_mol, 2, nBits=1024, useChirality=False,
                                                                      useFeatures=True).ToBitString()]
                    mol_input_fps.append(input_fp)
                substructure_matrix = np.array(mol_input_fps, dtype=np.float64)
                substructure_matrix = substructure_matrix[:, np.sum(substructure_matrix, axis=0) != 0]
                mol_num, substructure_num = substructure_matrix.shape
                substructure_links = []
                for mol in range(mol_num):
                    for i in range(substructure_num):
                        if substructure_matrix[mol, i] == 1:
                            substructure_links.append([mol, mol_num + i])
                substructure_links = [item + [substr_weight] for item in substructure_links]
                mol_interact_relations = [item + [1 / substr_weight] for item in mol_interact_relations]
                links = mol_interact_relations + substructure_links
                mat_nodes = list(itertools.chain.from_iterable(links))
                mat_nodes = set(mat_nodes)
                mat_nodes = {np.int32(node) for node in mat_nodes}
                mat_size = np.int32(max(mat_nodes) + 1)
                network = np.zeros((mat_size, mat_size))
                for item in links:
                    network[np.int32(item[0]), np.int32(item[1])] = item[2]
                network = network + network.T
                score = np.dot(np.dot(network, network), network)
                degree = np.tile(np.sum(network, 0), (np.size(network, 0), 1))
                degree_multiply = degree * degree.T
                score = np.divide(score, 10 * degree_multiply)
                score[np.isnan(score)] = 0
                score[np.isinf(score)] = 0
                mol_index_df = pd.DataFrame(mol_index, columns=['Structure', 'Index'])
                mol_syno_df = pd.DataFrame(mol_syno, columns=['Synonym', 'Structure'])
                mol_index_df = pd.merge(mol_syno_df, mol_index_df, on='Structure', how='inner')
                target_score = score[:, input_smi_index]
                target_score_df = pd.DataFrame({'Score': target_score})
                target_score_df['Index'] = range(np.size(target_score, axis=0))
                target_score_df = pd.merge(target_score_df, mol_index_df, on='Index', how='inner')
                target_score_df = target_score_df.sort_values(by='Score', ascending=False)
                target_score_df['Rank'] = target_score_df['Score'].rank(ascending=False, method='dense')
                target_score_df['Rank'] = target_score_df['Rank'].astype(int)
                target_hits_df = target_score_df[target_score_df['Rank'] <= np.int32(hit_num)]
                if input_smi_index < mol_index_df.shape[0]:
                    target_hits_df = target_hits_df.drop(columns=['Score', 'Index'],
                                                         axis=1)
                    target_ligand_index = np.where(
                        network[:mol_index_df.shape[0], input_smi_index] == 1 / substr_weight)
                    target_ligand = mol_index_df.iloc[
                        target_ligand_index].copy()
                    target_ligand['Rank'] = 'Known'
                    target_ligand = target_ligand.drop(columns='Index',
                                                       axis=1)
                    target_hits_df = target_hits_df[~target_hits_df['Structure'].isin(
                        target_ligand['Structure'])]
                    target_hits_df = pd.concat([target_ligand, target_hits_df], join='outer')
                else:
                    target_hits_df = target_hits_df.drop(columns=['Score', 'Index'], axis=1)
                target_hits_df.to_csv('./results/SMINBR.csv', index=False)
            elif Model == '3D-SMINBR':
                substr_weight = 0.12
                mol_input_fps = mol_fps_3d
                if input_smi_index >= len(mol_index):
                    input_mol = Chem.MolFromSmiles(can_input_smi)
                    input_mol.SetProp('_Name', str(input_smi_index))
                    Chem.AddHs(input_mol)
                    AllChem.EmbedMolecule(input_mol, randomSeed=42, useExpTorsionAnglePrefs=True,
                                          useBasicKnowledge=True)
                    prop = AllChem.MMFFGetMoleculeProperties(input_mol)
                    mmff = AllChem.MMFFGetMoleculeForceField(input_mol, prop)
                    mmff.Minimize()
                    fprint_params = fprint_params
                    input_fp = fprints_from_mol(input_mol, fprint_params=fprint_params)
                    input_fp = mean(input_fp)
                    input_fp = input_fp.to_vector(sparse=0)
                    mol_input_fps.append(input_fp)
                substructure_matrix = np.array(mol_input_fps, dtype=np.float64)
                substructure_matrix = substructure_matrix[:, np.sum(substructure_matrix, axis=0) != 0]
                mol_num, substructure_num = substructure_matrix.shape
                substructure_links = []
                for mol in range(mol_num):
                    for i in range(substructure_num):
                        if substructure_matrix[mol, i] == 1:
                            substructure_links.append([mol, mol_num + i])
                substructure_links = [item + [substr_weight] for item in substructure_links]
                mol_interact_relations = [item + [1 / substr_weight] for item in mol_interact_relations]
                links = mol_interact_relations + substructure_links
                mat_nodes = list(itertools.chain.from_iterable(links))
                mat_nodes = set(mat_nodes)
                mat_nodes = {np.int32(node) for node in mat_nodes}
                mat_size = np.int32(max(mat_nodes) + 1)
                network = np.zeros((mat_size, mat_size))
                for item in links:
                    network[np.int32(item[0]), np.int32(item[1])] = item[2]
                network = network + network.T
                degree = np.tile(np.sum(network, 0), (np.size(network, 0), 1))
                degree = degree.T
                trans_mat = np.divide(network, degree)
                unit_mat = np.eye(np.size(trans_mat, 0))
                score = unit_mat
                step_i = 0
                while step_i < 3:
                    score = np.dot(trans_mat.T, score)
                    step_i += 1
                score = score + score.T
                score[np.isnan(score)] = 0
                score[np.isinf(score)] = 0
                mol_index_df = pd.DataFrame(mol_index, columns=['Structure', 'Index'])
                mol_syno_df = pd.DataFrame(mol_syno, columns=['Synonym', 'Structure'])
                mol_index_df = pd.merge(mol_syno_df, mol_index_df, on='Structure',
                                        how='inner')
                target_score = score[:, input_smi_index]
                target_score_df = pd.DataFrame({'Score': target_score})
                target_score_df['Index'] = range(np.size(target_score, axis=0))
                target_score_df = pd.merge(target_score_df, mol_index_df, on='Index',
                                           how='inner')
                target_score_df = target_score_df.sort_values(by='Score',
                                                              ascending=False)
                target_score_df['Rank'] = target_score_df['Score'].rank(ascending=False,
                                                                        method='min')
                target_score_df['Rank'] = target_score_df['Rank'].astype(
                    int)
                target_hits_df = target_score_df[
                    target_score_df['Rank'] <= np.int32(hit_num)]
                if input_smi_index < mol_index_df.shape[0]:
                    target_hits_df = target_hits_df.drop(columns=['Score', 'Index'],
                                                         axis=1)
                    target_ligand_index = np.where(
                        network[:mol_index_df.shape[0], input_smi_index] == 1 / substr_weight)
                    target_ligand = mol_index_df.iloc[
                        target_ligand_index].copy()
                    target_ligand['Rank'] = 'Known'
                    target_ligand = target_ligand.drop(columns='Index',
                                                       axis=1)
                    target_hits_df = target_hits_df[
                        ~target_hits_df['Structure'].isin(
                            target_ligand['Structure'])]
                    target_hits_df = pd.concat([target_ligand, target_hits_df], join='outer')
                else:
                    target_hits_df = target_hits_df.drop(columns=['Score', 'Index'], axis=1)
                target_hits_df.to_csv('./results/3D-SMINBR.csv', index=False)
            file_name = Model
            def mol_to_excel(file):
                header = ['Synonym', 'SMILES', 'Rank']
                item_style = {
                    'align': 'center',
                    'valign': 'vcenter',
                    'top': 2,
                    'left': 2,
                    'right': 2,
                    'bottom': 2,
                    'text_wrap': 1
                }
                header_style = {
                    'bold': 1,
                    'valign': 'vcenter',
                    'align': 'center',
                    'top': 2,
                    'left': 2,
                    'right': 2,
                    'bottom': 2
                }
                excel_file_path = os.path.join("results", f'{file}.xlsx')
                workbook = xlsxwriter.Workbook(excel_file_path)
                ItemStyle = workbook.add_format(item_style)
                HeaderStyle = workbook.add_format(header_style)
                worksheet = workbook.add_worksheet()
                worksheet.set_column('A:A', 38)
                worksheet.set_column('B:B', 40)
                worksheet.set_column('C:C', 10)
                for ix_, i in enumerate(header):
                    worksheet.write(0, ix_, i, HeaderStyle)
                df = pd.read_csv(os.path.join("results", f"{file}.csv"))
                for i in range(df.shape[0]):
                    synonym = df.iloc[i, 0]
                    structure_smi = df.iloc[i, 1]
                    rank = df.iloc[i, 2]
                    img_data_structure = BytesIO()
                    c_structure = Chem.MolFromSmiles(structure_smi)
                    img_structure = Draw.MolToImage(c_structure)
                    img_structure.save(img_data_structure, format='PNG')
                    worksheet.set_row(i + 1, 185)
                    worksheet.insert_image(i + 1, 0, 'f',
                                           {'x_scale': 0.9, 'y_scale': 0.8, 'image_data': img_data_structure,
                                            'positioning': 1})
                    worksheet.write(i + 1, 1, structure_smi, ItemStyle)
                    worksheet.write(i + 1, 2, rank, ItemStyle)
                workbook.close()
            mol_to_excel(file_name)
            filename = f"{Model}.xlsx"
            source_path = os.path.abspath(os.path.join("results", filename))
            desktop_path = os.path.expanduser("~/Desktop")
            destination_path = os.path.join(desktop_path, filename)
            try:
                shutil.copy2(source_path, destination_path)
                messagebox.showinfo("Success",
                                    f"Recommendation process complete. Excel file saved on desktop: {destination_path}")
            except Exception as e:
                messagebox.showinfo("Error", f"Error: {e}")
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {str(e)}")
class DDCpredict:
    def __init__(self, master):
        self.master = master
        master.title('Predict Drug-drug Cocrystal')
        master.config(background='#FFFFFF')
        master.resizable(False, False)
        self.model_var = tk.StringVar(value="Please select the model TDDCNBR or NTDDCNBR")
        self.input_smi_var = tk.StringVar(value="")
        self.DDC_num_var = tk.IntVar(value=50)

        image = Image.open("./Data/page3.png")
        img = image.resize((1025, 100))
        self.my_img = ImageTk.PhotoImage(img)
        self.background_label = Label(master, image=self.my_img)
        self.background_label.grid(row=0, column=0, columnspan=50)

        self.canvas = tk.Canvas(master, width=1020, height=420, bg="#FEFEFE")
        self.canvas.grid(row=1, column=0, rowspan=8, columnspan=50)

        self.lab_line1 = tk.Label(master, text="The recommendation model you can chose to perform prediction:",
                                  font=('Arial', 12, "bold"), fg='#233C6F', bg='#FFFFFF')
        self.lab_line1.grid(row=1, column=0, pady=5, columnspan=50)


        self.model_combobox = ttk.Combobox(master, textvariable=self.model_var, values=["TDDCNBR", "NTDDCNBR"],
                                           state='readonly', font=('Arial', 12), background='#FFFFFF',
                                           foreground='#000000', width=68)
        self.model_combobox.grid(row=2, column=0, pady=5, columnspan=50)


        enter_smiles = '''TDDCNBR: Please input a SMILES string for a drug to predict:
NTDDCNBR: Please input SMILES strings for some drugs to predict:'''
        self.lab_line2 = tk.Label(master, text=enter_smiles, font=('Arial', 12, "bold"), fg='#233C6F', bg='#FFFFFF')
        self.lab_line2.grid(row=3, column=0, pady=5, columnspan=50)

        self.input_smi_entry = tk.Text(master, height=8, font=('Arial', 12), wrap="none", bg='#FFFFFF', fg="#000000", width=70)
        self.input_smi_entry.grid(row=4, column=0, pady=5, columnspan=50)

        self.yscrollbar = tk.Scrollbar(master, orient='vertical')
        self.yscrollbar.grid(row=4, column=46, sticky='ens')
        self.yscrollbar.config(command=self.input_smi_entry.yview)
        self.xscrollbar = tk.Scrollbar(master, orient='horizontal')
        self.xscrollbar.grid(row=4, column=11, sticky='ews', columnspan=36)
        self.xscrollbar.config(command=self.input_smi_entry.xview)
        self.input_smi_entry.config(yscrollcommand=self.yscrollbar.set)
        self.input_smi_entry.config(xscrollcommand=self.xscrollbar.set)

        self.lab_line3 = tk.Label(master, text="The number of recommendations you want to get:",
                                  font=('Arial', 12, "bold"), fg='#233C6F', bg='#FFFFFF')
        self.lab_line3.grid(row=6, column=0, pady=5, columnspan=50)

        self.DDC_num_entry = tk.Entry(master, textvariable=self.DDC_num_var, font=('Arial', 12), bg='#FFFFFF', fg="#000000", width=70)
        self.DDC_num_entry.grid(row=7, column=0, pady=5, columnspan=50)

        self.run_button = tk.Button(master, text="Submit", command=self.run_recommendation, font=('Arial', 13, "bold"),
                                    bg='#FFFFFF', relief=tk.RAISED, width=20)
        self.run_button.grid(row=8, column=0, padx=10, pady=5, columnspan=50)



        title_sminbr = """

        Targeted Drug-drug Cocrystal Network-Based Recommendation (TDDCNBR)"""
        self.lab_long = tk.Label(master, text=title_sminbr,
                                 font=('Arial', 11, "bold"), fg='#000000', bg='#FFFFFF')
        self.lab_long.grid(row=9, column=0, columnspan=50)

        SMINBR_text = """
    TDDCNBR and NTDDCNBR are both built on the existing 3D-Substructure-Molecular-Interaction Network. TDDCNBR filters the predicted results of 3D-SMINBR, retaining only the drug coformers to recommend drug-drug cocrystals. The advantage of the model is that it can purposefully recommend drug-drug cocrystals for a drug, the drawback is that the recommended drug coformers are limited to the drug molecules within the network. Furthermore, the calculation only considers the binding between the input molecule and other molecules in the network, without taking into account whether other drug molecules outside the network could potentially form cocrystals. This model is suitable for targeted drug-drug cocrystals  recommendations, such as finding drug coformers that can form cocrystal with aspirin.
        """
        self.lab_long = tk.Label(master, text=SMINBR_text,
                                 font=('Arial', 10), fg='#1C1C1C', bg='#FFFFFF', width=130, justify='left',
                                 wraplength=1000)
        self.lab_long.grid(row=10, column=0, columnspan=50)

        self.lab_long = tk.Label(master,
                                 text="Non Targeted Drug-drug Cocrystal Network-Based Recommendation (NTDDCNBR)",
                                 font=('Arial', 11, "bold"), fg='#000000', bg='#FFFFFF')
        self.lab_long.grid(row=11, column=0, columnspan=50)
        SMINBR3d_text = """
    NTDDCNBR  replaces the original single input channel of 3D-SMINBR with multiple input channels. All input molecules must be drug molecules, allowing for both internal and external drug molecules to be considered. The model computes the number and score of third-order pathways between the input drug molecules, thus suggesting the molecules most probable to form drug-drug co-crystals. This model is ideal for large-scale screening or focusing on specific types of drug-drug cocrystals discovery. For instance, to find antihypertensive drug-drug cocrystals other than sacubitril-valsartan, inputting various antihypertensive drugs into the network model would enable the model to identify combinations of drugs with higher potential for cocrystal formation.
        """
        self.lab_long = tk.Label(master, text=SMINBR3d_text,
                                 font=('Arial', 10), fg='#1C1C1C', bg='#FFFFFF', width=130, justify='left',
                                 wraplength=1000)
        self.lab_long.grid(row=12, column=0, columnspan=50)

        self.return_button = Button(master, text='GO to CCF Model', command=self.return_to_ccf, font=('Arial', 12),
                                    bg='#FFFFFF', relief=tk.RAISED)
        self.return_button.grid(row=2, column=47, pady=20, sticky=W)
    def return_to_ccf(self):
        self.master.withdraw()
        root = tk.Toplevel(self.master)
        app = CCFpredict(root)

    def run_recommendation(self):
        try:
            Model = self.model_var.get()
            DDC_num = self.DDC_num_var.get()
            mol_interact_relations = pickle.load(open('./Data/data1_molecular_link_index.pkl', 'rb'))
            mol_index = pickle.load(open('./Data/data2_smiles_index.pkl', 'rb'))
            mol_syno = pickle.load(open('./Data/data3_synonym_smiles.pkl', 'rb'))
            can_mol_index = pickle.load(open('./Data/data4_smiles_index.pkl', 'rb'))
            mol_fps = pickle.load(open('./Data/data5_fingerprint_3d.pkl', 'rb'))
            hit_num = len(mol_fps)
            substr_weight = 0.12
            fprint_params = {'bits': 1024,
                             'level': -1,
                             'first': 1,
                             'radius_multiplier': 1.189,
                             'stereo': True,
                             'counts': False,
                             'include_disconnected': True,
                             'rdkit_invariants': False,
                             'remove_duplicate_substructs': True,
                             'exclude_floating': True}
            item_style = {
                'align': 'center',
                'valign': 'vcenter',
                'top': 2,
                'left': 2,
                'right': 2,
                'bottom': 2,
                'text_wrap': 1
            }
            header_style = {
                'bold': 1,
                'valign': 'vcenter',
                'align': 'center',
                'top': 2,
                'left': 2,
                'right': 2,
                'bottom': 2
            }
            if Model == 'TDDCNBR':
                input_smi = self.input_smi_entry.get("1.0", tk.END).splitlines()
                if len(input_smi) > 1:
                    messagebox.showinfo("warning", "Only SMILES of one molecule can be entered !")
                    return
                input_smi = input_smi[0]
                can_mol_dict = dict(can_mol_index)
                input_mol = Chem.MolFromSmiles(input_smi)
                can_input_smi = Chem.MolToSmiles(input_mol, isomericSmiles=False)
                input_smi_index = can_mol_dict.get(can_input_smi, len(can_mol_index))
                mol_input_fps = mol_fps
                if input_smi_index >= len(mol_index):
                    input_mol = Chem.MolFromSmiles(can_input_smi)
                    input_mol.SetProp('_Name', str(input_smi_index))
                    Chem.AddHs(input_mol)
                    AllChem.EmbedMolecule(input_mol, randomSeed=42, useExpTorsionAnglePrefs=True,
                                          useBasicKnowledge=True)
                    prop = AllChem.MMFFGetMoleculeProperties(input_mol)
                    mmff = AllChem.MMFFGetMoleculeForceField(input_mol, prop)
                    mmff.Minimize()
                    fprint_params = fprint_params
                    input_fp = fprints_from_mol(input_mol, fprint_params=fprint_params)
                    input_fp = mean(input_fp)
                    input_fp = input_fp.to_vector(sparse=0)
                    mol_input_fps.append(input_fp)
                substructure_matrix = np.array(mol_input_fps, dtype=np.float64)
                substructure_matrix = substructure_matrix[:, np.sum(substructure_matrix, axis=0) != 0]
                mol_num, substructure_num = substructure_matrix.shape
                substructure_links = []
                for mol in range(mol_num):
                    for i in range(substructure_num):
                        if substructure_matrix[mol, i] == 1:
                            substructure_links.append([mol, mol_num + i])
                substructure_links = [item + [substr_weight] for item in substructure_links]
                mol_interact_relations = [item + [1 / substr_weight] for item in mol_interact_relations]
                links = mol_interact_relations + substructure_links
                mat_nodes = list(itertools.chain.from_iterable(links))
                mat_nodes = set(mat_nodes)
                mat_nodes = {np.int32(node) for node in mat_nodes}
                mat_size = np.int32(max(mat_nodes) + 1)
                network = np.zeros((mat_size, mat_size))
                for item in links:
                    network[np.int32(item[0]), np.int32(item[1])] = item[2]
                network = network + network.T
                degree = np.tile(np.sum(network, 0), (np.size(network, 0), 1))
                degree = degree.T
                trans_mat = np.divide(network, degree)
                unit_mat = np.eye(np.size(trans_mat, 0))
                score = unit_mat
                step_i = 0
                while step_i < 3:
                    score = np.dot(trans_mat.T, score)
                    step_i += 1
                score = score + score.T
                score[np.isnan(score)] = 0
                score[np.isinf(score)] = 0
                mol_index_df = pd.DataFrame(mol_index, columns=['Structure', 'Index'])
                mol_syno_df = pd.DataFrame(mol_syno, columns=['Synonym', 'Structure'])
                mol_index_df = pd.merge(mol_syno_df, mol_index_df, on='Structure',
                                        how='inner')
                target_score = score[:, input_smi_index]
                target_score_df = pd.DataFrame({'Score': target_score})
                target_score_df['Index'] = range(np.size(target_score, axis=0))
                target_score_df = pd.merge(target_score_df, mol_index_df, on='Index',
                                           how='inner')
                target_score_df = target_score_df.sort_values(by='Score',
                                                              ascending=False)
                target_score_df['Rank'] = target_score_df['Score'].rank(ascending=False,
                                                                        method='min')
                target_score_df['Rank'] = target_score_df['Rank'].astype(int)
                target_hits_df = target_score_df[
                    target_score_df['Rank'] <= np.int32(hit_num)]

                if input_smi_index < mol_index_df.shape[0]:
                    target_hits_df = target_hits_df.drop(columns=['Score', 'Index'],
                                                         axis=1)
                    target_ligand_index = np.where(
                        network[:mol_index_df.shape[0], input_smi_index] == 1 / substr_weight)
                    target_ligand = mol_index_df.iloc[
                        target_ligand_index].copy()
                    target_ligand['Rank'] = 'Known'
                    target_ligand = target_ligand.drop(columns='Index',
                                                       axis=1)
                    target_hits_df = target_hits_df[
                        ~target_hits_df['Structure'].isin(
                            target_ligand['Structure'])]
                    target_hits_df = pd.concat([target_ligand, target_hits_df], join='outer')
                else:
                    target_hits_df = target_hits_df.drop(columns=['Score', 'Index'], axis=1)
                reference_file = pd.read_csv('Data/data_drug.csv', encoding='gbk')
                reference_values = set(reference_file['SMILES'])
                filtered_df = target_hits_df[target_hits_df['Structure'].isin(reference_values)]
                filtered_df.to_csv('./results/TDDCNBR.csv', index=False, encoding='gbk')
                file_name = 'TDDCNBR'
                def mol_to_excel(file):
                    header = ['Synonym', 'SMILES', 'Rank', 'Drug']
                    excel_file_path = os.path.join("results", f'{file}.xlsx')
                    workbook = xlsxwriter.Workbook(excel_file_path)
                    ItemStyle = workbook.add_format(item_style)
                    HeaderStyle = workbook.add_format(header_style)
                    worksheet = workbook.add_worksheet()
                    worksheet.set_column('A:A', 38)
                    worksheet.set_column('B:B', 40)
                    worksheet.set_column('C:C', 10)
                    for ix_, i in enumerate(header):
                        worksheet.write(0, ix_, i, HeaderStyle)
                    df = pd.read_csv(os.path.join("results", f"{file}.csv"), encoding='gbk').head(DDC_num)
                    df['Drug'] = df.iloc[:, 0]
                    for i in range(df.shape[0]):
                        synonym = df.iloc[i, 0]
                        structure_smi = df.iloc[i, 1]
                        rank = df.iloc[i, 2]
                        img_data_structure = BytesIO()
                        c_structure = Chem.MolFromSmiles(structure_smi)
                        img_structure = Draw.MolToImage(c_structure)
                        img_structure.save(img_data_structure, format='PNG')
                        worksheet.set_row(i + 1, 185)
                        worksheet.insert_image(i + 1, 0, 'f',
                                               {'x_scale': 0.9, 'y_scale': 0.8, 'image_data': img_data_structure,
                                                'positioning': 1})
                        worksheet.write(i + 1, 1, structure_smi, ItemStyle)
                        worksheet.write(i + 1, 2, rank, ItemStyle)
                        worksheet.write(i + 1, 3, synonym, ItemStyle)
                    workbook.close()
                mol_to_excel(file_name)
            elif Model == 'NTDDCNBR':
                mol_interact_relations_copy = mol_interact_relations.copy()
                input_all_smi = self.input_smi_entry.get("1.0", tk.END).splitlines()
                if len(input_all_smi) <= 1:
                    messagebox.showinfo("warning", "Please enter the SMILES of two or more molecules !")
                    return
                can_mol_dict = dict(can_mol_index)
                mol_input_fps = mol_fps
                molecules_index = []
                n = len(can_mol_index)
                for input_smi in input_all_smi:
                    input_mol = Chem.MolFromSmiles(input_smi)
                    can_input_smi = Chem.MolToSmiles(input_mol, isomericSmiles=False)
                    input_smi_index = can_mol_dict.get(can_input_smi, n)
                    if input_smi_index == n:
                        n += 1
                    molecules_index.append(input_smi_index)
                    if input_smi_index >= len(mol_index):
                        input_mol = Chem.MolFromSmiles(can_input_smi)
                        input_mol.SetProp('_Name', str(input_smi_index))
                        Chem.AddHs(input_mol)
                        AllChem.EmbedMolecule(input_mol, randomSeed=42, useExpTorsionAnglePrefs=True,
                                              useBasicKnowledge=True)
                        prop = AllChem.MMFFGetMoleculeProperties(input_mol)
                        mmff = AllChem.MMFFGetMoleculeForceField(input_mol, prop)
                        mmff.Minimize()
                        fprint_params = fprint_params
                        input_fp = fprints_from_mol(input_mol, fprint_params=fprint_params)
                        input_fp = mean(input_fp)
                        input_fp = input_fp.to_vector(sparse=0)
                        mol_input_fps.append(input_fp)
                substructure_matrix = np.array(mol_input_fps, dtype=np.float64)
                substructure_matrix = substructure_matrix[:, np.sum(substructure_matrix, axis=0) != 0]
                mol_num, substructure_num = substructure_matrix.shape
                substructure_links = []
                for mol in range(mol_num):
                    for i in range(substructure_num):
                        if substructure_matrix[mol, i] == 1:
                            substructure_links.append([mol, mol_num + i])
                substructure_links = [item + [substr_weight] for item in substructure_links]
                mol_interact_relations = [item + [1 / substr_weight] for item in mol_interact_relations]
                links = mol_interact_relations + substructure_links
                mat_nodes = list(itertools.chain.from_iterable(links))
                mat_nodes = set(mat_nodes)
                mat_nodes = {np.int32(node) for node in mat_nodes}
                mat_size = np.int32(max(mat_nodes) + 1)
                network = np.zeros((mat_size, mat_size))
                for item in links:
                    network[np.int32(item[0]), np.int32(item[1])] = item[2]
                network = network + network.T
                degree = np.tile(np.sum(network, 0), (np.size(network, 0), 1))
                degree = degree.T
                trans_mat = np.divide(network, degree)
                unit_mat = np.eye(np.size(trans_mat, 0))
                score = unit_mat
                step_i = 0
                while step_i < 3:
                    score = np.dot(trans_mat.T, score)
                    step_i += 1
                score = score + score.T
                score[np.isnan(score)] = 0
                score[np.isinf(score)] = 0
                list1 = molecules_index
                list2 = input_all_smi
                molecules_dict = dict(zip(list1, list2))
                def calculate_scores(score_matrix, nodes):
                    node_count = len(nodes)
                    scores = []
                    for i in range(node_count):
                        for j in range(i + 1, node_count):
                            node1, node2 = nodes[i], nodes[j]
                            scores.append(
                                (
                                    node1,
                                    node2,
                                    score_matrix[node1, node2]
                                )
                            )
                    scores.sort(key=lambda x: x[2], reverse=True)
                    return scores
                molecules_scores = calculate_scores(score, molecules_index)
                df = pd.DataFrame(molecules_scores, columns=['molecule1', 'molecule2', 'score'])
                df['rank'] = range(1, len(df) + 1)
                all_combinations = []
                for index, row in df.iterrows():
                    combination_1 = tuple([row['molecule1'], row['molecule2']])
                    combination_2 = tuple([row['molecule2'], row['molecule1']])
                    all_combinations.append(combination_1)
                    all_combinations.append(combination_2)
                list_file = mol_interact_relations_copy
                list_file = [tuple(inner_list) for inner_list in list_file]
                known_combinations = list_file
                df['status'] = 'unknown'
                for index, row in df.iterrows():
                    if (row['molecule1'], row['molecule2']) in known_combinations or (
                            row['molecule2'], row['molecule1']) in known_combinations:
                        df.at[index, 'status'] = 'known'
                df['molecule1'] = df['molecule1'].map(molecules_dict)
                df['molecule2'] = df['molecule2'].map(molecules_dict)
                df.to_csv('./results/NTDDCNBR.csv', index=False)
                file_name = 'NTDDCNBR'
                def mol_to_excel(file):
                    header = ['molecule1', 'molecule2', 'score', 'rank', 'status', 'SMILES1', 'SMILES2']
                    excel_file_path = os.path.join("results", f'{file}.xlsx')
                    workbook = xlsxwriter.Workbook(excel_file_path)
                    ItemStyle = workbook.add_format(item_style)
                    HeaderStyle = workbook.add_format(header_style)
                    worksheet = workbook.add_worksheet()
                    worksheet.set_column('A:B', 38)
                    worksheet.set_column('C:C', 10)
                    for ix_, i in enumerate(header):
                        worksheet.write(0, ix_, i, HeaderStyle)
                    df = pd.read_csv(os.path.join("results", f"{file}.csv"), encoding='gbk').head(DDC_num)
                    for i in range(df.shape[0]):
                        molecule1 = df.iloc[i, 0]
                        molecule2 = df.iloc[i, 1]
                        score = df.iloc[i, 2]
                        rank = df.iloc[i, 3]
                        status = df.iloc[i, 4]
                        SMILES1 = df.iloc[i, 0]
                        SMILES2 = df.iloc[i, 1]
                        img_data_molecule1 = BytesIO()
                        img_data_molecule2 = BytesIO()
                        c_molecule1 = Chem.MolFromSmiles(molecule1)
                        c_molecule2 = Chem.MolFromSmiles(molecule2)
                        img_molecule1 = Draw.MolToImage(c_molecule1)
                        img_molecule2 = Draw.MolToImage(c_molecule2)
                        img_molecule1.save(img_data_molecule1, format='PNG')
                        img_molecule2.save(img_data_molecule2, format='PNG')
                        worksheet.set_row(i + 1, 185)
                        worksheet.insert_image(i + 1, 0, 'f',
                                               {'x_scale': 0.9, 'y_scale': 0.8, 'image_data': img_data_molecule1,
                                                'positioning': 1})
                        worksheet.insert_image(i + 1, 1, 'f',
                                               {'x_scale': 0.9, 'y_scale': 0.8, 'image_data': img_data_molecule2,
                                                'positioning': 1})
                        worksheet.write(i + 1, 2, score, ItemStyle)
                        worksheet.write(i + 1, 3, rank, ItemStyle)
                        worksheet.write(i + 1, 4, status, ItemStyle)
                        worksheet.write(i + 1, 5, SMILES1, ItemStyle)
                        worksheet.write(i + 1, 6, SMILES2, ItemStyle)
                    workbook.close()
                mol_to_excel(file_name)
            filename = f"{Model}.xlsx"
            source_path = os.path.abspath(os.path.join("results", filename))
            desktop_path = os.path.expanduser("~/Desktop")
            destination_path = os.path.join(desktop_path, filename)
            try:
                shutil.copy2(source_path, destination_path)
                messagebox.showinfo("Success",
                                    f"Recommendation process complete. Excel file saved on desktop: {destination_path}")
            except Exception as e:
                messagebox.showinfo("Error", f"Error: {e}")
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {str(e)}")
def main():
    root = tk.Tk()
    app = WelcomeWindow(root)
    root.mainloop()
if __name__ == "__main__":
    main()
