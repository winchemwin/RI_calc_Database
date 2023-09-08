import customtkinter as ctk
import tkinter as tk
from tkinter import ttk
from tkinter import StringVar, messagebox
from turtle import onclick
from tkinter import filedialog
import numpy as np
import pandas as pd
from numpy import matlib
from tkinter.ttk import Entry
import psi4
from  rdkit import rdBase, Chem
from rdkit.Chem import AllChem
import re

# Smiles file select by filedialog 
def data_import(self):
    type=[('smi file','*.smi'), ('text file', '*.txt')]
    file=filedialog.askopenfilename(filetypes=type)
    filename_sv.set(file)
    global fname_smi
    fname_smi=filename_sv.get()

    return fname_smi

def start_calc_button(self):
    Labelex_1=ctk.CTkLabel(master=root, text='Calculation was started',font=("Times",12))
    Labelex_1.place(x=50, y=370)
    
    root.after(100, calc_run)
    
def Geom_Opt():
    # Geometry optimization
    # getting meth(calculation method), func(function), base(basis set) 
    meth=method_sv.get()
    func=function_sv.get()
    base=baseset_sv.get()
    
    psi4.set_options({'geom_maxiter': 50,
    'optking__g_convergence': 'gau_loose',})
    
    # geometry optimize calculation (method in HF and MP2, function in DFT)
    if meth=='HF':
        psi4.optimize(meth+'/'+base, molecule=molgeom)
        
    elif meth=='DFT':
        
        psi4.optimize(func+'/'+base, molecule=molgeom)
        
    elif meth=='MP2':
        psi4.optimize(meth+'/'+base, molecule=molgeom)
    
    optimized_geom = molgeom.save_string_xyz_file()
                  
def Vib_Calc():
    # Frequency calculation
    meth=method_sv.get()
    func=function_sv.get()
    base=baseset_sv.get()
    
    # Frequency calculation (method in HF and MP2, function in DFT)
    if meth=='HF':        
        energy, wfn = psi4.frequency(meth+'/'+base, molecule=molgeom, return_wfn=True)
        
    elif meth=='DFT':
        energy, wfn = psi4.frequency(func+'/'+base, molecule=molgeom, return_wfn=True)
        
    elif meth=='MP2':
        energy, wfn = psi4.frequency(meth+'/'+base, molecule=molgeom, return_wfn=True)
    
    # making frequency list

    freqlist = np.array(wfn.frequencies())
    freqposi=freqlist[(freqlist>1)&(freqlist<5000)]
    global freqnega
    freqnega=freqlist[(freqlist<0)]

    return freqposi, freqnega

def Pol_calc():
    meth_pol=method_sv_pol.get()
    base=baseset_sv.get()
          
    polar = psi4.properties(meth_pol+'/'+base,properties=['polarizability'], molecule=molgeom)
        

def calc_run():
    global fname_smi
    chr=charge.get()
    mul=multi.get()
    
    with open (fname_smi) as f:
        readsmiles=f.read()
    
    smiles=readsmiles.split()
    
    freq_results=[]

    smiIndex=[]
    
    polar_list=[]
    volume_list=[]
    refIndex_list=[]
      
    #  Conformer search by rdkit then make xyz file
    for i, smi in enumerate(smiles): 

        mol = Chem.MolFromSmiles(smi)
        mol_H = Chem.AddHs(mol)
            
        confs = AllChem.EmbedMultipleConfs(mol_H, 10, pruneRmsThresh=1)
        prop = AllChem.MMFFGetMoleculeProperties(mol_H)

        energy=[]
        
        for conf in confs:
            mmff = AllChem.MMFFGetMoleculeForceField(mol_H, prop,confId=conf)
            mmff.Minimize()
            energy.append((mmff.CalcEnergy(), conf))
    
            conflist = np.array(energy)
            sortconf=conflist[np.argsort(conflist[:,0]) ]

            stconfid=sortconf[0,1]
            stconfid2=int(stconfid)
            stconfgeom=mol_H.GetConformer(stconfid2)

        xyz = f'{chr} {mul}'
        for atom, (x,y,z) in zip(mol_H.GetAtoms(), stconfgeom.GetPositions()):
            xyz += '\n'
            xyz += '{}\t{}\t{}\t{}'.format(atom.GetSymbol(), x, y, z)
                  
        #  Setting input file
        global molgeom
        molgeom = psi4.geometry(xyz)
               
    # Set calculation method
        thre=threads.get()
        mem=memory.get()
        psi4.set_num_threads(nthread=thre)
        psi4.set_memory(mem*1000000)

        psi4.core.clean()
        # Set output files
        psi4.set_output_file(str(i)+smi+'-molgeom.txt')

        # global root, Labelex_2
        Labelex_2=ctk.CTkLabel(master=root, text='No'+str(i)+'  '+smi+' Calc was started',font=("Times",10))
        Labelex_2.place(x=50, y=390)

        Labelex_2.update()
        
        try:
            Geom_Opt()
            root.update()
        except  psi4.driver.p4util.exceptions.OptimizationConvergenceError as err:
            print (err)
       
            continue
        
        optimized_geom = molgeom.save_string_xyz_file()
        with open (str(i)+smi+'-optgeom.xyz','w') as f:
            f.write(optimized_geom)
        
        Vib_Calc()
        root.update()
        if freqnega:
            freq_data='True'         
        else:
            freq_data='False'
            
        freq_results.append(freq_data)  
        
        Pol_calc()
        root.update()
        
        path = '/home/winmori/'+str(i)+smi+'-molgeom.txt'
        with open (path) as f:
            alltext= f.readlines()
        endgyou = len(alltext)
        endtxt = alltext[endgyou-1].strip()
        print (endtxt)

        polar_base=re.findall(r"[\d\.\d]+",endtxt)
        
        polar_result=float(polar_base[1])
        
        # Molvolume 
        try:
            AllChem.EmbedMolecule(mol_H)
            rdvdwvol=AllChem.ComputeMolVolume(mol_H)
        except ValueError:
            pass
        
	# Refindex calculation
        vol=(1/(rdvdwvol*(10**-30)))
        fai=(4*3.14159/3)*vol*(polar_result*0.148185*10**-30)
        refindex=((1+2*fai)/(1-fai))**0.5


        refIndex_list.append(refindex)
        polar_list.append(polar_result)
        volume_list.append(rdvdwvol)
        
        smiIndex.append(smi)
        
        Data_RI=pd.DataFrame(list(zip(smiIndex, polar_list, volume_list, refIndex_list)) , columns=['Smiles','polarizability/a.u', 'VdwVolume/Angstrom3/mol', 'Refractive Index']) 
       
        Data_RI['FreqNega']=freq_results
    
        Data_RI.to_csv('Calcd_RI_List.csv')     
    
    Labelex_4=ctk.CTkLabel(master=root, text='Calculation was finished', font=("Times",12,"bold"))
    Labelex_4.place(x=50, y=410)
                
# finish program

def scry_finish():
    exit()        
           
#　Custon Tkinter main 

ctk.set_appearance_mode('dark')
# global root
root = ctk.CTk()
root.title("Reflactive Index List Calculator") 
root.geometry('700x450')

# Select main molecule file
Label1=ctk.CTkLabel(master=root, fg_color='blue',corner_radius=10, text='Smiles File Select',font=("Arial",16,"bold","italic"))
Label1.place(x=20, y=20)

filename_sv = tk.StringVar()
filenameEntry = ctk.CTkEntry(master=root, width=500,  textvariable=filename_sv)
filenameEntry.place(x=20, y= 50)

Button1 = ctk.CTkButton(master=root, text='Select',width=30,font=("Times",14,"bold"))
Button1.bind("<Button-1>", data_import) 
Button1.place(x=600, y=50)



# Calculation Run
Button2_3=ctk.CTkButton(master=root, text='Quantum Calculation Start',width=50, font=("Times",14,"bold"))
Button2_3.bind("<Button-1>", start_calc_button) 
Button2_3.place(x=70, y=330)


# Select calculation methods
Label3=ctk.CTkLabel(master=root, fg_color='blue',corner_radius=10, text='Quantum Calculation Method', font=("Arial",16,"bold","italic"))
Label3.place(x=20, y=100)

Label3_1=ctk.CTkLabel(master=root, text='Method',font=("Times",12))
Label3_1.place(x=70, y=140)
method_sv = tk.StringVar()
methods=('HF','DFT', 'MP2')
comboBox3_1=ctk.CTkComboBox(master=root, height=5, width=80, border_color='blue',state='readonly', values=methods, variable=method_sv)
comboBox3_1.place(x=70, y=160)

Label3_2=ctk.CTkLabel(master=root, text='Function',font=("Times",12))
Label3_2.place(x=170, y=140)
function_sv = tk.StringVar()
functions=('','b3lyp','cam-b3lyp', 'edf2','m06', 'pbe','wb97x-d')
comboBox3_2=ctk.CTkComboBox(master=root, height=5, width=80, border_color='blue',state='readonly', values=functions, variable=function_sv)
comboBox3_2.place(x=170, y=160)

Label3_3=ctk.CTkLabel(master=root, text='Basis set',font=("Times",12))
Label3_3.place(x=270, y=140)
baseset_sv = tk.StringVar()
base_sets=('3-21g','6-31g', '6-31g(d)','6-311g', 'aug-cc-pvtz')
comboBox3_3=ctk.CTkComboBox(master=root, height=5, width=80,border_color='blue', state='readonly', values=base_sets, variable=baseset_sv)
comboBox3_3.place(x=270, y=160)

# Select options
Label4=ctk.CTkLabel(master=root, text='Options', font=("Times",14,"bold"))
Label4.place(x=70, y=180)

Label4_1=ctk.CTkLabel(master=root, text='Tread',font=("Times",12))
Label4_1.place(x=80, y=200)
threads=tk.IntVar(value=2)
textBox1_1=ctk.CTkEntry(master=root, width=50, textvariable=threads)
textBox1_1.place(x=80, y=220)

Label4_2=ctk.CTkLabel(master=root, text='Memory/MB',font=("Times",12))
Label4_2.place(x=180, y=200)
memory=tk.IntVar(value=500)
textBox1_2=ctk.CTkEntry(master=root, width=50, textvariable=memory)
textBox1_2.place(x=180, y=220)

Label4_3=ctk.CTkLabel(master=root, text='Charge',font=("Times",12))
Label4_3.place(x=80, y=250)
charge=tk.IntVar(value=0)
textBox2_1=ctk.CTkEntry(master=root, width=50, textvariable=charge)
textBox2_1.place(x=80, y=270)

Label4_4=ctk.CTkLabel(master=root, text='Multiplicity',font=("Times",12))
Label4_4.place(x=180, y=250)
multi=tk.IntVar(value=1)
textBox2_2=ctk.CTkEntry(master=root, width=50, textvariable=multi)
textBox2_2.place(x=180, y=270)

# select Polar calc method

Label5=ctk.CTkLabel(master=root, fg_color='blue',corner_radius=10, text='Polarizability Calc Set', font=("Arial",12,"bold","italic"))
Label5.place(x=400, y=100)

Label5_1=ctk.CTkLabel(master=root, text='Calculation method',font=("Times",12))
Label5_1.place(x=420, y=130)
method_sv_pol = tk.StringVar()
methods=('cc2','ccsd')
comboBox5_1=ctk.CTkComboBox(master=root, height=5, width=80, border_color='blue',state='readonly', values=methods, variable=method_sv_pol)
comboBox5_1.place(x=420, y=160)


# Finish program
Label6=ctk.CTkLabel(master=root, text='Finish the program')
Label6.place(x=550, y=370)
Button3 = ctk.CTkButton(master=root, text='Quit',width=30, command=scry_finish,font=("Times",14,"bold"))
Button3.place(x=550, y=390)


root.mainloop()
