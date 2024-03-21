#!/usr/bin/env python
import argparse
import os
import shutil
import subprocess
import csv
import numpy as np
import pandas as pd
import glob
import re

###################################################################################################
#print('FlexiSimCont: \nUsage: FlexiSimCont.py '
#      '[-contact ]    <Return 1 if the distance between atoms is less than 4.5 Å.Return the average value within the analyzed time interval (0 < value < 1) when analyzing the trajectories. Residue information in the format: pro_chain1,pro_number1,pro_atomtype1,...>'
#      '[-hbond ]      <Return 1 if the hydrogen bond present. Return the average value within the analyzed time interval (0 < value < 1) when analyzing the trajectories. Residue information in the format: pro_chain1,pro_number1,pro_atomtype1,lig_atomtype1,...>'
#      '[-hbond_res ]  <Return 1 if the hydrogen bond present. Return the average value within the analyzed time interval (0 < value < 1) when analyzing the trajectories. Residue information in the format: pro_chain1,pro_number1,pro_atomtype1,lig_atomtype1,...>'
#      '[-dockpose ]   <Enable DockPose mode. Analyze the non-covalent interactions in docking pose.>'
#      '[-t]           <Time for analysis, e.g. 1-200, defaulting to last 5 ns of MD simulation trajectories> '
#      '[-workdir]     <The work directory must be set to the same directory as FlexiSimScr.py, defaulting to the current directory if not specified> '
#      '[-num ]        <Number of ligand, e.g. 1, defaulting to the total number of small molecules>'
#print ('Command: \npython FlexiSimCont.py -num 1 -t 100 -gpuid 0 &')
###################################################################################################

#### 1. Parse input line & Preparation ####

# Parse input line
parser = argparse.ArgumentParser()
parser.add_argument("-dockpose", action="store_true", help="Enable DockPose mode")
parser.add_argument("-t", type=str, required=False, help="Time_for_analysis", default="10")
parser.add_argument("-pbsa", action="store_true", help="Enable PB model")
parser.add_argument("-gbsa", action="store_true", help="Enable GB model ")
parser.add_argument("-receptor", type=str, required=False, help="Residue information in the format: pro_chain1,pro_number1,pro_atomtype1,...")
parser.add_argument("-ligand", type=str, required=False, help="Residue information in the format: pro_chain1,pro_number1,pro_atomtype1,...")
parser.add_argument("-workdir", type=str, required=False, help="Work directory")
parser.add_argument("-num", type=str, required=False,help="Number of ligands, default set to the total number of small molecules.")
args = parser.parse_args()

num = int(args.num) if args.num is not None else 1 # Number setting part1

workdir = args.workdir if args.workdir else None # Work directory setting

if workdir is None:
    print('Work directory not specified, setting to current directory.')
    workdir = os.path.abspath(os.getcwd())

# Preparation 
os.chdir(workdir)
complex_md_dir = os.path.join(workdir, 'complex_md')

if not os.path.exists(complex_md_dir):
    print("Error: 'complex_md' directory not found. Exiting the program.")
    exit()
else:
    os.chdir(complex_md_dir)

# User needs to set this manually:
subprocess.run(['/bin/bash', '-c', 'source /home/soft/amber22/amber.sh'])

print ('Preparation has been completed.') 

#### 2. Performing calculations individually for each ligand ####

###################################################################################################

# read_csv
def read_csv(file_path):
    data = []
    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            data.append(row)
    return np.array(data)

# write_csv
def write_csv(file_path, data):
    with open(file_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(data)

# read_mmpbsa
def read_mmpbsa():
    try:
        with open('mmpbsa_raw1.csv', 'r') as file:
            reader = csv.reader(file)
            mmpbsa_data = list(reader)
        return np.array(mmpbsa_data)
    except FileNotFoundError:
        return None

# write_mmpbsa
def write_mmpbsa(data):
    with open('mmpbsa_raw1.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(data)

#calculate_simulation_time
def calculate_simulation_time(file_path):
    try:
        with open(file_path, 'r') as file:
            content = file.read()
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return None

    nstlim_match = re.search(r'nstlim\s*=\s*(\d+)', content)
    dt_match = re.search(r'dt\s*=\s*([0-9.]+)', content)

    if nstlim_match and dt_match:
        nstlim_0 = int(nstlim_match.group(1))
        dt_0 = float(dt_match.group(1))
        time = nstlim_0 * dt_0/1000
        return time
    else:
        print("Failed to extract nstlim and dt from the file.")
        return None

# count_mol2_files()
def count_ligand_files():
    ligand_files = [file for file in os.listdir() if file.startswith('lig') and file.endswith('.mol2')]
    return len(ligand_files)

###################################################################################################
# setting number to 'total number of given ligands'
if args.num is None:
    lig_folders = [folder for folder in os.listdir() if os.path.isdir(folder) and folder.startswith('lig')]

    number_ligand = len(lig_folders) 
    num = number_ligand

# loop part
for stage in range(1, int(num) + 1):
  
    stage = str(stage)
    lig_dir = os.path.join(os.getcwd(), f'lig{stage}')

    if not os.path.exists(lig_dir):
        print(f"Error: 'lig{stage}' directory does not exist. Exiting the program.")
        exit()
    else:
        os.chdir(lig_dir)

    #### 2.1 Preparation for receptor&ligand&complex
        
    #Pdb files for receptor&ligand&complex
    
    mmpbsa_folder = 'mmpbsa'
    os.makedirs(mmpbsa_folder, exist_ok=True)

    shutil.copy('LIG.prmtop', os.path.join(mmpbsa_folder, 'ligand.prmtop'))
    shutil.copy('LIG.inpcrd', os.path.join(mmpbsa_folder, 'ligand.inpcrd'))
    shutil.copy('LIG_raw.pdb', os.path.join(mmpbsa_folder, 'ligand.pdb'))
    shutil.copy('LIG.frcmod', os.path.join(mmpbsa_folder, 'LIG.frcmod'))
    shutil.copy('LIG.lib', os.path.join(mmpbsa_folder, 'LIG.lib'))

    shutil.copy('comp_raw2.pdb', os.path.join(mmpbsa_folder, 'complex.pdb'))
    
    os.chdir(mmpbsa_folder)

    with open('complex.pdb', 'r') as input_file, open('receptor.pdb', 'w') as output_file:
        for line in input_file:
            if 'LIG ' not in line:
                output_file.write(line)
  
    # Prmtop&inpcrd files for receptor&ligand&complex
    leap = f'''
    source leaprc.protein.ff19SB
    source leaprc.gaff2
    source leaprc.water.opc
    loadamberparams frcmod.ionslm_1264_opc
    loadoff LIG.lib
    loadamberparams LIG.frcmod

    com = loadpdb complex.pdb
    saveamberparm com com.prmtop com.inpcrd
    recep = loadpdb receptor.pdb
    saveamberparm recep rec.prmtop rec.inpcrd
    lig = loadpdb ligand.pdb
    saveamberparm lig lig.prmtop lig.inpcrd
    quit'''
        
    with open('leap.in', 'w') as f_out:
        f_out.write(leap)

    with open('leap.log', 'w') as log_file:
        process = subprocess.Popen(['/home/soft/amber22/bin/tleap', '-s', '-f', 'leap.in'], stdout=log_file, stderr=subprocess.PIPE)

    process.communicate()

    # nc file for complex
    mdrun_in_path = '../mdrun.in'  
    simulation_time = calculate_simulation_time(mdrun_in_path)

    time_for_analysis = int(args.t) if args.t else None # Time setting
    if time_for_analysis is None or time_for_analysis <= 0:
        print('Analysis time is not provided or invalid, setting to default: 5 ns.')
        time_for_analysis = 5

    simulation_frames = simulation_time * 100  
    analysis_frames = time_for_analysis * 100 # time_for_analysis = int(args.t)
    start_frame = simulation_frames - analysis_frames + 1
    interval = 1
    end_frame = simulation_frames
    end_frame_pbsa = end_frame/interval

    start_frame = int(start_frame)
    end_frame = int(end_frame)
    end_frame_pbsa = int(end_frame_pbsa)
    
    # parameters setting
    cation_name = "Na+" 
    anion_name = "Cl-"
    #The types of ions that can be used include:
    #Li+,Na+,K+,Rb+,Cs+,Tl+,Cu+,Ag+,NH4,HE+,HZ+,H3O+
    #F-,Cl-,Br-,I-,
    #Be2+,Cu2+,Ni2+,Pt2+,Zn2+,Co2+,Pd2+,Ag2+,Cr2+,Fe2+,Mg2+,V2+,
    #Mn2+,Hg2+,Cd2+,Yb2+,Ca2+,Sn2+,Pb2+,Eu2+,Sr2+,Sm2+,Ba2+,Ra2+
    #Al3+,Fe3+,Cr3+,In3+,Tl3+,Y3+,La3+,Ce3+,Pr3+,Nd3+,Sm3+,Eu3+,
    #Gd3+,Tb3+,Dy3+,Er3+,Tm3+,Lu3+
    #Hf4+,Zr4+,Ce4+,U4+,Pu4+,Th4+

    # DockPose mode
    if args.dockpose: 
        cpp = f'''parm ../comp_sol.prmtop
        trajin ../comp_sol.inpcrd
        center !(:WAT,{cation_name},{anion_name})
        image center
        strip :WAT,{cation_name},{anion_name}
        trajout com.nc
        go''' 
     
    else:
        cpp = f'''
        parm ../comp_sol.prmtop
        trajin ../mdrun.nc {start_frame} {end_frame} {interval}
        center !(:WAT,{cation_name},{anion_name})
        image center
        strip :WAT,{cation_name},{anion_name}
        trajout com.nc
        go'''
        
    with open('cpp.in', 'w') as f_out:
        f_out.write(cpp)

    with open('cpp.log', 'w') as log_file:
        process = subprocess.Popen(['cpptraj', '-i', 'cpp.in'], stdout=log_file, stderr=subprocess.PIPE)

    process.communicate()

    print (f'{stage}.1 The input files for MMPBSA calculation have been generated.') 

    #### 2.2 calculation for MMPBSA
    if args.dockpose: 
        end_frame_pbsa=1

    pbsa = f'''
Input file for running PB and GB
&general
  startframe=1, endframe={end_frame_pbsa}, verbose=1,
/
&gb
  igb=5, saltcon=0.150,
/
&pb
  istrng=0.100,
/\n'''
 
    gbsa = f'''
Input file for running PB and GB
&general
  startframe=1, endframe={end_frame_pbsa}, verbose=1,
/
&gb
  igb=5, saltcon=0.150,
/\n'''
        
    with open('pbsa.in', 'w') as f_out:
        f_out.write(pbsa)

    with open('gbsa.in', 'w') as f_out:
        f_out.write(gbsa)

    if args.gbsa: 

        command1 = "$AMBERHOME/bin/MMPBSA.py -O -i gbsa.in -o FINAL_RESULTS_MMGBSA.dat -cp com.prmtop -rp rec.prmtop -lp lig.prmtop -y com.nc"
        process = subprocess.Popen(command1, shell=True, stderr=subprocess.PIPE)
        process.communicate()

    if args.pbsa:
        command2 = "$AMBERHOME/bin/MMPBSA.py -O -i pbsa.in -o FINAL_RESULTS_MMPBSA.dat -do FINAL_DECOMP_MMPBSA_pbsa.dat -cp com.prmtop -rp rec.prmtop -lp lig.prmtop -y com.nc"
        process = subprocess.Popen(command2, shell=True, stderr=subprocess.PIPE)
        process.communicate()

    # ID
    id_data = stage
    
    mmpbsa_data = read_mmpbsa()

    if mmpbsa_data is not None:
        new_column = [f"{id_data}" for _ in range(len(mmpbsa_data))]
        mmpbsa_data = np.column_stack((mmpbsa_data, new_column))
    else:
        mmpbsa_data = [[f"{id_data}"]]

    write_mmpbsa(mmpbsa_data)   

    # Title
    with open('../LIG_raw.mol2', 'r') as file:
        lines = file.readlines()

    if len(lines) >= 2:
        title_data = lines[1].strip()
    else:
        print("Error: File 'LIG_raw.mol2' doesn't have at least two lines.")
    
    mmpbsa_data = read_mmpbsa()

    if mmpbsa_data is not None:
        new_column = [f"{title_data}" for _ in range(len(mmpbsa_data))]
        mmpbsa_data = np.column_stack((mmpbsa_data, new_column))
    else:
        mmpbsa_data = [[f"{title_data}"]]
    
    write_mmpbsa(mmpbsa_data)
    
    # Get delta total&Std. dev.
    
    if args.gbsa:
        file_path = 'FINAL_RESULTS_MMGBSA.dat'
        with open(file_path, 'r') as file:

            for line in file:
                if line.startswith('DELTA TOTAL'):
                    data = line.split()
                    delta_total = float(data[2])
                    std_dev = float(data[3])
                    break
        
        delta_total = round(delta_total,2)
        std_dev = round(std_dev,2)

        mmpbsa_data = read_mmpbsa()

        if mmpbsa_data is not None:
            new_column_contact = [f"{delta_total}" for _ in range(len(mmpbsa_data))]
            mmpbsa_data = np.column_stack((mmpbsa_data, new_column_contact))
        else:
            mmpbsa_data = [[f"{delta_total}"]]
            
        write_mmpbsa(mmpbsa_data)

        mmpbsa_data = read_mmpbsa()

        if mmpbsa_data is not None:
            new_column_contact = [f"{std_dev}" for _ in range(len(mmpbsa_data))]
            mmpbsa_data = np.column_stack((mmpbsa_data, new_column_contact))
        else:
            mmpbsa_data = [[f"{std_dev}"]]
            
        write_mmpbsa(mmpbsa_data)

    if args.pbsa:
        file_path = 'FINAL_RESULTS_MMPBSA.dat'
        with open(file_path, 'r') as file:
        
            for line in file:
                if line.startswith('DELTA TOTAL'):
                    data = line.split()
                    delta_total = float(data[2])
                    std_dev = float(data[3])
                    break

        delta_total = round(delta_total,2)
        std_dev = round(std_dev,2)

        mmpbsa_data = read_mmpbsa()

        if mmpbsa_data is not None:
            new_column_contact = [f"{delta_total}" for _ in range(len(mmpbsa_data))]
            mmpbsa_data = np.column_stack((mmpbsa_data, new_column_contact))
        else:
            mmpbsa_data = [[f"{delta_total}"]]
            
        write_mmpbsa(mmpbsa_data)

        mmpbsa_data = read_mmpbsa()

        if mmpbsa_data is not None:
            new_column_contact = [f"{std_dev}" for _ in range(len(mmpbsa_data))]
            mmpbsa_data = np.column_stack((mmpbsa_data, new_column_contact))
        else:
            mmpbsa_data = [[f"{std_dev}"]]
            
        write_mmpbsa(mmpbsa_data)

    print (f'{stage}.2 The mmpbsa_raw1.csv files has been generated.') 

    #### 2.3 mmpbsa_raw.csv 
    source_csv_file = './mmpbsa_raw1.csv'
    destination_csv_file = '../../mmpbsa_raw.csv'

    with open(source_csv_file, 'r') as source_file:
        last_line = source_file.readlines()[-1]

    with open(destination_csv_file, 'a') as dest_file: # Append the last line to the destination CSV file
        dest_file.write(last_line)

    print(f"{stage}.3 Data has been appended to mmpbsa_raw.csv.")

    #### 2.4 Go back to the parent directory. ####
    os.chdir(complex_md_dir)

#### 3. get mmpbsa: mmpbsa.csv ####

title_raw_path = 'title_raw.csv'
if os.path.exists(title_raw_path):
    title_raw_data = read_csv(title_raw_path)
else:
    title_raw_data = []

# title of row
id_row = 'ID'
title_row = 'Title'
delta_total_row = 'Delta_Total'
std_dev_row = 'Std_Dev'

row_to_append = [id_row, title_row, delta_total_row,std_dev_row] 

title_raw_data.append(row_to_append)

write_csv(title_raw_path, title_raw_data) # Write to title_raw.csv

# mmpbsa_raw.csv + title_raw.csv --> mmpbsa.csv
title_raw_path = 'title_raw.csv'
mmpbsa_raw_path = 'mmpbsa_raw.csv'
mmpbsa_path = 'mmpbsa.csv'

# Read the content of both files
with open(title_raw_path, 'r') as title_raw_file:
    title_raw_content = title_raw_file.read()

with open(mmpbsa_raw_path, 'r') as mmpbsa_raw_file:
    mmpbsa_raw_content = mmpbsa_raw_file.read()

# Concatenate the content and write to the output file
with open(mmpbsa_path, 'w') as mmpbsa_file:
    mmpbsa_file.write(title_raw_content + mmpbsa_raw_content)

# Delete raw files.
pattern = '*raw.csv'
[os.remove(file) and print(f"Deleted: {file}") for file in glob.glob(pattern)]

print("All data has been saved in mmpbsa.csv.")

###################################################################################################
# List of Scripts for Flexible Simulation and Stable Analysis Post Screening
#
#1. FlexiSimScr.py:     Flexible Simulation for Screening
#2. FlexiSimComp.py:    Flexible Simulation for complex
#3. FlexiSimPro.py:     Flexible Simulation for protein
#4. FlexiSimCont.py:    Calculate contact/stablity after Docking/MD simulation
#5. FlexiSimPBSA.py:    Calculate binding free energy after Docking/MD simulation

###################################################################################################

#1## Match the "time" parameter with the FlexSimScr.py script.
#
#2## Specify the analysis duration using the "t" parameter for the last 't' nanoseconds of the trajectory.
#
#3## Keep a consistent format for "contact" and "hbond" parameters. 
#### Each group includes chain, residue number, atom name of protein atom 1, and atom name of small molecule atom 2. Add another set for additional criteria.
#
#4## Enabling the "dockpose" parameter activates DockPose mode, analyzing the interactions between proteins and small molecules in the docking poses.

#example1:
# nohup python FlexiSimPBSA.py -dockpose -gbsa
#example2:
# nohup python FlexiSimPBSA.py -dockpose -gbsa
#current:
# nohup python FlexiSimPBSA.py -t 5 -gbsa&