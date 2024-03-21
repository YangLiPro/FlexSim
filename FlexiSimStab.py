#!/usr/bin/env python
import argparse
import os
import shutil
import subprocess
import csv
import numpy as np
import pandas as pd
from rdkit import Chem
import glob
import re

###################################################################################################
#print('FlexiSimStab: \nUsage: FlexiSimStab.py '
#      '[-contact ]    <Return 1 if the distance between atoms is less than 4.5 Å.Return the average value within the analyzed time interval (0 < value < 1) when analyzing the trajectories. Residue information in the format: pro_chain1,pro_number1,pro_atomtype1,...>'
#      '[-hbond ]      <Return 1 if the hydrogen bond present. Return the average value within the analyzed time interval (0 < value < 1) when analyzing the trajectories. Residue information in the format: pro_chain1,pro_number1,pro_atomtype1,lig_atomtype1,...>'
#      '[-hbond_res ]  <Return 1 if the hydrogen bond present. Return the average value within the analyzed time interval (0 < value < 1) when analyzing the trajectories. Residue information in the format: pro_chain1,pro_number1,pro_atomtype1,lig_atomtype1,...>'
#      '[-dockpose ]   <Enable DockPose mode. Analyze the non-covalent interactions in docking pose.>'
#      '[-t]           <Time for analysis, e.g. 1-200, defaulting to last 5 ns of MD simulation trajectories> '
#      '[-workdir]     <The work directory must be set to the same directory as FlexiSimScr.py, defaulting to the current directory if not specified> '
#      '[-num ]        <Number of ligand, e.g. 1, defaulting to the total number of small molecules>'
#print ('Command: \npython FlexiSimStab.py -num 1 -t 100 -gpuid 0 &')
###################################################################################################

#### 1. Parse input line & Preparation ####

# Parse input line
parser = argparse.ArgumentParser()
parser.add_argument("-contact", type=str, required=False, help="Residue information in the format: pro_chain1,pro_number1,pro_atomtype1,...")
parser.add_argument("-hbond", type=str, required=False, help="Residue information in the format: pro_chain1,pro_number1,pro_atomtype1,...")
parser.add_argument("-hbond_res", type=str, required=False, help="Residue information in the format: pro_chain1,pro_number1,...")
parser.add_argument("-dockpose", action="store_true", help="Enable DockPose mode")
parser.add_argument("-t", type=str, required=False, help="Time_for_analysis", default="10")
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

# find_N_atom_coordinates
def find_N_atom_coordinates(protein_pdb, chain, residue_number):

    n_atom_coordinates = []

    with open(protein_pdb, 'r') as file:

        for line in file:

            if line.startswith('ATOM') and line[21] == chain and int(line[22:26]) == residue_number and line[12:16].strip() == 'N':
                x_coord = float(line[30:38])
                y_coord = float(line[38:46])
                z_coord = float(line[46:54])
                n_atom_coordinates.append((x_coord, y_coord, z_coord))
    return n_atom_coordinates

# find_atom_coordinates
def find_atom_coordinates(protein_pdb, chain, residue_number, atom_name):

    atom_coordinates = []

    with open(protein_pdb, 'r') as file:
        for line in file:
            if line.startswith('ATOM') and line[21] == chain and int(line[22:26]) == residue_number and line[12:16].strip() == atom_name:
                x_coord = float(line[30:38])
                y_coord = float(line[38:46])
                z_coord = float(line[46:54])
                atom_coordinates.append((x_coord, y_coord, z_coord))

    return atom_coordinates

# find_lig_atom_coordinates
def find_lig_atom_coordinates(ligand_pdb, atom_name):
    lig_atom_coordinates = []

    with open(ligand_pdb, 'r') as file:
        for line in file:
            if (line.startswith('ATOM') or line.startswith('HETATM')) and line[12:16].strip() == atom_name:
                x_coord = float(line[30:38])
                y_coord = float(line[38:46])
                z_coord = float(line[46:54])
                lig_atom_coordinates.append((x_coord, y_coord, z_coord))

    return lig_atom_coordinates

# find_residue_number
def find_residue_number(comp_pdb, n_atom_coords):

    residue_numbers = []

    with open(comp_pdb, 'r') as file:

        for line in file:

            if line.startswith('ATOM'):
                x_coord = float(line[30:38])
                y_coord = float(line[38:46])
                z_coord = float(line[46:54])

                if any((x_coord, y_coord, z_coord) == coord for coord in n_atom_coords):
                    residue_number = int(line[22:26])
                    residue_numbers.append(residue_number)
    return residue_numbers

# find_atom_indices
def find_atom_indices(comp_pdb, atom_coords):
    atom_indices = []

    with open(comp_pdb, 'r') as file:
        for line in file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                x_coord = float(line[30:38])
                y_coord = float(line[38:46])
                z_coord = float(line[46:54])

                if any((x_coord, y_coord, z_coord) == coord for coord in atom_coords):
                    atom_index = int(line[6:11])
                    atom_indices.append(atom_index)

    return atom_indices

# convert_to_csv
def convert_to_csv(input_file, output_file):
    # Read the data from the input file, skipping the first line
    with open(input_file, 'r') as infile:
        lines = infile.readlines()[1:]

    # Process the data and write to CSV
    data = [line.split() for line in lines]
    with open(output_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerows(data)

def convert_to_csv_title(input_file, output_file):
    # Read the data from the input file, skipping the first line
    with open(input_file, 'r') as infile:
        lines = infile.readlines()[0:]

    # Process the data and write to CSV
    data = [line.split() for line in lines]
    with open(output_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerows(data)

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

# calculate_mean
def calculate_mean(data, column):
    values = data[:, column].astype(float)
    return np.mean(values)

# calculate_sum
def calculate_sum(data, column):
    values = data[:, column].astype(float)
    return np.sum(values)

# write_properties
def write_properties_four(mean_rmsd, mean_surf, mean_contact, sum_hbond):
    with open('properties_raw1.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([mean_rmsd, mean_surf, mean_contact, sum_hbond])

# read_properties
def read_properties():
    try:
        with open('properties_raw1.csv', 'r') as file:
            reader = csv.reader(file)
            properties_data = list(reader)
        return np.array(properties_data)
    except FileNotFoundError:
        return None

# write_properties
def write_properties(data):
    with open('properties_raw1.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(data)
    
# Function to concatenate two CSV files
def concatenate_csv_horizontally(input_files, output_file):
    with open(output_file, 'w', newline='') as output_csv:
        writer = csv.writer(output_csv)
        
        # Write the headers from all input files
        headers = []
        for input_file in input_files:
            with open(input_file, 'r') as input_csv:
                header = next(csv.reader(input_csv))
                headers.extend(header)
        
        # Write the combined header
        writer.writerow(headers)
        
        # Write the content from all input files
        for input_file in input_files:
            with open(input_file, 'r') as input_csv:
                # Skip the header line
                next(input_csv)
                
                # Write the content
                for line in input_csv:
                    writer.writerow(line.strip().split(','))

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
        #print("Failed to extract nstlim and dt from the file.")
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

    #### 2.1 Preparation for CPPTRAJ
        
    # User needs to set this manually:
    c_salt = 0.15
    cation_name = "Na+" 
    anion_name = "Cl-" # c_salt, cation_name, anion_name must be set to the same time as FlexiSimScr.py

    resid = f"!:WAT,{cation_name},{anion_name}"
    distance = 4.5 

    # write cpp.in file
    data_csv_dir = 'data_csv'
    shutil.rmtree(data_csv_dir, ignore_errors=True)
    os.makedirs(data_csv_dir)

    mdrun_in_path = 'mdrun.in'  
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

    dt = interval * 0.1

    # DockPose mode
    if args.dockpose: 
        cpp = f'''
    parm ../comp_sol.prmtop
    trajin ../comp_sol.inpcrd
    reference ../comp_sol.inpcrd [DockPose]
    rms Protein ref [DockPose] ({resid})&(!@H=) out all_rmsd{stage}.csv 
    rms Ligand ref [DockPose] :1&(!@H=) out all_rmsd{stage}.csv nofit
    molsurf surf1 :1 probe 1 out all_surf{stage}.csv
    nativecontacts :1&!@H= (!:WAT,Na+,Cl-)&(!:1)&(!@H=) distance {distance} out all_contact{stage}.csv
    hbond donormask :1 acceptormask ({resid})&(!:1) avgout avghb2.csv out all_hbond{stage}.csv 
    hbond donormask ({resid})&(!:1) acceptormask :1 avgout avghb1.csv out all_hbond{stage}.csv''' 
     
    else:
        cpp = f'''
    parm ../comp_sol.prmtop
    trajin ../mdrun.nc {start_frame} {end_frame} {interval}
    reference ../comp_sol.inpcrd [DockPose]
    rms Protein ref [DockPose] ({resid})&(!@H=) time {dt} out all_rmsd{stage}.csv 
    rms Ligand ref [DockPose] :1&(!@H=) out all_rmsd{stage}.csv nofit
    molsurf surf1 :1 probe 1 time {dt} out all_surf{stage}.csv
    nativecontacts :1&!@H= (!:WAT,Na+,Cl-)&(!:1)&(!@H=) distance {distance} time {dt} out all_contact{stage}.csv
    hbond donormask :1 acceptormask ({resid})&(!:1) time {dt} avgout avghb1.csv out all_hbond{stage}.csv 
    hbond donormask ({resid})&(!:1) acceptormask :1 time {dt} avgout avghb1.csv out all_hbond{stage}.csv'''
        
    # User-defined contact part in cpp.in
    if args.contact:
        contact_info = args.contact.split(',')
        contact_pairs = [contact_info[i:i+3] for i in range(0, len(contact_info), 3)]
        number1 = 1
        
        for contact_pair in contact_pairs:
            chain1 = contact_pair[0]
            residue1 = int(contact_pair[1])
            atom1 = contact_pair[2]

            atom_coords1 = find_atom_coordinates('../protein_raw.pdb', chain1, residue1,atom1)

            print(f'{residue1}&{atom1}    ', end='', file=open('./data_csv/contact_title.csv', 'a'))

            atom_indices1 = find_atom_indices('comp_raw2.pdb',atom_coords1)
            atom_indices1_str = str(atom_indices1).strip('[]')

            cpp += f'''
    nativecontacts @{atom_indices1_str} :1&!@H= distance {distance} out contact{stage}_{number1}.csv'''
            
            number1 = number1 + 1

    # User-defined Hbond part in cpp.in
    if args.hbond:
        hbond_info = args.hbond.split(',')
        hbond_pairs = [hbond_info[i:i+3] for i in range(0, len(hbond_info), 3)]
        number2 = 1

        for hbond_pair in hbond_pairs:
            
            chain2 = hbond_pair[0]
            residue2 = int(hbond_pair[1])
            atom2 = hbond_pair[2]

            atom_coords2 = find_atom_coordinates('../protein_raw.pdb', chain2, residue2,atom2)

            print(f'{residue2}&{atom2}    ', end='', file=open('./data_csv/hbond_title.csv', 'a'))

            atom_indices3 = find_atom_indices('comp_raw2.pdb',atom_coords2)
            atom_indices3_str = str(atom_indices3).strip('[]')
            
            cpp += f'''
    hbond donormask @{atom_indices3_str} acceptormask :1 avgout avghb2.csv out hbond{stage}_{number2}.csv
    hbond donormask :1 acceptormask @{atom_indices3_str} avgout avghb2.csv out hbond{stage}_{number2}.csv'''
            number2 = number2 + 1
            
    # User-defined Hbond_res part in cpp.in
    if args.hbond_res:
        hbond_res_info = args.hbond_res.split(',')
        hbond_res_pairs = [hbond_res_info[i:i+2] for i in range(0, len(hbond_res_info), 2)]
        number4 = 1

        for hbond_res_pair in hbond_res_pairs:
            
            chain4 = hbond_res_pair[0]
            residue4 = int(hbond_res_pair[1])

            atom_coords4 = find_N_atom_coordinates('../protein_raw.pdb', chain4, residue4)

            print(f'{residue4}    ', end='', file=open('./data_csv/hbond_res_title.csv', 'a'))

            hbond_res_id4 = find_residue_number('comp_raw2.pdb',atom_coords4)
            hbond_res_id4_str = str(hbond_res_id4).strip('[]')
            
            cpp += f'''
    hbond donormask :{hbond_res_id4_str} acceptormask :1 avgout avghb3.csv out hbond_res{stage}_{number4}.csv
    hbond donormask :1 acceptormask @{hbond_res_id4_str} avgout avghb3.csv out hbond_res{stage}_{number4}.csv'''
            number4 = number4 + 1

    # continue cpp.in
    cpp += f'''
    go
    quit'''

    with open('cpp.in', 'w') as f_out:
        f_out.write(cpp)

    # Calculated the csv_data
    shutil.copy('cpp.in', os.path.join(data_csv_dir, 'cpp.in'))
    os.chdir(data_csv_dir)

    with open('cpp.log', 'w') as log_file:
        process = subprocess.Popen(['cpptraj', '-i', 'cpp.in'], stdout=log_file, stderr=subprocess.PIPE)

    process.communicate()

    print (f'{stage}.1 The output files of Cpp.in files has been generated.') 

    #### 2.2 Get properties: properties_raw1.csv
    # RMSD + SURF + CONTACT + HBOND
    convert_to_csv(f'all_rmsd{stage}.csv', 'all_rmsd.csv')
    convert_to_csv(f'all_surf{stage}.csv', 'all_surf.csv')
    convert_to_csv(f'all_contact{stage}.csv', 'all_contact.csv')
    convert_to_csv(f'all_hbond{stage}.csv', 'all_hbond.csv')

    all_rmsd_data = read_csv(f'all_rmsd.csv') 
    all_surf_data = read_csv(f'all_surf.csv')
    all_contact_data = read_csv(f'all_contact.csv')
    all_hbond_data = read_csv(f'all_hbond.csv')

    mean_rmsd_all = calculate_mean(all_rmsd_data, 2)
    mean_surf_all = calculate_mean(all_surf_data, 1)
    mean_contact_all = calculate_mean(all_contact_data, 2)
    mean_hbond_data1 = calculate_mean(all_hbond_data, 1)
    mean_hbond_data2 = calculate_mean(all_hbond_data, 2)
    sum_hbond_all = mean_hbond_data1 + mean_hbond_data2

    mean_rmsd_all = round(mean_rmsd_all, 2)
    mean_surf_all = round(mean_surf_all, 2)
    mean_contact_all = round(mean_contact_all, 2)
    sum_hbond_all = round(sum_hbond_all, 2)

    write_properties_four(mean_rmsd_all, mean_surf_all, mean_contact_all, sum_hbond_all)

    # ID
    id_data = stage
    
    properties_data = read_properties()

    if properties_data is not None:
        new_column = [f"{id_data}" for _ in range(len(properties_data))]
        properties_data = np.column_stack((properties_data, new_column))
    else:
        properties_data = [[f"{id_data}"]]

    write_properties(properties_data)   

    # Title
    with open('../LIG_raw.mol2', 'r') as file:
        lines = file.readlines()

    if len(lines) >= 2:
        title_data = lines[1].strip()
    else:
        print("Error: File 'LIG_raw.mol2' doesn't have at least two lines.")
    
    properties_data = read_properties()

    if properties_data is not None:
        new_column = [f"{title_data}" for _ in range(len(properties_data))]
        properties_data = np.column_stack((properties_data, new_column))
    else:
        properties_data = [[f"{title_data}"]]
    
    write_properties(properties_data)
    
    # User-defined contact
    if args.contact:
        for contact_i in range(1,number1): 

            convert_to_csv(f'contact{stage}_{contact_i}.csv', f'contact_{contact_i}.csv')
            contact_data_users = read_csv(f'contact_{contact_i}.csv')

            mean_contact_users = calculate_mean(contact_data_users, 2)
            mean_contact_users = round(mean_contact_users, 2)

            properties_data = read_properties()

            if properties_data is not None:
                new_column_contact = [f"{mean_contact_users}" for _ in range(len(properties_data))]
                properties_data = np.column_stack((properties_data, new_column_contact))
            else:
                properties_data = [[f"{mean_contact_users}"]]

            
            write_properties(properties_data)

    # User-defined Hbond
    if args.hbond:
        for hbond_i in range(1,number2):

            convert_to_csv(f'hbond{stage}_{hbond_i}.csv', f'hbond{stage}_{hbond_i}_raw.csv')
            hbond_data_user = read_csv(f'hbond{stage}_{hbond_i}_raw.csv')        
            properties_data = read_properties()

            mean_hbond_1_user = calculate_mean(hbond_data_user, 1)
            mean_hbond_2_user = calculate_mean(hbond_data_user, 2)
            sum_hbond_users = mean_hbond_1_user + mean_hbond_2_user
            sum_hbond_users = round(sum_hbond_users, 2)

            if properties_data is not None:
                new_column_hbond = [f"{sum_hbond_users}" for _ in range(len(properties_data))]
                properties_data = np.column_stack((properties_data, new_column_hbond))
            else:
                properties_data = [[f"{sum_hbond_users}"]]

            write_properties(properties_data)

    # User-defined Hbond_res
    if args.hbond_res:
        for hbond_res_i in range(1,number4):

            convert_to_csv(f'hbond_res{stage}_{hbond_res_i}.csv', f'hbond_res{stage}_{hbond_res_i}_raw.csv')
            hbond_res_data_user = read_csv(f'hbond_res{stage}_{hbond_res_i}_raw.csv')        
            properties_data = read_properties()

            mean_hbond_res_1_user = calculate_mean(hbond_res_data_user, 1)
            mean_hbond_res_2_user = calculate_mean(hbond_res_data_user, 2)
            sum_hbond_res_users = mean_hbond_res_1_user + mean_hbond_res_2_user
            sum_hbond_res_users = round(sum_hbond_res_users, 2)

            if properties_data is not None:
                new_column_hbond = [f"{sum_hbond_res_users}" for _ in range(len(properties_data))]
                properties_data = np.column_stack((properties_data, new_column_hbond))
            else:
                properties_data = [[f"{sum_hbond_res_users}"]]

            write_properties(properties_data)

    print (f'{stage}.2 The properties_raw1.csv files has been generated.') 

    #### 2.3 properties_raw.csv 
    source_csv_file = './properties_raw1.csv'
    destination_csv_file = '../../properties_raw.csv'

    with open(source_csv_file, 'r') as source_file:
        last_line = source_file.readlines()[-1]

    with open(destination_csv_file, 'a') as dest_file: # Append the last line to the destination CSV file
        dest_file.write(last_line)

    print(f"{stage}.3 Data has been appended to properties_raw.csv.")

    #### 2.4 Go back to the parent directory. ####
    os.chdir(complex_md_dir)

#### 3. get properties: properties.csv ####
# Get titles: titles_raw1.csv 
if args.contact:
    source_path1 = f"{complex_md_dir}/lig1/data_csv/contact_title.csv"
    destination_path1 = f"{complex_md_dir}/contact_title_raw.csv"
    shutil.copyfile(source_path1, destination_path1)
    contact_title_data = read_csv('contact_title_raw.csv')

if args.hbond:
    source_path2 = f"{complex_md_dir}/lig1/data_csv/hbond_title.csv"
    destination_path2 = f"{complex_md_dir}/hbond_title_raw.csv"
    shutil.copyfile(source_path2, destination_path2)
    hbond_title_data = read_csv('hbond_title_raw.csv')
    contact_values = [item.strip() for sublist in contact_title_data for item in sublist[0].split()]
    hbond_values = [item.strip() for sublist in hbond_title_data for item in sublist[0].split()]

if args.hbond_res:
    source_path3 = f"{complex_md_dir}/lig1/data_csv/hbond_res_title.csv"
    destination_path3 = f"{complex_md_dir}/hbond_res_title_raw.csv"
    shutil.copyfile(source_path3, destination_path3)
    hbond_res_title_data = read_csv('hbond_res_title_raw.csv')
    hbond_res_values = [item.strip() for sublist in hbond_res_title_data for item in sublist[0].split()]

title_raw_path = 'title_raw.csv'
if os.path.exists(title_raw_path):
    title_raw_data = read_csv(title_raw_path)
else:
    title_raw_data = []

# title of row
mean_rmsd = 'Mean_RMSD'
mean_surf = 'Mean_Surf'
mean_contact = 'Mean_Contact_Count'
mean_hbond = 'Mean_Hbond_Count'
id_row = 'ID'
title_row = 'Title'

row_to_append = [mean_rmsd, mean_surf, mean_contact, mean_hbond, id_row, title_row] 

if args.contact:
    row_to_append = row_to_append + contact_values

if args.hbond:
        row_to_append = row_to_append + hbond_values  

if args.hbond_res:
        row_to_append = row_to_append + hbond_res_values 

title_raw_data.append(row_to_append)

write_csv(title_raw_path, title_raw_data) # Write to title_raw.csv

# properties_raw.csv + title_raw.csv --> properties.csv
title_raw_path = 'title_raw.csv'
properties_raw_path = 'properties_raw.csv'
properties_path = 'properties.csv'

# Read the content of both files
with open(title_raw_path, 'r') as title_raw_file:
    title_raw_content = title_raw_file.read()

with open(properties_raw_path, 'r') as properties_raw_file:
    properties_raw_content = properties_raw_file.read()

# Concatenate the content and write to the output file
with open(properties_path, 'w') as properties_file:
    properties_file.write(title_raw_content + properties_raw_content)

# Delete raw files.
pattern = '*raw.csv'
[os.remove(file) and print(f"Deleted: {file}") for file in glob.glob(pattern)]

print("All data has been saved in properties.csv.")

###################################################################################################
# List of Scripts for Flexible Simulation and Stable Analysis Post Screening
#
#1. FlexiSimScr.py:     Flexible Simulation for Screening
#2. FlexiSimComp.py:    Flexible Simulation for complex
#3. FlexiSimPro.py:     Flexible Simulation for protein
#4. FlexiSimStab.py:    Calculate contact/stablity after Docking/MD simulation
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


# command in current folder:
#example1:
# nohup python FlexiSimStab.py -num 1 -t 5 -contact Z,694,NE2,Z,694,NE2 -hbond Z,694,NE2,Z,694,NE2 -hbond_res Z,694
#example2:
# nohup python FlexiSimStab.py -dockpose -contact R,329,ND2,R,121,OD2 -hbond R,329,ND2,R,121,OD2 -hbond_res R,329,R,121 

