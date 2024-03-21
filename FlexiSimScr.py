#!/usr/bin/env python
import argparse
import os
import shutil
import subprocess
from rdkit import Chem

###################################################################################################
#print('FlexiSimScr: \nUsage: FlexiSimScr.py '
#      '[-pro ] <Protein in .pdb format, e.g. protein.pdb> '
#      '[-lig ] <Ligand in .mol2 or .sdf format, e.g. top_molecules.mol2, preference is given to mol2 input>'
#      '[-time] <Simulation time, e.g. 100, defaulting to 100 (ns)> '
#      '[-gpuid]    <GPUID, defaulting to 0>')
#      '[-resid ]   <constrained residue, e.g. 1-200, defaulting to residues of protein & ligand>'
#      '[-workdir]  <Work directory, e.g. /home/yangli02, defaulting to the current directory> '
#      '[-ssbond]   <Residues for forming Disulfide bonds, e.g. B,71,B,85>
#print ('Command: \npython FlexiSimScr.py -pro protein.pdb -lig top_molecules.mol2'
#      ' -ssbond B,71,B,85 -num 1 -time 100 -gpuid 0 &')
###################################################################################################
       
#### 1. Parse input line ####
parser = argparse.ArgumentParser()
parser.add_argument("-pro", type=str, required=True, help="A protein file in PDB format.") 
parser.add_argument("-lig", type=str, required=True,help="Ligand files in MOL2 or SDF format. Preference is given to mol2 input.")
parser.add_argument("-num", type=str, required=False,help="Number of ligands, default set to the total number of small molecules.")
parser.add_argument("-time", type=str, required=False, help="Simulation time", default="100")
parser.add_argument("-gpuid", type=str, required=False, help="GPU ID", default="0")
parser.add_argument("-resid", type=str, required=False, help="Constrained residues, defaulting to residues of protein & ligand.")
parser.add_argument("-workdir", type=str, required=False, help="Work directory, defaulting to the current directory.")
parser.add_argument("-ssbond", type=str, help="Residue information in the format: chain1,residue_number1,chain2,residue_number2")
args = parser.parse_args()

# Number setting part1
num = int(args.num) if args.num is not None else 1

# Time setting
time = int(args.time) if args.time else None
if time is None or time <= 0:
    print('Simulation time is not provided or invalid, setting to default: 100 ns.')
    time = 100

# GPUID setting
GPUID = int(args.gpuid) if args.gpuid else None
if GPUID is None:
    print('GPUID is not provided, setting to default: 0.')
    GPUID = 0

# Constrained residues
if args.resid is None:
    print('Constrained residues is not provided')

# Work directory setting
workdir = args.workdir if args.workdir else None

if workdir is None:
    print('Work directory not specified, setting to current directory.')
    workdir = os.path.abspath(os.getcwd())


#### 2. Work directory & Preparation ####
###################################################################################################
# separate_mol2: Save each molecule as a separate .sdf file.
def separate_mol2(input_file): 
    global number_lig

    with open(input_file, 'r') as file:
        content = file.read().split('@<TRIPOS>MOLECULE')[1:] 

    number_lig = len(content)

    for i, mol in enumerate(content, start=1):
        
        mol_block = '@<TRIPOS>MOLECULE\n' + mol.strip()      
        mol2_filename = f'lig{i}.mol2'

        with open(mol2_filename, 'w') as output_file:
            output_file.write(mol_block)

# separate_sdf
def separate_sdf(input_file):
    global number_lig

    with open(input_file, 'r') as file:
        content = file.read().split('$$$$\n')[:-1] 

    number_lig = len(content)

    for i, mol in enumerate(content, start=1):

        # Reconstruct the molecule block including the footer
        mol_block = mol.strip() + '$$$$\n'
        sdf_filename = f'lig{i}.sdf'

        # Save each molecule as a separate .sdf file
        with open(f'lig{i}.sdf', 'w') as output_file:
            output_file.write(mol_block)

        # convert to .mol2 file
        subprocess.run(['obabel', '-isdf', sdf_filename, '-omol2', '-O', f'lig{i}.mol2'])

    # assign_chain_if_missing

# assign_chain_if_missing
def assign_chain_if_missing(protein_pdb):
    has_chain = False
    with open(protein_pdb, 'r') as file:
        lines = file.readlines()

    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            if line[21].strip():  # Checks if the chain field is not empty
                has_chain = True
                break

    if not has_chain:
        new_lines = []
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                new_line = line[:21] + 'Z' + line[22:]
                new_lines.append(new_line)
            else:
                new_lines.append(line)

        with open(protein_pdb, 'w') as file:
            file.writelines(new_lines)
            print("Protein chain set to 'Z'.")

# remove_chain_info
def remove_chain_info(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    updated_lines = []
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            updated_line = line[:72] + ' ' + line[73:]
            updated_lines.append(updated_line)
        else:
            updated_lines.append(line)

    with open(file_path, 'w') as file:
        file.writelines(updated_lines)

###################################################################################################

os.chdir(workdir)
complex_md_dir = os.path.join(os.getcwd(), 'complex_md')

if not os.path.exists(complex_md_dir):
    os.makedirs(complex_md_dir)

protein_path = os.path.join(workdir, args.pro)
shutil.copy(protein_path, os.path.join(complex_md_dir, args.pro))
ligand_path = os.path.join(workdir, args.lig)
shutil.copy(ligand_path, os.path.join(complex_md_dir, args.lig))
os.chdir(complex_md_dir)
subprocess.run(['/bin/bash', '-c', 'source /home/soft/amber22/amber.sh'])

#Save each molecule as a separate .mol2 file.
if args.lig.endswith('.mol2'):
    separate_mol2(args.lig) # from the '.mol2' file

if args.lig.endswith('.sdf'):
    separate_sdf(args.lig) # from the '.sdf' file

if args.num is None:
    num = number_lig
    print('Ligand number is not specified, setting to default: Total number of ligands in file.')

# Preparation for protein
with open(str(args.pro), 'r') as file:
    lines = file.readlines()

filtered_lines = [line for line in lines if 'CONECT' not in line and 'MASTER' not in line]

modified_lines = [line.replace('END', 'TER') for line in filtered_lines]
modified_lines = [line.replace('CA  NMA', 'C   NME') for line in modified_lines]
modified_lines = [line.replace('NMA', 'NME') for line in modified_lines]

final_lines = [line for line in modified_lines if line.startswith(('ATOM', 'TER', 'END'))]

with open('protein_raw.pdb', 'w') as file:
    file.writelines(final_lines)

file_path = 'protein_raw.pdb' # change the residue number

with open(file_path, 'r') as file:
    lines = file.readlines()

new_lines = []

for line in lines:

    if line.startswith('ATOM'):
        atom_number = int(line[6:11].strip())
        residue_name = line[17:20].strip()
        chain_id = line[21:22]
        residue_number = int(line[22:26].strip())
        insert_code = line[26:27]

        if residue_name == 'NME' and insert_code == 'A':
            new_residue_number = residue_number + 1
            new_line = f"{line[:22]}{new_residue_number:4} {line[27:]}"  
            new_lines.append(new_line)
        else:
            new_lines.append(line)
    else:
        new_lines.append(line)

with open(file_path, 'w') as file:
    file.writelines(new_lines)

# assign Z for chain if missing
assign_chain_if_missing('protein_raw.pdb')
remove_chain_info('protein_raw.pdb')

print ('File preparation has been completed.') # check 1 #

#### 3. Performing calculations individually for each ligand ####

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

# replace_residue_names
def replace_residue_names(pdb_file, residue_numbers, new_residue_name):

    new_lines = []

    with open(pdb_file, 'r') as file:

        for line in file:

            if line.startswith('ATOM') and int(line[22:26]) in residue_numbers and line[17:20] == 'CYS':
                line = line[:17] + new_residue_name + line[20:]
            new_lines.append(line)

    with open(pdb_file, 'w') as file:
        file.writelines(new_lines)
###################################################################################################
 
for stage in range(1, int(num) + 1):
  
    stage = str(stage)
    lig_dir = f'lig{stage}'
    os.makedirs(lig_dir, exist_ok=True)

    shutil.copy(f'lig{stage}.mol2', f'{lig_dir}/LIG_raw.mol2')
    os.chdir(lig_dir)
    
    # To create LIG.lib and LIG.frcmod
    subprocess.run(['/home/soft/amber22/bin/antechamber', '-i', 'LIG_raw.mol2', '-fi', 'mol2', '-o', 'LIG_raw1.mol2', '-fo', 'mol2', '-c', 'bcc', '-s', '2', '-at', 'gaff2', '-dr', 'no'], stdout=open('antechamber.log', 'w'), stderr=subprocess.PIPE)
    subprocess.run(['/home/soft/amber22/bin/parmchk2', '-i', 'LIG_raw1.mol2', '-f', 'mol2', '-o', 'LIG.frcmod'])

    with open('LIG_raw1.mol2', 'r') as file: 
        lines = file.readlines()

    filtered_lines = [line for line in lines if 'CONECT' not in line and 'MASTER' not in line]
    modified_lines = [line.replace(' UNL', ' LIG').replace(' MOL', ' LIG').replace(' UNK', ' LIG').replace('CL', 'Cl').replace('BR', 'Br') for line in filtered_lines]

    with open('LIG_raw2.mol2', 'w') as file:
        file.writelines(modified_lines)

    # To create LIG_raw.pdb
    leap_lig = '''
    source leaprc.protein.ff19SB
    source leaprc.gaff2
    LIG = loadmol2 LIG_raw2.mol2
    loadamberparams LIG.frcmod
    saveoff LIG LIG.lib
    saveamberparm LIG LIG.prmtop LIG.inpcrd
    savepdb LIG LIG_raw.pdb
    quit
    '''

    with open('leap_lig.in', 'w') as f_out:
        f_out.write(leap_lig)

    subprocess.run(['/home/soft/amber22/bin/tleap', '-s', '-f', 'leap_lig.in'], stdout=open('leap_lig.log', 'w'), stderr=subprocess.PIPE)
    print(f'{stage}.1 Ligand{stage} preparation has been completed.') # check 2 #
    
    # To create the top&crd file for compound
    with open('LIG_raw.pdb', 'r') as file: 
        lines = file.readlines()

    modified_lines = [line.replace('END', 'TER').replace('CL', 'Cl').replace('BR', 'Br') for line in lines if 'CONECT' not in line]
    
    with open('LIG_raw1.pdb', 'w') as file: # related to the ligand
        file.writelines(modified_lines)
        
    with open('LIG_raw1.pdb', 'r') as lig_file, open('../protein_raw.pdb', 'r') as protein_file, open('comp_raw.pdb', 'w') as output_file:
        lig_lines = lig_file.readlines()

        output_file.writelines(lig_lines)
        output_file.write('\n')  

        for line in protein_file:
            output_file.write(line)

    with open('comp_raw.pdb', 'r') as file: # related to the compound
        lines = file.readlines()

    filtered_lines = [line for line in lines if line.strip() and line.strip().split()[-1] != 'H']

    with open('comp_raw.pdb', 'w') as file:
        file.writelines(filtered_lines)
    
    with open('pdb4amb.log', 'w') as log_file:
        subprocess.run(['/home/soft/amber22/bin/pdb4amber', '-y', '-i', 'comp_raw.pdb', '-o', 'comp_raw1.pdb', '-l', 'pdb4amb.log'], stdout=log_file)

    complex = "comp_raw1.pdb"
   
    print(f'{stage}.2 Protein{stage} preparation has been completed.') # check 3 #

    # add water
    leap_cal = '''\
    source leaprc.protein.ff19SB
    source leaprc.gaff2
    source leaprc.water.opc
    loadamberparams frcmod.ionslm_1264_opc
    loadamberparams LIG.frcmod
    loadoff LIG.lib
    com = loadpdb {}
    solvatebox com OPCBOX 12
    charge com
    quit
    '''.format(complex)

    with open('leap_cal.in', 'w') as f_out:
        f_out.write(leap_cal)

    with open('leap_cal.log', 'w') as log_file:
        process = subprocess.Popen(['/home/soft/amber22/bin/tleap', '-s', '-f', 'leap_cal.in'], stdout=log_file, stderr=subprocess.PIPE)

    process.communicate()

    # cal water num
    with open('leap_cal.log', 'r') as log_file:
        lines = log_file.readlines()

    residues_line = None
    charge_line = None

    for line_number, line in enumerate(lines, start=1):

        if 'residues' in line:
            residues_line = line_number
        elif 'Total perturbed charge' in line:
            charge_line = line_number

    # parameters setting
    c_salt = 0.15
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

    if residues_line and charge_line:
        wat_num = int(lines[residues_line - 1].split()[1]) 
        charge = float(lines[charge_line - 1].split()[3])  
        cation_num = int(wat_num * c_salt * 18 / 1000) 
        anion_num = int(wat_num * c_salt * 18 / 1000 + charge) # parameters setting: If multivalent cations are used, multiply the valence state by the wat_num.

    # create the top&crd file
    leap_comp_sol = f'''
    source leaprc.protein.ff19SB
    source leaprc.gaff2
    source leaprc.water.opc
    loadamberparams frcmod.ionslm_1264_opc
    loadamberparams LIG.frcmod
    loadoff LIG.lib
    com = loadpdb {complex}
    '''

    # ssbond
    if args.ssbond:
        ssbond_info = args.ssbond.split(',')
        ssbond_pairs = [ssbond_info[i:i+4] for i in range(0, len(ssbond_info), 4)]

        for ssbond_pair in ssbond_pairs:
            chain1 = ssbond_pair[0]
            residue1 = int(ssbond_pair[1])
            chain2 = ssbond_pair[2]
            residue2 = int(ssbond_pair[3])
        #   print(chain1, residue1,chain2,residue2)
            n_atom_coords1 = find_N_atom_coordinates('../protein_raw.pdb', chain1, residue1)
            n_atom_coords2 = find_N_atom_coordinates('../protein_raw.pdb', chain2, residue2)
        #   print(n_atom_coords1,n_atom_coords2)
            residue_numbers1 = find_residue_number(complex, n_atom_coords1)
            residue_numbers2 = find_residue_number(complex, n_atom_coords2)
        #   print(residue_numbers1,residue_numbers2)
            replace_residue_names('comp_raw1.pdb', residue_numbers1, 'CYX')
            replace_residue_names('comp_raw1.pdb', residue_numbers2, 'CYX')

        #   Add bond information to leap_comp_sol
            leap_comp_sol += f"bond com.{residue_numbers1[0]}.SG com.{residue_numbers2[0]}.SG\n"

    leap_comp_sol += f'''savepdb com comp_raw2.pdb
    solvatebox com OPCBOX 12
    addions com {cation_name} {cation_num}
    addions com {anion_name} {anion_num}
    saveamberparm com comp_sol.prmtop comp_sol.inpcrd
    savepdb com comp_sol.pdb
    quit
    '''

    with open('leap_comp_sol.in', 'w') as f_out:
        f_out.write(leap_comp_sol)

    with open('leap_comp_sol.log', 'w') as log_file:
        process = subprocess.Popen(['/home/soft/amber22/bin/tleap', '-s', '-f', 'leap_comp_sol.in'], stdout=log_file, stderr=subprocess.PIPE)

    process.communicate()
    
    print(f'{stage}.3 Complex{stage} preparation has been completed.') # check 4 #

#### 4. Establishing the MD input file  ####
    
###################################################################################################
    def write_input_files(file_contents):

        for file_name, content in file_contents.items():

            with open(file_name, 'w') as f_out:
                f_out.write(content)
###################################################################################################
                
    # Constrain residues setting
    resid = int(args.resid) if args.resid else None

    if resid is None:
        print('Constrain residues not provided, setting to default: protein & ligand')
        resid = f"!:WAT,{cation_name},{anion_name}"

    # parameters setting
    step_value = int(time) * 500000
    cut_value = 12.0
    temp0_value = 310.0
    dt_value = 0.002
    ntpr_value = 5000
    ntwx_value = 5000
    ntwr_value = 5000
    
    file_contents = {
        'min.in': (
            "minimise ras-raf\n &cntrl\n  imin=1,maxcyc=10000,ncyc=5000,\n"
            f"  cut={cut_value},ntb=1,ntc=1,ntf=1,ntpr=500,\n /\n"
        ),
        'heat.in': (
            "heat 100k,310k\n &cntrl\n  imin=0,irest=0,ntx=1,nstlim=1000000,dt=0.001,\n"
            f"  ig=-1,ntc=1,ntf=1,cut={cut_value},ntb=1,ntpr=5000,\n"
            f"  ntwx=5000,ntwr=5000,ntt=3,gamma_ln=2.0,tempi=0.0,temp0={temp0_value},\n"
            f"  ntr=1,restraintmask='{resid}',restraint_wt=50.0,\n /\n"
            " &wt TYPE='TEMP0', istep1=0, istep2=250000,\n  value1=0.1, value2=100.0, /\n"
            f" &wt TYPE='TEMP0', istep1=250001, istep2=500000,\n  value1=100.0, value2={temp0_value}, /\n"
            " &wt TYPE='END' /\n"
        ),
        'density1.in': (
            f"NVT\n &cntrl\n  imin=0,irest=0,ntx=1,nstlim=150000,dt={dt_value},\n"
            f"  iwrap=1,ntc=2,ntf=2,cut={cut_value}, ntb=2, ntp=1, taup=1.0,\n"
            f"  ntpr={ntpr_value},ntwx={ntwx_value},ntwr={ntwr_value},ntt=3, gamma_ln=2.0,\n"
            f"  temp0={temp0_value},ntr=1,restraintmask='{resid}',restraint_wt=50.0,\n /\n"
        ),
        'density2.in': (
            f"NVT\n &cntrl\n  imin=0,irest=0,ntx=1,nstlim=350000,dt={dt_value},\n"
            f"  iwrap=1,ntc=2,ntf=2,cut={cut_value}, ntb=2, ntp=1, taup=1.0,\n"
            f"  ntpr={ntpr_value},ntwx={ntwx_value},ntwr={ntwr_value},ntt=3, gamma_ln=2.0,\n"
            f"  temp0={temp0_value},ntr=1,restraintmask='{resid}',restraint_wt=50.0,\n /\n"
        ),
        'equil1.in': (
            f"NVT\n &cntrl\n  imin=0,irest=0,ntx=1,nstlim=1000000,dt={dt_value}\n"
            f"  iwrap=1,ntc=2,ntf=2,cut={cut_value}, ntb=2, ntp=1, taup=2.0,\n"
            f"  ntpr={ntpr_value},ntwx={ntwx_value},ntwr={ntwr_value},ntt=3, gamma_ln=2.0,\n"
            f"  temp0={temp0_value},ntr=1,restraintmask='{resid}',restraint_wt=50.0,\n /\n"
        ),
        'equil2.in': (
            f"NVT\n &cntrl\n  imin=0,irest=0,ntx=1,nstlim=500000,dt={dt_value},\n"
            f"  iwrap=1,ntc=2,ntf=2,cut={cut_value}, ntb=2, ntp=1, taup=2.0,\n"
            f"  ntpr={ntpr_value},ntwx={ntwx_value},ntwr={ntwr_value},ntt=3, gamma_ln=2.0,\n"
            f"  temp0={temp0_value},ntr=1,restraintmask='{resid}',restraint_wt=5.0,\n /\n"
        ),
        'equil3.in': (
            f"NVT\n &cntrl\n  imin=0,irest=0,ntx=1,nstlim=500000,dt={dt_value},\n"
            f"  iwrap=1,ntc=2,ntf=2,cut={cut_value}, ntb=2, ntp=1, taup=2.0,\n"
            f"  ntpr={ntpr_value},ntwx={ntwx_value},ntwr={ntwr_value},ntt=3, gamma_ln=2.0,\n"
            f"  temp0={temp0_value},ntr=1,restraintmask='{resid}',restraint_wt=0.5,\n /\n"
        ),
        'mdrun.in': (
            f"Prod\n &cntrl\n  imin=0,irest=0,ntx=1,nstlim={step_value},dt={dt_value},\n"
            f"  iwrap=1,ntc=2,ntf=2,cut={cut_value}, ntb=2, ntp=1, taup=2.0,\n"
            "  ntpr=5000, ntwx=5000,ntwr=5000,ntt=3, gamma_ln=2.0,\n"
            f"  temp0={temp0_value}, ig=-1,\n /\n"
        )
    }

    write_input_files(file_contents)
    print(f'{stage}.4 Complex{stage} input files have been generated.') # check 5#

    # mdrun bashfile 
    mdrun_commands = [
        f'pmemd.cuda -O -i min.in -o min.out -p comp_sol.prmtop -c comp_sol.inpcrd -r min.rst -ref ref.inpcrd',
        f'pmemd.cuda -O -i heat.in -o heat.out -p comp_sol.prmtop -c min.rst -r heat.rst -x heat.nc -ref min.rst',
        f'pmemd.cuda -O -i density1.in -o density1.out -p comp_sol.prmtop -c heat.rst -r density1.rst -x density1.nc -ref heat.rst',
        f'pmemd.cuda -O -i density2.in -o density2.out -p comp_sol.prmtop -c density1.rst -r density2.rst -x density2.nc -ref density1.rst',
        f'pmemd.cuda -O -i equil1.in -o equil1.out -p comp_sol.prmtop -c density2.rst -r equil1.rst -x equil1.nc -ref density2.rst',
        f'pmemd.cuda -O -i equil2.in -o equil2.out -p comp_sol.prmtop -c equil1.rst -r equil2.rst -x equil2.nc -ref equil1.rst',
        f'pmemd.cuda -O -i equil3.in -o equil3.out -p comp_sol.prmtop -c equil2.rst -r equil3.rst -x equil3.nc -ref equil2.rst',
        f'pmemd.cuda -O -i mdrun.in  -o mdrun.out  -p comp_sol.prmtop -c equil3.rst -r mdrun.rst  -x mdrun.nc  -ref equil3.rst'
    ]

    with open('mdrun.sh', 'w') as f_out:
        f_out.write('source /home/soft/amber22/amber.sh\n')
        f_out.write(f'export CUDA_VISIBLE_DEVICES={GPUID}\n')

        for mdrun_command in mdrun_commands:
            f_out.write(f'{mdrun_command}\n')

    subprocess.run(['chmod', '+x', "mdrun.sh"])


#    mdrun_script_path = 'mdrun.sh'
#    subprocess.run(['nohup', 'sh', mdrun_script_path, '&'], check=True)

    # Wait for the completion of the mdrun.sh script
#    subprocess.run(['wait'], check=True)
#    print (f'{stage}.5 MD simulation for complex{stage} has been completed.') # check 6 #
    
    os.chdir(complex_md_dir)

print ('MD simulation for complex has been completed.')
exit()

###################################################################################################
# List of Scripts for Flexible Simulation and Stable Analysis Post Screening
#
#1. FlexiSimScr.py:     Flexible Simulation for Screening
#2. FlexiSimComp.py:    Flexible Simulation for complex
#3. FlexiSimPro.py:     Flexible Simulation for protein
#4. FlexiSimCont.py:    Calculate contact/stablity after Docking/MD simulation
#5. FlexiSimPBSA.py:    Calculate binding free energy after Docking/MD simulation

###################################################################################################
#### 
#1## The "protein_for_FlexiSimScr.pdb" file contains proteins without chain information, which is a special case. 
#### Proteins from PDB bank or processed by Schrödinger are suitable.
#
#2## The "ligands_for_FlexiSimScr.mol2" file includes three small molecules with post-docking conformations, each with different charges for testing.
#### Mol2 files exported from Schrödinger's docking results in the pv file can be used here.
#
#3## The "ssbond" is provided here for testing. If the input protein lacks chain information, use 'Z' as a placeholder. 
#### If chain information is available, use the corresponding chain identifier.
#
#4## The code for conducting MD run on 'line 550' has been commented out in advance for use during testing or when exclusively analyzing docking modes. 
#### It is recommended for users to uncomment and submit the task after testing.

# command in current folder:
# example1:
# nohup python FlexiSimComp.py -pro protein_for_FlexiSimScr.pdb -lig ligands_for_FlexiSimScr.mol2 -time 5 -gpuid 0 --ssbond Z,1320,Z,1082 &
# example2:
# nohup python FlexiSimScr.py -pro protein.pdb -lig ligands.mol2 -time 5 -gpuid 0 &