import sys
import numpy as np

lammps_dump_file = open("output_ef.dump", "r")

line1 = ''
lattice = [0,0,0]
pos_x = []
pos_y = []
pos_z = []
f_x = []
f_y = []
f_z = []
epa = []
energy = 0
idx = 0
for line in lammps_dump_file:
    if 'ITEM: NUMBER' in line:
        line1 = 'ITEM: NUMBER'
        continue
    if 'ITEM: NUMBER' in line1:
        cfg_size = int(line)
        line1 = ''
    if 'ITEM: BOX' in line:
        line1 = 'ITEM: BOX'
        continue
    if 'ITEM: BOX' in line1:
        idx1 = 0
        for x in line.split():
            if (idx1 == 1):
                lattice[idx] = float(x)
            idx1 += 1
        idx += 1
        if (idx == 3): 
            line1 = ''
            idx = 0
        continue
    if 'ITEM: ATOMS' in line:
        line1 = 'ITEM: ATOMS'
        continue
    if 'ITEM: ATOMS' in line1:
        idx1 = 0
        for x in line.split():
            if (idx1 == 2): 
                pos_x.append(float(x))
            if (idx1 == 3): 
                pos_y.append(float(x))
            if (idx1 == 4): 
                pos_z.append(float(x))
            if (idx1 == 5): 
                f_x.append(float(x))
            if (idx1 == 6): 
                f_y.append(float(x))
            if (idx1 == 7): 
                f_z.append(float(x))
            if (idx1 == 8):
                epa.append(float(x))
            idx1 += 1
        idx += 1
        if (idx == cfg_size):
            line1 = ''
            
lammps_dump_file.close()
 
mlip_configuration_efs = open("output_ef.cfg", "w")
 
mlip_configuration_efs.write("BEGIN_CFG\n")
mlip_configuration_efs.write(" Size\n")
mlip_configuration_efs.write("    "+str(cfg_size)+"\n")
mlip_configuration_efs.write(" Supercell\n")
mlip_configuration_efs.write("    "+str(lattice[0])+"    0    0\n")
mlip_configuration_efs.write("    0    "+str(lattice[1])+"    0\n")
mlip_configuration_efs.write("    0    0    "+str(lattice[2])+"\n")
mlip_configuration_efs.write(" AtomData:  id type    cartes_x    cartes_y    cartes_z    fx    fy    fz\n")
for i in range (0,cfg_size):
    mlip_configuration_efs.write("    "+str(i+1)+"    0    "+str(pos_x[i])+"    "+str(pos_y[i])+"    "+str(pos_z[i])+"    "+str(f_x[i])+"    "+str(f_y[i])+"    "+str(f_z[i])+"\n")
    energy += epa[i]
mlip_configuration_efs.write(" Energy\n")
mlip_configuration_efs.write("    "+str(energy)+"\n")
mlip_configuration_efs.write("END_CFG\n")
mlip_configuration_efs.write("\n")

mlip_configuration_efs.close()
