# vasp2pdb
Usage: vasp2pdb.py [POSCAR] [SETTINGS] [OUTPUT]

This program converts VASP POSCAR files to PDB (protain data bank) format.

Atoms in POSCAR are split into molecules (molecules are defined by LT-format file. LT-format file can be generated using moltemplate software https://www.moltemplate.org/) and then converted to PDB format.

Each line of SETTINGS file should contain pair of (.lt filename, number of that molecule).
Example:
  ethanol.lt 4
  water.lt 4

this file denotes that there are 4 ethanol molecules and 4 water molecules in the cell.
