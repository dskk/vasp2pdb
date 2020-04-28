import sys
import numpy as np
from copy import deepcopy
from operator import attrgetter
from pymatgen.io.vasp.inputs import Poscar

class UnionFind:
    def __init__(self,n):
        self.parent=[-1 for i in range(n)]
    def root(self,x):
        if self.parent[x]<0:
            return x
        else:
            self.parent[x]=self.root(self.parent[x])
            return self.parent[x]
    def unite(self, x, y):
        x=self.root(x)
        y=self.root(y)
        if x!=y:
            if(self.parent[y]<self.parent[x]):
                x,y=y,x
            self.parent[x]+=self.parent[y]
            self.parent[y]=x
            return False
        return True
    def same(self, x, y):
        return self.root(x)==self.root(y)
    def size(self, x):
        return -self.parent[self.root(x)]

class Atom:
    def __init__(self):
        self.position=None #np.array
        self.element=None #str e.g. "C"
        self.label=None #str e.g. "C11"
        self.degree=0 #int

class MolGraph:
    def __init__(self):
        self.atom_list=[] # list of Atom()
        self.label_to_ind_dic={}
        self.graph=None # graph[i][j]==1 if a bond exists between atom_i(==atom_list[i]) and atom_j. ==0 elsewise.
        self.atom_order=[]
    def dump(self):
        n=len(self.atom_list)
        print("Label        X              Y              Z        Deg  ",end="")
        for j in range(n):
            print(j%10,end=" ")
        print()
        for i in range(n):
            a=self.atom_list[i]
            if a.label==None:
                print("None    ",end="")
            else:
                print(f"{a.label: <8}",end="")
            print(f"{a.position[0]:10.5f}",end="     ")
            print(f"{a.position[1]:10.5f}",end="     ")
            print(f"{a.position[1]:10.5f}",end="     ")
            print(f"{a.degree: <4}",end="")
            for j in range(n):
                print("#" if self.graph[i][j] else "･",end=" ")
            print(f"{i%10}")

def permtest(k, mg1, mg2):
    """
        mg1: Test MolGraph. Whilc backtracking, ordering of mg1.atom_list will be changed.
        mg2: Reference MolGraph. Will not be changed.
        k: Subgraph made of first k atoms in mg1.atom_list is isomorphic to that of mg2.
    """
    n=len(mg1.atom_list)
    if k>=n:
        for i in range(n):
            mg1.atom_list[i].label=mg2.atom_list[i].label
            mg1.atom_order=[i for i in mg2.atom_order]
        return True
    for i in range(k,n):
        # mg1.atom_list[:k] is already constructed. Pick a candidate (mg.atom_list[i]) from mg1.atom_list[k:] and put it at [k].
        if mg1.atom_list[i].degree!=mg2.atom_list[k].degree or mg1.atom_list[i].element!=mg2.atom_list[k].element: #degree and element must be equal.
            continue
        mg1.atom_list[k],mg1.atom_list[i]=mg1.atom_list[i],mg1.atom_list[k] #Swap atom_k and atom_i
        for j in range(n): #Swap atom_k and atom_i
            mg1.graph[j][k],mg1.graph[j][i]=mg1.graph[j][i],mg1.graph[j][k] #Swap atom_k and atom_i
        mg1.graph[k],mg1.graph[i]=mg1.graph[i],mg1.graph[k] #Swap atom_k and atom_i
        isok=True
        for j in range(k+1):
            if mg1.graph[k][j] != mg2.graph[k][j]: #for directed graphs, check if mg1.graph[j][k]!=mg2.graph[j][k]
                isok=False
                break
        if isok and permtest(k+1,mg1,mg2):
            return True
        #atom_i could not satisfy constraint. Swap again to bring it back to the initial state
        #Swap atom_k and atom_i
        mg1.atom_list[k],mg1.atom_list[i]=mg1.atom_list[i],mg1.atom_list[k]
        for j in range(n):
            mg1.graph[j][k],mg1.graph[j][i]=mg1.graph[j][i],mg1.graph[j][k]
        mg1.graph[k],mg1.graph[i]=mg1.graph[i],mg1.graph[k]
    return False

def isomorph(mg1, mg2):
    n=len(mg1.atom_list)
    if n!=len(mg2.atom_list):
        return False
    return permtest(0, mg1, mg2) # return if labels are successfully assigned

if len(sys.argv)!=4:
    print("Usage: vasp2pdb.py [POSCAR] [SETTINGS] [OUTPUT]")
    sys.exit(0)
poscar_file=sys.argv[1]
setting_file=sys.argv[2]
out_file=sys.argv[3]

bond_count=dict()
bonds=[]
lt_molgraph_list=[]
with open(setting_file) as f:
    for line in f:
        filename, num_molecule = line.split()
        num_molecule = int(num_molecule)
        mg=MolGraph()
        with open(filename) as ff:
            for line in ff:
                line=line.split()
                if len(line)==0:
                    continue #skip if line[0] does not exist
                if "$atom:" in line[0]:
                    at=Atom()
                    at.label=line[0][6:]
                    at.element=at.label.strip("1234567890")
                    at.position=np.array([float(line[-3]),float(line[-2]),float(line[-1])])
                    mg.atom_list.append(at)
                    mg.label_to_ind_dic[at.label]=len(mg.atom_list)-1 #reconstruct later (after sort)
                    mg.atom_order.append(len(mg.atom_list)-1)
                elif "$bond:" in line[0]:
                    label_1=line[1][6:]
                    label_2=line[2][6:]
                    bonds.append([label_1,label_2]) # e.g. ["C01","C02"]
                    mg.atom_list[mg.label_to_ind_dic[label_1]].degree+=1
                    mg.atom_list[mg.label_to_ind_dic[label_2]].degree+=1
            n=len(mg.atom_list)
            mg.graph=[[0 for i in range(n)] for j in range(n)]
            for b in bonds:
                atom_index_1=mg.label_to_ind_dic[b[0]]
                atom_index_2=mg.label_to_ind_dic[b[1]]
                mg.graph[atom_index_1][atom_index_2]=1
                mg.graph[atom_index_2][atom_index_1]=1
            mindeg=n
            for a in mg.atom_list:
                mindeg=min(mindeg,a.degree)
            for i in range(n):
                if mg.atom_list[i].degree==mindeg:
                    mg.atom_list[i],mg.atom_list[0]=mg.atom_list[0],mg.atom_list[i]
                    mg.atom_order[i],mg.atom_order[0]=mg.atom_order[0],mg.atom_order[i]
                    for k in range(n):
                        mg.graph[k][i],mg.graph[k][0]=mg.graph[k][0],mg.graph[k][i]
                    mg.graph[i],mg.graph[0]=mg.graph[0],mg.graph[i]
                    break
            for i in range(1,n-1):
                for j in range(i,n):
                    if sum(mg.graph[j][:i])>0:
                        mg.atom_list[i],mg.atom_list[j]=mg.atom_list[j],mg.atom_list[i]
                        mg.atom_order[i],mg.atom_order[j]=mg.atom_order[j],mg.atom_order[i]
                        for k in range(n):
                            mg.graph[k][i],mg.graph[k][j]=mg.graph[k][j],mg.graph[k][i]
                        mg.graph[i],mg.graph[j]=mg.graph[j],mg.graph[i]
                        break
            for i in range(n):
                mg.label_to_ind_dic[mg.atom_list[i].label]=i
            lt_molgraph_list.append(mg)
            temp_bond_count=dict()
            for bond in bonds:
                elem1=bond[0].strip("1234567890")
                elem2=bond[1].strip("1234567890")
                if elem1<=elem2:
                    bond_type=elem1+","+elem2
                else:
                    bond_type=elem2+","+elem1
                if bond_type not in temp_bond_count.keys():
                    temp_bond_count[bond_type]=0
                temp_bond_count[bond_type]+=1
            for key in temp_bond_count.keys():
                if key not in bond_count.keys():
                    bond_count[key]=0
                bond_count[key]+=temp_bond_count[key]*num_molecule


poscar=Poscar.from_file(poscar_file).structure
for pos in poscar:
    pos.to_unit_cell(in_place=True)

mesh_length=2
xmin=0; ymin=0; zmin=0
xmax=0; ymax=0; zmax=0
for i in (0,1):
    for j in (0,1):
        for k in (0,1):
            x,y,z=poscar.lattice.get_cartesian_coords([i,j,k])
            xmin=min(xmin,x)
            xmax=max(xmax,x)
            ymin=min(ymin,y)
            ymax=max(ymax,y)
            zmin=min(zmin,z)
            zmax=max(zmax,z)

nx=int((xmax-xmin)/mesh_length)
ny=int((ymax-ymin)/mesh_length)
nz=int((zmax-zmin)/mesh_length)
grid=[[[[] for k in range(nz)] for j in range(ny)] for i in range(nx)]
for ind in range(len(poscar)):
    pos=poscar[ind]
    xind=int((pos.x-xmin)/(xmax-xmin)*nx)
    yind=int((pos.y-ymin)/(ymax-ymin)*ny)
    zind=int((pos.z-zmin)/(zmax-zmin)*nz)
    if xind==nx: xind=0
    if yind==ny: yind=0
    if zind==nz: zind=0
    grid[xind][yind][zind].append(ind)

bondset=set()
for ind in range(len(poscar)):
    pos=poscar[ind]
    xind=int((pos.x-xmin)/(xmax-xmin)*nx)
    yind=int((pos.y-ymin)/(ymax-ymin)*ny)
    zind=int((pos.z-zmin)/(zmax-zmin)*nz)
    if xind==nx: xind=0
    if yind==ny: yind=0
    if zind==nz: zind=0
    for i in (-1,0,1):
        for j in (-1,0,1):
            for k in (-1,0,1):
                neigh_xind=xind+i
                neigh_yind=yind+j
                neigh_zind=zind+k
                if neigh_xind==nx: neigh_xind=0
                if neigh_yind==ny: neigh_yind=0
                if neigh_zind==nz: neigh_zind=0
                for ind2 in grid[neigh_xind][neigh_yind][neigh_zind]:
                    if ind<ind2:
                        bondset.add((ind,ind2))
                    elif ind>ind2:
                        bondset.add((ind2,ind))

dists_cand=[]
for bond in bondset:
    i, j = bond
    dists_cand.append([poscar[i].distance(poscar[j]),bond])
dists_cand.sort()
dists=[]
for item in dists_cand:
    i, j = item[1]
    elem1=poscar[i].species.elements[0].symbol
    elem2=poscar[j].species.elements[0].symbol
    if elem1<=elem2:
        bond_type=elem1+","+elem2
    else:
        bond_type=elem2+","+elem1
    if bond_type in bond_count.keys() and bond_count[bond_type]>0:
        bond_count[bond_type]-=1
        dists.append(item)

#structure.lattice.alpha
retval=[]

atoms=[]
for i in range(len(poscar)):
    pos=poscar[i]
    at=Atom()
    at.position=pos.coords
    at.element=pos.species.elements[0].symbol
    at.poscar_index=i #set new attribute
    at.edge_to=[] #set new attribute
    atoms.append(at)
for item in dists:
    i,j=item[1]
    atoms[i].edge_to.append(j)
    atoms[j].edge_to.append(i)
for at in atoms:
    at.degree=len(at.edge_to)

uf=UnionFind(len(poscar))
for item in dists:
    i, j = item[1]
    uf.unite(i,j)
root_index_set=set()
for i in range(len(poscar)):
    root_index_set.add(uf.root(i))
for root_index in root_index_set:
    mg=MolGraph()
    poscar_index_to_atomlist_index=dict()
    for i in range(len(poscar)):
        if uf.root(i)==root_index:
            mg.atom_list.append(atoms[i])
    mg.atom_list.sort(key=attrgetter("degree"))
    for i in range(len(mg.atom_list)):
        at=mg.atom_list[i]
        poscar_index_to_atomlist_index[at.poscar_index]=i
    mg.graph=[[0 for i in range(len(mg.atom_list))] for j in range(len(mg.atom_list))]
    for i in range(len(mg.atom_list)):
        at=mg.atom_list[i]
        for to in at.edge_to:
            to=poscar_index_to_atomlist_index[to]
            mg.graph[i][to]=1
            mg.graph[to][i]=1
    retval.append(mg)

vasp_molgraph_list=retval
for vmg in vasp_molgraph_list:
    #print("before")
    #vmg.dump()
    #print()
    isok=False
    for lmg in lt_molgraph_list:
        isok=isomorph(vmg,lmg)
        if isok:
            break
    assert isok, "The VASP molecule is not isomorphic to any of the molecules given in .lt files"
    #print("after")
    #vmg.dump()
    #print()

with open(out_file,"w") as fo:
    fo.write("HEADER \n")
    fo.write("TITLE     Built with vasp2pdb                                            \n")
    fo.write("CRYST1")
    fo.write(f"{poscar.lattice.a:9.3f}")
    fo.write(f"{poscar.lattice.b:9.3f}")
    fo.write(f"{poscar.lattice.c:9.3f}")
    fo.write(f"{poscar.lattice.alpha:7.2f}")
    fo.write(f"{poscar.lattice.beta:7.2f}")
    fo.write(f"{poscar.lattice.gamma:7.2f}")
    fo.write(" ")
    fo.write("P 1        ")
    fo.write("   1")
    fo.write("\n")
    atom_index=1
    for mol_index in range(len(vasp_molgraph_list)):
        mg=vasp_molgraph_list[mol_index]
        for i in range(len(mg.atom_list)):
            a=mg.atom_list[mg.atom_order.index(i)]
            x,y,z=a.position
            fo.write("ATOM  ") # 1-6 Begin ATOM record
            fo.write(f"{atom_index: >5}") # 7-11 Atom serial number
            fo.write(" ") #12 Spacer
            fo.write(f"{a.label: <4}") #13-16 #Atom name
            fo.write(" ") #17 Alternate location indicator
            fo.write("MOL") #18-20 Residue name
            fo.write(" ") # 21 Spacer
            fo.write("A") # 22 Chain identifier
            fo.write(f"{mol_index+1: >4}") #23-26 Residue sequence number
            fo.write(" ") #Code for insertion of residues
            fo.write("   ")
            fo.write(f"{x:8.3f}{y:8.3f}{z:8.3f}") #31-54 Orthogonal x,y,z coordinates in Angstroms
            fo.write("  1.00") #Occupancy
            fo.write("  0.00") #Temperature factor
            fo.write("          ")
            fo.write(f" {a.element}") #Element symbol, right-justified
            #fo.write("  ") #Charge on the atom (Can be omitted?)
            fo.write("\n") # End of the record
            atom_index+=1
    fo.write("END\n")
