from glob import glob
import pdbreader
#compound_name = input("Please write compound code")
compound_name = "17T_"
pdb_path = glob("Compounds/*")[0]
print(pdb_path)
pdb = pdbreader.read_pdb(pdb_path)

atoms = pdb["ATOM"]
n_atoms = len(atoms)

bonds = pdb["CONECT"]
pairbonds = []
anglebonds = []
torsionbonds =[]
#for bond in bonds:
    #print(bond)
print(bonds)
print(bonds[0])
exit()

n_bonds = len(bonds)
print(atoms)
print(bonds)
print("N atoms: ", n_atoms)
print("N bonds: ", n_bonds)
#print(pdb)



