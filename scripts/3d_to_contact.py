'''
3d_to_contact.py
################
This is a Python 3 script to load .pdb files and convert their coordinates to 
an RR format contact map file.
'''
###############################################################################
# LIBRARIES
###############################################################################
import numpy as np

###############################################################################
# FUNCTIONS
###############################################################################
def load_fasta_sequence(file):
	fp = open(file,'r')
	s = list(fp)
	fp.close()
	return s[1].rstrip()

def load_pdb_file(file):
	c,p,x,y,z = [],[],[],[],[]
	fp = open(file,'r')
	for line in fp:
		if line[0:4] == 'ATOM':
			if line[17:20].strip()=='GLY':
				if line[12:16].strip()=='CA':
					c.append(line[21])
					p.append(line[22:26].replace(' ',''))
					x.append(float(line[30:38]))
					y.append(float(line[38:46]))
					z.append(float(line[46:54]))
			else:
				if line[12:16].strip()=='CB':
					c.append(line[21])
					p.append(line[22:26].replace(' ',''))
					x.append(float(line[30:38]))
					y.append(float(line[38:46]))
					z.append(float(line[46:54]))
	fp.close()
	if len(set(c))>1:
		indices = np.core.defchararray.add(c,p)
		result  = np.column_stack((indices,x,y,z))
	else:
		result = np.column_stack((p,x,y,z))
	return result

def create_RR(coords):
	temp = []
	for i in range(coords.shape[0]):
		d = np.linalg.norm(coords[i,1:].astype(np.float)-coords[:i+1,1:].astype(np.float),axis=1)
		for j in range(d.shape[0]):
			if d[j] < 8.0 and d[j] > 0.0:
				if (i-j) > 6:
					temp.append([coords[j,0],coords[i,0],0,8,1.0])
	RR = np.array(temp)
	return RR		

def write_RR_file(file,target,RR,sequence):
	fp = open(file,'w')
	fp.write('PFRMAT RR\n')
	fp.write('TARGET '+target+'\n')
	# fp.write('AUTHOR 1234-5678-9000\n')
	# fp.write('REMARK None\n')
	fp.write('METHOD Own Python mplementation\n')
	fp.write('MODEL 1\n')
	splitted_sequence = [sequence[i:i+50] for i in range(0,len(sequence),50)]
	for s in splitted_sequence:
		fp.write(s+'\n')
	for r in RR:
		fp.write(r[0]+' '+r[1]+' '+r[2]+' '+r[3]+' '+r[4]+'\n')
	fp.close()
	return

###############################################################################
# MAIN
###############################################################################
def main():
	seq_easy = load_fasta_sequence('../fasta/T0951_easy/T0951.fasta')
	seq_medi = load_fasta_sequence('../fasta/T0866/T0866.fasta')
	print("SEQUENCE LENGTHS")
	print("Easy:   ",len(seq_easy))
	print("Medium: ",len(seq_medi))

	pdb_easy = load_pdb_file('../pdb/native/T0951_easy/5z82.pdb')
	pdb_medi = load_pdb_file('../pdb/native/T0866/5uw2.pdb')
	print("PDB BACKBONES (C-beta)")
	print("PDB Easy:   ",pdb_easy.shape)
	print("PDB Medium: ",pdb_medi.shape)

	RR_easy = create_RR(pdb_easy)
	RR_medi = create_RR(pdb_medi)
	print("RR FILE FOR <8A AND 6 RESIDUES APART")
	print("RR Easy:   ",RR_easy.shape)
	print("RR Medium: ",RR_medi.shape)

	write_RR_file('../contact_maps/own/true_T0951_easy.RR','T0951',RR_easy,seq_easy)
	write_RR_file('../contact_maps/own/true_T0866.RR','T0966',RR_medi,seq_medi)

if __name__ == '__main__':
	main()