#####################################################################
#contact_map_eval.py
'''
Python3 script to compute several metrics across a set of given 
residue-residue contact maps.
'''
#####################################################################
# LIBRARIES
#####################################################################
import numpy as np

#####################################################################
# FUNCTIONS
#####################################################################
def load_fasta_sequence(file):
	fp = open(file,'r')
	s = list(fp)
	fp.close()
	return s[1].rstrip()

def splitted_sequence(sequence):
	s = [sequence[i:i+50] for i in range(0,len(sequence),50)]
	return s

def load_pdb_file(file,sequence):
	c,p,x,y,z = [],[],[],[],[]
	fp = open(file,'r')
	for line in fp:
		if line[0:4] == 'ATOM':
			if line[17:20].strip()=='GLY':
				if line[12:16].strip()=='CA':
					c.append(line[21])
					p.append(line[22:26].replace(' ',''))
					x.append(line[30:38])
					y.append(line[38:46])
					z.append(line[46:54])
			else:
				if line[12:16].strip()=='CB':
					c.append(line[21])
					p.append(line[22:26].replace(' ',''))
					x.append(line[30:38])
					y.append(line[38:46])
					z.append(line[46:54])
	fp.close()
	result = np.full((len(sequence),4),'999')
	if len(set(c))>1:
		results = []
		temp1 = np.column_stack((p,x,y,z))
		for chain in sorted(set(c)):
			temp2 = temp1[np.where(np.array(c)==chain)]
			for i in temp2:
				result[int(i[0])-1] = i
			results.append(result)
			result = np.full((len(sequence),4),'999')
		return np.array(results) # return one array per chain.
	else:
		for i in range(len(p)):
			result[int(p[i])-1] = (p[i],x[i],y[i],z[i]) 
		return result

def load_RR(file,splitted_sequence):
	'''
	Open an .rr file following the CASP RR format, and return its 
	listed contact pairs and probabilities.
	'''
	start = 0
	temp = []
	fp = open(file,'r')
	for line in fp:
		l = line.rstrip()
		if l!='':
			if l[0:3]!='END':
				if l[0:6] not in {'PFRMAT','TARGET','AUTHOR','METHOD','MODEL '}:
					if l not in set(splitted_sequence):
						s = l.split()
						temp.append([s[0],s[1],s[-1]])
	c = np.array(temp)
	print(c.shape)
	return c

def distance_matrix(coords):
	'''
	Create a distance matrix for a native structure's
	coordinates (previously retrieved from its PDB file).
	'''
	# d = np.ndarray((coords.shape[0],coords.shape[0]))
	d = np.ndarray((coords.shape[0],coords.shape[0]))
	for i in range(coords.shape[0]):
		x = np.linalg.norm(
			coords[:i+1,1:].astype(np.float)-
			coords[i,1:].astype(np.float),axis=1
			)
		d[i,:i+1] = x
	dt = d.T
	d = dt + d
	for j in range(0,d.shape[0]):
		d[j,j] = 999.0
		z = d[j,j+1:j+6]
		z[:] = 999.0
		z = d[j+1:j+6,j]
		z[:] = 999.0
	print(d.shape)
	return (d < 8.0)

def precision(TP,FP):
	'''
	Same as accuracy in older CASP.
	'''
	return TP/(TP+FP)

def recall(TP,FN):
	'''
	Same as coverage since TP+FN = D True
	'''
	return TP/(TP+FN)

def f_score(p,r):
	return 2*p*r/(p+r)

def impR(TP,FP,T,N):
	'''
	Improvement over random.
	'''
	p = TP/(TP+FP)
	return p/(T/N)

def mcc(TP,FP,FN,TN):
	'''
	Matthew's Correlation Coefficient(MCC).
	'''
	num = TP*TN - FP*FN
	den = np.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
	return num/den

def evaluate_RR_file(RR,D):
	TP,FP = 0,0
	in_RR = RR.shape[0]
	for i in range(RR.shape[0]):
		if D[int(RR[i,0]),int(RR[i,1])]:
			TP += 1
		else:
			FP += 1
	T = np.sum(D==1)
	F = np.sum(D==0)
	N = T + F
	FN = T - TP
	TN = F - FP
	print(TP,FP,FN,TN)

#####################################################################
# MAIN
#####################################################################
def main():
	seq1 = load_fasta_sequence('../fasta/T0951_easy/T0951.fasta')
	seq2 = load_fasta_sequence('../fasta/T0866/T0866.fasta')

	sseq1 = splitted_sequence(seq1)
	sseq2 = splitted_sequence(seq2)

	RR1 = load_RR('../contact_maps/predicted/CMApro/T0951.rr',sseq1)
	RR2 = load_RR('../contact_maps/predicted/CMApro/T0866.rr',sseq2)

	pdb1 = load_pdb_file('../pdb/native/T0951_easy/5z82.pdb',seq1)
	pdb2 = load_pdb_file('../pdb/native/T0866/5uw2.pdb',seq2)

	print(pdb1.shape)
	print(pdb2.shape)

	D1 = distance_matrix(pdb1)
	# D2 = distance_matrix(pdb2)

	# evaluate_RR_file(RR1,D1)

if __name__ == '__main__':
	main()