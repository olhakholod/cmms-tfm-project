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
					x.append(line[30:38].strip())
					y.append(line[38:46].strip())
					z.append(line[46:54].strip())
			else:
				if line[12:16].strip()=='CB':
					c.append(line[21])
					p.append(line[22:26].replace(' ',''))
					x.append(line[30:38].strip())
					y.append(line[38:46].strip())
					z.append(line[46:54].strip())
	fp.close()
	# result = np.full((len(sequence),4),"--------")
	if len(set(c))>1:
		results = []
		temp1 = np.column_stack((p,x,y,z))
		for chain in sorted(set(c)):
			temp2 = temp1[np.where(np.array(c)==chain)]
			# for i in temp2:
			# 	result[int(i[0])-1] = i
			# results.append(result)
			# result = np.full((len(sequence),4),"--------")
			results.append(temp2)
		return np.array(results) # return one array per chain.
	else:
		# for i in range(len(p)):
			# result[int(p[i])-1] = (p[i],x[i],y[i],z[i])
		result = np.column_stack((p,x,y,z))
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
	return c

def contact_matrix(coords,sequence):
	s = len(sequence)
	'''
	Create a distance matrix for a native structure's
	coordinates (previously retrieved from its PDB file).
	'''
	d = np.ndarray((s,s))
	# d = np.ndarray((coords.shape[0],coords.shape[0]))
	for i in range(coords.shape[0]):
		x = np.linalg.norm(coords[:i+1,1:].astype(np.float)-coords[i,1:].astype(np.float),axis=1)
		p = coords[:i+1,0].astype(np.int)
		# d[i,:i+1] = x
		d[p[i],p[:i+1]] = x
	dt = d.T
	d = dt + d
	for j in range(d.shape[0]):
		d[j,j] = 9999.999
		z = d[j,j+1:j+6]
		z[:] = 9999.999
		z = d[j+1:j+6,j]
		z[:] = 9999.999
	r1 = d > 0.0
	r2 = d < 8.0
	return (r1 & r2)

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
		if D[int(RR[i,0])-1,int(RR[i,1])-1]:
			TP += 1
		else:
			FP += 1
	T = np.sum(D==1)
	F = np.sum(D==0)
	N = T + F
	FN = T - TP
	TN = F - FP
	print("TP\t FP\t FN\t TN")
	print(TP,'\t',FP,'\t',FN,'\t',TN)
	print("--------------------------------------")
	print("PRECISION:         ",precision(TP,FP))
	print("RECALL/COVERAGE:   ",recall(TP,FN))
	print("F-SCORE:           ",f_score(precision(TP,FP),recall(TP,FN)))
	print("IMPROVEMENT OVER R:",impR(TP,FP,T,N))
	print("MCC:               ",mcc(TP,FP,FN,TN))

#####################################################################
# MAIN
#####################################################################
def main():
	seq1  = load_fasta_sequence('../fasta/T0951_easy/T0951.fasta')
	seq2  = load_fasta_sequence('../fasta/T0866/T0866.fasta')

	pdb1 = load_pdb_file('../pdb/native/T0951_easy/5z82.pdb',seq1)
	pdb2 = load_pdb_file('../pdb/native/T0866/5uw2.pdb',seq2)
	pdb2_a = pdb2[0]
	pdb2_b = pdb2[1]
	pdb2_c = pdb2[2]

	D1  = contact_matrix(pdb1,seq1)
	D2A = contact_matrix(pdb2_a,seq2)
	D2B = contact_matrix(pdb2_b,seq2)
	D2C = contact_matrix(pdb2_c,seq2)

	sseq1 = splitted_sequence(seq1)
	sseq2 = splitted_sequence(seq2)

	RR1_CMA = load_RR('../contact_maps/predicted/CMApro/T0951.rr',sseq1)
	RR2_CMA = load_RR('../contact_maps/predicted/CMApro/T0866.rr',sseq2)
	RR1_Raptor = load_RR('../contact_maps/predicted/RaptorX/T0951.rr',sseq1)
	RR2_Raptor = load_RR('../contact_maps/predicted/RaptorX/T0866.rr',sseq2)
	RR1_SVMcon = load_RR('../contact_maps/predicted/SVMcon/T0951.rr',sseq1)
	RR2_SVMcon = load_RR('../contact_maps/predicted/SVMcon/T0866.rr',sseq2)


	print("======================================")
	print(" Results for CMAPro -- TARGET: T0951")
	print("======================================")
	evaluate_RR_file(RR1_CMA,D1)

	print("======================================")
	print(" Results for CMAPro -- TARGET: T0866")
	print("======================================")
	evaluate_RR_file(RR2_CMA,D2A)

if __name__ == '__main__':
	main()