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
import matplotlib.pyplot as plt
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
				if l[0:6] not in {'PFRMAT','TARGET','AUTHOR','METHOD','MODEL ','REMARK','FRMAT ','SEQRES'}:
					if l not in set(splitted_sequence):
						s = l.split()
						if s[0] == 'CONTC':
							temp.append([s[1],s[2],s[-1]])
						else:
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
	# d = np.triu(d,6)
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
	Same as 'accuracy' in older CASPs.
	'''
	return np.round(TP/(TP+FP),3)

def recall(TP,FN):
	'''
	Same as coverage since TP+FN = D True
	'''
	return np.round(TP/(TP+FN),3)

def f_score(p,r):
	if (p+r)==0:
		return 0
	return np.round(2*p*r/(p+r),3)

def impR(TP,FP,T,N):
	'''
	Improvement over random.
	'''
	p = TP/(TP+FP)
	return np.round(p/(T/N),3)

def mcc(TP,FP,FN,TN):
	'''
	Matthew's Correlation Coefficient(MCC).
	'''
	num = TP*TN - FP*FN
	den = np.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
	return np.round(num/den,3)

def eval_subset(RR,D,T,F):
	sub_tp,sub_fp = 0,0
	for i in range(RR.shape[0]):
		if D[int(RR[i,0])-1,int(RR[i,1])-1]:
			sub_tp += 1
		else:
			sub_fp += 1
	p = precision(sub_tp,sub_fp)
	r = np.round(100*recall(sub_tp,T - sub_tp),3)
	f = f_score(p,r)
	i = impR(sub_tp,sub_fp,T,F+T)
	return p,r,f,i

def evaluate_RR_file(RR,D,coords,seq):
	# Calculate mismatches between sequence and PDB. Skip missing residues in PDB
	new_seq = list('-'*len(seq))
	for i in coords:
		new_seq[int(i[0])] = seq[int(i[0])]
	missing = [i for i in range(len(new_seq)) if new_seq[i]=='-']

	# Calculate TP, FP between RR file and distance matrix
	d_missing = 0
	rr_not_missing = []
	TP,FP = 0,0
	for i in range(RR.shape[0]):
		if int(RR[i,0]) not in missing and int(RR[i,1]) not in missing:
			if D[int(RR[i,0])-1,int(RR[i,1])-1]:
				TP += 1
			else:
				FP += 1
			rr_not_missing.append(i)
		else:
			d_missing +=1
	
	T = np.sum(np.triu(D,1))
	F = np.sum(np.triu(np.invert(D),1))	- d_missing
	FN = T - TP
	TN = F - FP

	# Summary
	print("Predicted contacts in RR: ", RR.shape[0])
	print("Residues not in native:   ",len(missing))
	print("Contacts not in native:   ",d_missing)
	print("Remaining: ",RR.shape[0]-d_missing)
	print("MCC: ",mcc(TP,FP,FN,TN))
	print("--------------------------------------------------------------------")
	print("TP\t FP\t FN\t TN")
	print(TP,'\t',FP,'\t',FN,'\t',TN)
	print(T,'\t' '\t',F)
	print("--------------------------------------------------------------------")
	print("                 TOP-5\tL/10\tL/5\tL/2\tL\t2L")

	new_RR = RR[rr_not_missing]
	
	#SORTED BY PROBABILITIES
	sortedRR = np.argsort(new_RR[:,2].astype(np.float))
	top_5  = new_RR[sortedRR[-5:]]
	l_10   = new_RR[sortedRR[-round(len(seq)/10):]]
	l_5    = new_RR[sortedRR[-round(len(seq)/5):]]
	l_2    = new_RR[sortedRR[-round(len(seq)/2):]]
	top_l  = new_RR[sortedRR[-len(seq):]]
	top_2l = new_RR[sortedRR[-len(seq)*2:]]

	# RANDOM
	# top_5 = new_RR[np.random.choice(rr_not_missing,5)]
	# l_10  = new_RR[np.random.choice(rr_not_missing,round(len(seq)/10))]

	# RESULTS
	p_top_5, r_top_5,f_top_5,i_top_5 = eval_subset(top_5,D,T,F)
	p_l_10, r_l_10, f_l_10, i_l_10   = eval_subset(l_10,D,T,F)
	p_l_5, r_l_5, f_l_5, i_l_5 = eval_subset(l_5,D,T,F)
	p_l_2, r_l_2, f_l_2, i_l_2 = eval_subset(l_2,D,T,F)
	p_top_l, r_top_l, f_top_l, i_top_l = eval_subset(top_l,D,T,F)
	p_top_2l, r_top_2l, f_top_2l, i_top_2l = eval_subset(top_2l,D,T,F)

	print("PRECISION:       %.3f\t%3.3f\t%3.3f\t%3.3f\t%3.3f\t%3.3f" % (p_top_5,p_l_10,p_l_5,p_l_2,p_top_l,p_top_2l))
	print("RECALL/COVERAGE: %.3f\t%3.3f\t%3.3f\t%3.3f\t%3.3f\t%3.3f" % (r_top_5,r_l_10,r_l_5,r_l_2,r_top_l,r_top_2l))
	print("F-SCORE:         %.3f\t%3.3f\t%3.3f\t%3.3f\t%3.3f\t%3.3f" % (f_top_5,f_l_10,f_l_5,f_l_2,f_top_l,f_top_2l))
	print("ImpR             %.3f\t%3.3f\t%3.3f\t%3.3f\t%3.3f\t%3.3f" % (i_top_5,i_l_10,i_l_5,i_l_2,i_top_l,i_top_2l))

	p = [p_top_5,p_l_10,p_l_5,p_l_2,p_top_l,p_top_2l]
	r = [r_top_5,r_l_10,r_l_5,r_l_2,r_top_l,r_top_2l]
	f = [f_top_5,f_l_10,f_l_5,f_l_2,f_top_l,f_top_2l]
	i = [i_top_5,i_l_10,i_l_5,i_l_2,i_top_l,i_top_2l]
	return p,r,f,i

# def plots():
# 	subsets   = ['Top-5','L/10','L/5','L/2','L','2L']
# 	models    = ['CMApro','RaptorX','SVMcon','DCON2']
# 	fig,ax = plt.subplots()
# 	lines = ax.plot(names,precision)
# 	ax.legend()
# 	plt.savefig('../contact_maps/method_comparison.png')

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
	RR1_DNCON2 = load_RR('../contact_maps/predicted/DNCON2Contacts/T0951/T0951.dncon2.rr',sseq1)
	RR1_coneva = load_RR('../contact_maps/true/coneva/T0951.RR',sseq1)
	RR1_TRUE   = load_RR('../contact_maps/true/own/T0951_easy.rr',sseq1)

	precision,recall,fscores,impr = [],[],[],[]

	print("====================================================================")
	print(" Results for CMAPro -- TARGET: T0951")
	print("====================================================================")
	p,r,f,i = evaluate_RR_file(RR1_CMA,D1,pdb1,seq1)
	precision.append(p)
	recall.append(r)
	fscores.append(f)
	impr.append(i)

	print("====================================================================")
	print(" Results for RaptorX -- TARGET: T0951")
	print("====================================================================")
	p,r,f,i = evaluate_RR_file(RR1_Raptor,D1,pdb1,seq1)
	precision.append(p)
	recall.append(r)
	fscores.append(f)
	impr.append(i)

	print("====================================================================")
	print(" Results for SVMcon -- TARGET: T0951")
	print("====================================================================")
	p,r,f,i = evaluate_RR_file(RR1_SVMcon,D1,pdb1,seq1)
	precision.append(p)
	recall.append(r)
	fscores.append(f)
	impr.append(i)

	print("====================================================================")
	print(" Results for DNCON2 -- TARGET: T0951")
	print("====================================================================")
	p,r,f,i = evaluate_RR_file(RR1_DNCON2,D1,pdb1,seq1)
	precision.append(p)
	recall.append(r)
	fscores.append(f)
	impr.append(i)

	print("====================================================================")
	print(" Results for ConEVA -- TARGET: T0951")
	print("====================================================================")
	p,r,f,i = evaluate_RR_file(RR1_coneva,D1,pdb1,seq1)
	precision.append(p)
	recall.append(r)
	fscores.append(f)
	impr.append(i)

	print("====================================================================")
	print(" Results for TRUE -- TARGET: T0951")
	print("====================================================================")
	p,r,f,i = evaluate_RR_file(RR1_TRUE,D1,pdb1,seq1)
	precision.append(p)
	recall.append(r)
	fscores.append(f)
	impr.append(i)

	subsets   = ['Top-5','L/10','L/5','L/2','L','2L']
	models    = ['CMApro','RaptorX','SVMcon','DCON2','ConEVA','TRUE']
	for i in range(len(models)):
		plt.plot(subsets,precision[i],label=models[i])
	plt.legend()
	plt.savefig('../contact_maps/precision.png')

	plt.figure()
	for i in range(len(models)):
		plt.plot(subsets,recall[i],label=models[i])
	plt.legend()
	plt.savefig('../contact_maps/coverage.png')

	# print("====================================================================")
	# print(" Results for CMAPro -- TARGET: T0866")
	# print("====================================================================")
	# evaluate_RR_file(RR2_CMA,D2A,pdb2_a,seq2)

	# print("====================================================================")
	# print(" Results for RaptorX -- TARGET: T0951")
	# print("====================================================================")
	# evaluate_RR_file(RR2_Raptor,D2A,pdb2_a,seq2)

	# print("====================================================================")
	# print(" Results for SVMcon -- TARGET: T0951")
	# print("====================================================================")
	# evaluate_RR_file(RR2_SVMcon,D2A,pdb2_a,seq2)


if __name__ == '__main__':
	main()