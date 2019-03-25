'''
3d_to_contact.py
################
This is a Python 3 script to load .pdb files and convert their coordinates to a
 2d array of contact maps.
'''
###############################################################################
# LIBRARIES
###############################################################################
import numpy as np


###############################################################################
# FUNCTIONS
###############################################################################
def load_pdb_file(file):
	a,x,y,z,t = [],[],[],[],[]
	fp = open(file,'r')
	for line in fp:
		if line[0:4] == 'ATOM':
			if line[12:16].strip()=='CA': #Change to carbon-beta
				a.append(line[17:20])
				x.append(float(line[30:38]))
				y.append(float(line[38:46]))
				z.append(float(line[46:54]))
				t.append(float(line[60:66]))
	fp.close()
	new_a,new_x,new_y,new_z,new_t = [],[],[],[],[]
	j = 0
	for i in sequence:
		if i!='-':
			new_a.append(a[j])
			new_x.append(x[j])
			new_y.append(y[j])
			new_z.append(z[j])
			new_t.append(t[j])
			j += 1
		else:
			new_a.append('UNK')
			new_x.append(math.inf)
			new_y.append(math.inf)
			new_z.append(math.inf)
			new_t.append(0)
	c = np.array([[new_x[i],new_y[i],new_z[i]] for i in range(len(new_x))])
	return c

def create_RR(x,y,z):
	RR = np.zeros((len(x),5))
	return RR

###############################################################################
# MAIN
###############################################################################
def main():
	print("Hello World from 3d_to_contact.py!")

if __name__ == '__main__':
	main()