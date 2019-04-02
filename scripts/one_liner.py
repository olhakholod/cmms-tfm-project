# one_liner.py
# Read fasta file sequence and save as a single line element

def load_fasta_sequence(file):
	fp = open(file,'r')
	s = list(fp)
	fp.close()
	return s[1].rstrip()

def write_one_liner(file,sequence):
	fp = open(file,'w')
	fp.write(sequence)
	fp.close()

def main():
	seq2 = load_fasta_sequence('../fasta/T0866/T0866.fasta')
	seq1 = load_fasta_sequence('../fasta/T0951_easy/T0951.fasta')
	write_one_liner('../fasta/T0866/single_line_T0866.fasta',seq2)
	write_one_liner('../fasta/T0951_easy/single_line_T0951.fasta',seq1)

if __name__ == '__main__':
	main()