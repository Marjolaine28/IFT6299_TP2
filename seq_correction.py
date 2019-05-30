from collections import Counter
from itertools import repeat
from scipy.signal import argrelmin
import numpy as np
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--k", type=list, default=[12,13,14,15])
parser.add_argument("--files", type=list, default=["trimmed_chry5_S14_L001_R1_001_paired.fastq","trimmed_chry5_S14_L001_R2_001_paired.fastq","trimmed_chry5_S14_L001_R1_001_unpaired.fastq"])

args = parser.parse_args()

enc={'A':0, 'C':1, 'G':2, 'T':3}

def convert_to_bits_rcomp(seq,k,b=0):
	for nt in seq:
		b = (b >> 2 | ((3-enc[nt]) << (2*k-2)))
	return b

def convert_to_bits(seq,k,b=0):
	for nt in seq:
		b = (b << 2 | enc[nt]) & (4**k -1)
	return b

def spectrum(k,corrected = False):
	occ= np.zeros(4**k, dtype='int32')
	if corrected :
		files = [str(k)+"-mers_correction_"+file for file in args.files]
	else :
		files = args.files
	for fastq in files :
		print(str(k)+" : counting "+fastq)
		with open(fastq,'r') as reads:
			for i,read in enumerate(reads):
				if (i+1)%2==0 and (i+1)%4!=0 and 'N' not in read:
					kmer = convert_to_bits(seq=read[0:k],k=k)
					rcomp = convert_to_bits_rcomp(seq=read[0:k],k=k)
					canonical=min(kmer,rcomp,2)
					occ[canonical]+=1
					for nt in range(k,len(read)-1):
						kmer = convert_to_bits(seq=read[nt],k=k,b=kmer)
						rcomp = convert_to_bits(seq=read[nt],k=k,b=rcomp)
						canonical=min(kmer,rcomp)
						occ[canonical]+=1
	if not corrected :
		np.save('./counts_'+str(k)+'-mers.npy', occ)
	occ = list(filter(lambda o: o != 0, occ))
	counter = Counter(occ)
	if corrected :
		np.save('./spectra_corrected_'+str(k)+'-mers.npy',list(counter.items()))
	else :
		np.save('./spectra_'+str(k)+'-mers.npy',list(counter.items()))
	return len(occ)

def first_min_max(spectrum_file,k):
	spectrum = np.load(spectrum_file)
	spectrum = sorted(spectrum, key=lambda kmer: kmer[0]) #kmer[0] corresponds to the occurences
	i_min = argrelmin(np.array(spectrum))[0][0]
	min = spectrum[i_min]
	i_max = np.argmax(spectrum[i_min:],axis=0)[1]+i_min
	max = spectrum[i_max]
	np.save('./min-max_'+str(k)+'-mers.npy',[min,max])
	return min


def correct_errors(spectrum_file,counts_file,k):
	min_occ = first_min_max(spectrum_file,k)[0]
	occ=np.load(counts_file)
	c=0
	for fastq in args.files :
		print(str(k)+" : correction of "+fastq)
		with open(fastq,'r') as reads:
			for i,read in enumerate(reads):
				corrections={}
				if i%4 == 0:
					id=read
				elif (i+1)%2 == 0 and (i+1)%4!=0 and 'N' not in read:
					seq=read[:-1]
				elif (i+1)%4 == 0:
					qual = read[:-1]
					m = min(qual)
					if m < '0':
						q_mins = [i for i, j in enumerate(qual) if j == m]
						for q_min in q_mins:
							corr_nt = list("ATCG")
							corr_nt.remove(seq[q_min])
							neighbors=repeat(list(seq),3)
							valid = 0
							for start,stop in zip(list(range(k))[::-1],range(k)):
								kmer = seq[q_min-start:q_min+stop+1]
								rcomp = convert_to_bits(seq=kmer,k=k)
								kmer = convert_to_bits(seq=kmer,k=k)
								canonical=min(kmer,rcomp)
								valid += int(occ[canonical]>min_occ)
							for neighbor,nt in zip(neighbors,corr_nt):
								neighbor[q_min] = nt
								valid_c = 0
								for start,stop in zip(list(range(k))[::-1],range(k)):
									kmer = neighbor[q_min-start:q_min+stop+1]
									rcomp = convert_to_bits(seq=kmer,k=k)
									kmer = convert_to_bits(seq=kmer,k=k)
									canonical=min(kmer,rcomp)
									valid_c += int(occ[canonical]>min_occ)				
								corrections[(q_min,nt)] = valid_c - valid
						vals=list(corrections.values())
						if len(vals)!=0 and max(vals) > 0 :
							c+=1
							keys=list(corrections.keys())			
							pos, nt = keys[np.argmax(vals)]
							read_c = list(seq)
							read_c[pos] = nt
							seq=''.join(read_c)
					with open(str(k)+"-mers_correction_"+fastq,'a') as reads_c:
						reads_c.write(id)
						reads_c.write(seq+'\n')
						reads_c.write("+"+'\n')
						reads_c.write(qual+'\n')
	return c

if __name__ == "__main__":
	c = []
	for k in args.k:
		spectrum(k)
		c.append(correct_errors("spectra_"+str(k)+"-mers.npy","counts_"+str(k)+"-mers.npy",k))
		spectrum(k,corrected=True)
	np.save("stat_corrections.npy",c)

