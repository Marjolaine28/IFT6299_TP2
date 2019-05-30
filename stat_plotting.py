import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--k", type=list, default=[12,13,14,15])

args = parser.parse_args()

def number_corrections():
	c = np.load("stat_corrections.npy")
	plt.title("Corrections statistics")
	plt.xlabel("k")
	plt.ylabel("Number of corrections")
	plt.xticks(args.k)
	plt.scatter(args.k,c)
	plt.tight_layout()
	plt.savefig("stats_nc.png")
	plt.clf()


def number_kmers():
	n = []
	nc=[]
	for k in args.k:
		s = np.transpose(np.load("spectra_"+str(k)+"-mers.npy"))
		sc = np.transpose(np.load("spectra_corrected_"+str(k)+"-mers.npy"))
		n.append(np.sum(s[1]))
		nc.append(np.sum(sc[1]))
	print(nc,n)
	ratio = [a/b for a,b in zip(nc,n)]
	plt.title("Comparaison of the number of distinct k-mers")	
	plt.xlabel("k")
	plt.ylabel("after correction / before correction (ratio)")
	plt.xticks(args.k)
	plt.scatter(args.k,ratio)
	plt.tight_layout()
	plt.savefig("stats_nk.png")

if __name__ == "__main__":
	number_corrections()
	number_kmers()
