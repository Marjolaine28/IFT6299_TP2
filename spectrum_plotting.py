import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--k", type=str, default="12")


args = parser.parse_args()


if __name__ == "__main__":
	sc= np.transpose(np.load("spectrum_corrected_"+args.k+"-mers.npy"))
	s= np.transpose(np.load("spectrum_"+args.k+"-mers.npy"))
	max=np.load("min-max_"+args.k+"-mers.npy")[1]
	plt.xlim(0,max[0]*3)
	plt.ylim(0,max[1]*1.5)
	plt.scatter(s[0],s[1],s=3,label="original")
	plt.scatter(sc[0],sc[1],s=3,c="orange",label="corrected")
	plt.title(args.k+" mers frequency spectra (Agrobact. T. Chry5)")
	plt.xlabel('Occurence')
	plt.ylabel('Frequency')
	plt.legend()
	plt.tight_layout()
	plt.savefig("spectrum_"+args.k+"-mers.png")

