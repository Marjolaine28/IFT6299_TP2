# NGS error correction by k-mers frequency
## seq_correction.py :
Performs a one-base per read correction on a list of fastq files using a k-mers spectrum.

Exemple : seq_corrections --k=[12,13] --files=[R1.fastq,R2.fastq] returns 4 corrected fastq files (12-mers_correction_R1.fastq, 12-mers_correction_R2.fastq, 13-mers_correction_R1.fastq, 13-mers_correction_R2.fastq) 2 files reporting k-mers spectrum values (spectrum_12-mers.npy, spectrum_13-mers.npy), 2 files reporting occurence values per distinct k-mer (12-mers_counts.npy, 13-mers_counts.npy) and 2 files reporting the first extrema of the k-mers spectra (min-max_12-mers.npy, min-max_13-mers)

## spectrum_plotting.py :
Plots a spectrum of k-mers before and after the correction.

Exemple : spectrum_plotting.py --k=13 returns a png plot of the spectrum (13_mers_spectrum.png).

## stat_plotting.py :
Plots the statitics of correction.

Exemple : stat_corrections.py --k=[12,13] returns 1 png plot comparing the number of distinct k-mers before/after correction (stats_nk.png) and the number of correction (stats_nc.png), according to the k values given as parameter.
