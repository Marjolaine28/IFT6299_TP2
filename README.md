# NGS error correction by k-mers frequency
## seq_correction.py :
Performs a one-base per read correction on a list of fastq files using a k-mers spectrum.
exemple : seq_corrections --k=[12,13] --files=[R1.fastq,R2.fastq] returns 4 files : 12-mers_correction_R1.fastq, 12-mers_correction_R2.fastq, 13-mers_correction_R1.fastq, 13-mers_correction_R2.fastq.
## spectrum_plotting.py :
Plots a spectrum of k-mers before and after the correction
## stat_plotting.py :
Plots the statitics of correction (number of distinct k-mers before/after correction, number of correction)
