USAGE:
	perl Calculate_PC_Over_Time.pl Arg1 Arg2 Arg3 Arg4 Arg5

This script takes 5 input arguments (that follow the script name, and are separated by spaces):
1.)	Name of the FASTA file, which must also include the full file path if the script is not in the 
same directory as the file (e.g. /Users/Brian/SomeDir/File.fasta).

2.)	Name of the table that contains the collection dates of each isolate. This must be tab-delimited 
with the sequence name in the first column (same as in FASTA file) and the collection date in the second column.

3.)	The minor allele frequency (MAF) cutoff to use. Polymorphisms with frequencies below this number 
(or above 1-MAF) are excluded. This number should be between 0 and 0.5. I recommend nothing lower than 0.05
or 0.1, since these rare polymorphisms are not informative for recombination (likely too young to have
been involved in any transfers). You could run the script with several different parameter values to
see how results change.

4.)	The number of windows you wish to break the alignment into. Specify “1” to look at the whole 
alignment, “2” to analyze the alignment in two halves, etc. This number should be less than or equal 
to half the input sequence length, or mayhem will ensue.

5.)	Whether or not the analysis is done forward-in-time (F) or reverse-in-time (R). This determines 
whether 3 sequences with the oldest collection dates (F) or with the most recent collection dates (R) 
are used as the “anchor”. Then, the script cycles through all sequences not included in this trio and 
asks, for each pair of SNPs, whether the new fourth sequence contains a fourth haplotype, which is 
potential evidence of recombination (“potential” because it is also explained by recurrent mutation, 
but not likely). Suggestion: choose “F” or “R” to make the three “anchor” sequences have the same 
collection date, although this is not necessary. For instance, if you have at least 3 sequences from 
the most recent collection date but only 1 from the oldest date, choose “R”, or if you have at least 
3 sequences from the oldest collection date but only 1 from the most recent date, choose “F”. The choice 
really doesn’t matter much, but determines whether there is potentially temporal signal within the 
3 “anchor” sequences. If you use option "R", remember to analyze linkage backwards in time! That is,
linkage disequilibrium should decay as a function of time from the collection dates of the 3 "anchor"
sequences, which is backwards in time if option "R" is selected.

Input file in FASTA format, with sequence names following the “>” character, and the following line(s)
contains the sequence data, either all on the next line for that particular sequence or printed across
many lines, for example with 60 nucleotides per line. Either is fine. However, make sure the sequence
length for each individual is the same! These FASTA files may contain all sites OR only polymorphic
sites. No tabs or whitespaces after the sequence names. As it’s currently written, the code will
tolerate biallelic sites of any kind, and may be nucleotides (A,G,C,T) or 0’s and 1’s such as the
output from some simulators. However, ambiguous nucleotides (such as M, R, W, S, Y etc…) should
not be used. The code ignores all sites in which at least one individual has an “N” or “-“, and
all multiallelic sites in which more than 2 bases are present, although these are typically a 
small fraction of polymorphic sites in most species.

Make sure your file doesn’t have any hidden tabs or characters that may get introduced by other programs.  
I attempted to have the script remove carriage returns introduced by DOS/Windows, but not other hidden 
characters.

Please bear in mind that this script may take a bit of time (a few hours or more than a day, depending 
on dataset), because not only does it go through all pairs of polymorphisms, but does so n-3 times, where 
n in the sample size. If the script is taking a long time, you could sample a fraction of the sequences 
or a fraction of the sites. For instance, if you have 100 individuals all from the exact same collection 
date, I recommend downsampling individuals from this time point. You could also increase the number of 
windows, which I had to do for Helicobacter pylori (it has an insane number of polymorphisms). Lastly, 
you could increase the MAF cutoff (e.g. to 0.2), making the script go through fewer polymorphisms and 
only explore higher frequency polymorphisms.

The output table has four columns:
1.)	Window number (starts at 0, not at 1)
2.)	The name of the fourth sequence that was compared the “anchor” trio
3.)	The sample date for the fourth sequence
4.)	The fraction of polymorphism (SNP) pairs that exhibited fewer than 4 haplotypes when this fourth 
sequence was compared to the “anchor” trio. This is the fraction of SNP pairs that exhibit no evidence
 of recombination and are “pairwise compatible”, or compatible with a single phylogeny.

The output IS NOT a correlation analysis, but may be used in other programs such as R or Excel to look
at the relationship between the collection date of the fourth individual (3rd column) and pairwise 
compatibility (4th column).

This script could also potentially be used to look at the evolution of recombination hotspots through
time. For instance, if a set of contemporary sequences look hyper-recombinant in a particular region,
you could see if this signal is also present through time by using this script to look at this region,
either by making an alignment of this region by itself and feeding it to the program, or by choosing
a window size such that one of the windows covers this region. You would obviously need to do this
analysis on other regions that are not suspected of recombining more than average as a sort of
“control”.

If you run into problems of any kind, please feel free to contact me at brianjohnarnold@gmail.com. Thanks!!


