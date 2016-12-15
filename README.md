# BOND
Computation of perfect DNA oligonucleotides

==========================
# Compiling BOND components
==========================


BOND: The bond.cpp is the whole BOND program and to create the executable:

	g++ -O3 -fopenmp bond.cpp -o BOND


EVALUATION: pre_eval.cpp is the program that gets the output of bond (all detail i.e. melting temp., ...) and creates the oligos in plain format, i.e. only oligos in each line. To create the executable:
 
	g++ -O3 pre_eval.cpp -o PRE_EVAL

EVALUATION: eval.cpp is the evaluation program which is created as follows:

	g++ -O3 eval.cpp -o EVAL


MAXIMUM OLIGOS: The program to compute maximum number of oligos includes: main.cpp, genome.cpp, genome.h. To create the executable:

	g++ -O3 -fopenmp main.cpp genome.cpp -o MAX


============================================================
How to run BOND (our program for designing oligonucleotides)
============================================================

	./BOND_<platform> <inputSequences> <outputOligos> [-length] [-seqSim] [-maxMatch] [-maxGC] [-minGC] [-dimerLen] [-dimerSim] [-hairpinLen] [-minhpLoop] [-maxhpLoop] [-rangeTm] [-minTm] [-maxTm] [-oligCon] [-saltCon] [-secStr]  

Required parameters:

< inputSequencess > - the input file containing the gene sequences in FASTA format

< outputOligos > - the output file containing the designed oligos and additional information (gene number, oligo sequence, length, melting temperature, distance from 5' and 3' ends)


Optional parameters:

[-length]:     The length of oligo (Default value is 50)

[-seqSim]:     Sequence identity percentage (Default value is 75)

[-maxMatch]:   Maximum consecutive match (Default value is 15)

[-maxGC]:      Maximum GC-content percentage (Default value is 70)

[-minGC]:      Minimum GC-content percentage (Default value is 30)

[-dimerLen]:   The length of dimer window (Default value is 15)

[-dimerSim]:   Dimer stringency percentage (Default value is 85)

[-hairpinLen]: The length of hairpin stem (Default value is 6)

[-minhpLoop]:  Minimum length of hairpin loop (Default value is 1)

[-maxhpLoop]:  Maximum length of hairpin loop (Default value is 3)

[-minTm]:      Minimum melting temperature (Disabled in default mode) 

[-maxTm]:      Maximum melting temperature (Disabled in default mode)

[-rangeTm]:    The length of melting temperature interval (Default value is 10)

[-oligCon]:    Oligo concentration in nM (Default value is 1000)

[-saltCon]:    Salt concentration in mM (Default value is 75)

[-secStr]:     Enabling secondary structure checking


Example
=======

For the input dataset "ecoli.fsa" the commands for running BOND with default parameters on linux, mac, and windows are:

1. For Linux:
	./BOND_linux ecoli.fsa ecoli.bond

2. For Mac:
	./BOND_mac ecoli.fsa ecoli.bond

3. For Windows (64-bit):
	BOND_windows ecoli.fsa ecoli.bond


If the user would like to design oligos with these parameters:

- length: 60
- sequence identity percentage: 80
- maximum consecutive match: 16
- enabling secondary structure checking

then the linux command is:

	./BOND_linux ecoli.fsa ecoli.bond -length 60 -seqSim 80 -maxMatch 16 -secStr 



=============================================================================================================
How to run EVAL (the program for evaluating the specificity); PRE_EVAL prepares the oligo file for evaluation
=============================================================================================================

	./PRE_EVAL_<platform> <BONDOligos> <plainOligos>
	./EVAL_<platform> <inputSequences> <plainOligos> <outputReport>

Required parameters:

< BONDOligos > - the file containing the oligos designed by BOND

< plainOligos > - a file containing only the oligos, one per line, and no additional information

< inputSequencess > - the input file containing the gene sequences in FASTA format

< outputReport > - the output file containing the evaluation report (bad oligos, good oligos, and statistics)


Example
=======

To evaluate the oligos designed by BOND "ecoli.bond" for the dataset "ecoli.fsa" the linux commands are:

        ./PRE_EVAL_linux ecoli.bond ecoli.plain
        ./EVAL_linux ecoli.fsa ecoli.plain ecoli.report



=====================================================================================================
How to run MAX (the program for estimating the maximum number of oligos, needed to evaluate coverage)
=====================================================================================================

	./MAX_<platform> <inputSequences> <outputReport> 

Required parameters:

< inputSequencess > - the input file containing the gene sequences in FASTA format

< outputReport > - the output file containing the estimated maximum number of oligos


Example
=======

To estimate the maximum number of oligos for the dataset "ecoli.fsa" the linux command is:

        ./MAX_linux ecoli.fsa ecoli.max

# Datasets
- [arabidopsis.fsa](http://www.csd.uwo.ca/~ilie/BOND/arabidopsis.fsa)
- [bee.fsa](http://www.csd.uwo.ca/~ilie/BOND/bee.fsa)
- [celegans.fsa](http://www.csd.uwo.ca/~ilie/BOND/celegans.fsa)
- [chicken.fsa](http://www.csd.uwo.ca/~ilie/BOND/chicken.fsa)
- [drosophila.fsa](http://www.csd.uwo.ca/~ilie/BOND/drosophila.fsa)
- [ecoli.fsa](http://www.csd.uwo.ca/~ilie/BOND/ecoli.fsa)
- [human.fsa](http://www.csd.uwo.ca/~ilie/BOND/human.fsa)
- [maize.fsa](http://www.csd.uwo.ca/~ilie/BOND/maize.fsa)
- [mouse.fsa](http://www.csd.uwo.ca/~ilie/BOND/mouse.fsa)
- [plasmodium.fsa](http://www.csd.uwo.ca/~ilie/BOND/plasmodium.fsa)
- [rice.fsa](http://www.csd.uwo.ca/~ilie/BOND/rice.fsa)
- [yeast.fsa](http://www.csd.uwo.ca/~ilie/BOND/yeast.fsa)
- [zebrafish.fsa](http://www.csd.uwo.ca/~ilie/BOND/zebrafish.fsa)
- [mouserna.fsa](http://www.csd.uwo.ca/~ilie/BOND/mouserna.fsa)
- [mouse1421.fsa](http://www.csd.uwo.ca/~ilie/BOND/mouse1421.fsa)


# CITE

If you use BOND, please cite:

L. Ilie, H. Mohamadi, G.B. Golding, W.F. Smyth, [BOND: Basic OligoNucleotide Design, BMC Bioinformatics](http://link.springer.com/article/10.1186/1471-2105-14-69) 14 (2013) 69.

