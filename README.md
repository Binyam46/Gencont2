# Gencont2

1.	Introduction

a.	Purpose 

GENCONT-2 is a program to calculate optimal genetic contributions for selection candidates so that it maximizes the genetic gain whilst controlling the rate of inbreeding (managing genetic diversity). Inbreeding is managed by restricting the average relationship among selected parents because in a population inbreeding increases by 0.5 times the average relationship. GENCONT-2 uses a method that optimizes genetic contributions of selection candidates constraining on a predefined rate of inbreeding as described in Meuwissen (1997). GENCONT-2 is based on an iterative algorithm implementing the aforementioned method as presented in Dagnachew and Meuwissen (2014). 

b.	Features
GENCONT-2 consists of a single program and all the required information to run the program, input files and their layout is specified in a parameter file. GENCONT-2 also incorporates calculation of optimum genetic contributions in the case of population with overlapping generation as described in Meuwissen and Sonesson (1998) and breeding schemes with multiple selection stages. The latter feature allows optimizing more than one selection stages simultaneously (for example, pre-selection of young animals to enter test station and final selection of sire after progeny testing).   

The basic inputs are data file containing EBV of selection candidates and a pedigree file for computation of relationship among selection candidates. Some advance running options are described in the parameter file section. In the case of genomic selection, where animals are selected based on GEBVs instead of traditional EBVs, genomic relationship (G) matrix should be used in place of pedigree relationship (A) matrix (Sonesson, et.al., 2012). GENCONT-2 can take genomic relationship matrix (G) constructed outside the program and will use it to constrain relationship.  In some cases, some candidate groups are genotyped and have GEBVs, and other groups are ungenotyped and have traditional EBVs. In situation like this, a hybrid of genomic and pedigree relationship matrices (i.e. H matrix) can be constructed (Legarra et al., 2009) and could be used in GENCONT-2 to constrain inbreeding. 

2.	Availability
More to come …

3.	Installation and running GENCONT-2
GENCONT-2 is distributed as pre-compiled executable file. The main target operating system is Linux. In addition, executable files for Windows environments can be available. 

GENCONT-2 is expected to run form a command line interface, i.e. you need to open a ‘terminal’ window and type in the appropriate command, hitting return to start the program running. The general form of the command (under Linux) is:

	./gencont_2	parameter_file 	output_file
	
with the parameter_file is the name of the ‘parameter file’ and the output_file is the name of the ‘output file’ to which GENCONT-2 will write the output. If the output_file is not provided, the default output name is ‘gencont_2.out’. 

4.	Input files 

The minimal requirement of GENCONT-2 are parameter file, data file and pedigree file. 

Parameter file: 
	This file contains all information that specifies everything GENCONT-2 needs to know about the input files and their 	outline and user running preferences (see parameter file section).

Data file: 
	This file should at least contains ID and EBV or GEBV of the selection candidates. If more than one group (sex) of animals are included in the data file, then ‘sex’ of the 	individuals should be specified. Age-class and availability of the animals should be included in the file in the case of overlapping generation and multiple selection stages. Further, some maximum and minimum contributions that a specific animal can have, and fraction of offspring already allocated to a specific animal could be also included in the 	data file. The way to read the data file into GENCONT-2 is described under the parameter file section. 

Pedigree file:
	Almost in all cases (the exception is when ‘readin’ option is activated, see parameter file under ‘Amatrix’ keyword), a pedigree file is required to calculates relationship among selection candidates. GENCONT-2 requires the pedigree file to be 	renumbered from 1 to N. The individuals ID provided in the data file should match the ID in the pedigree file. 
