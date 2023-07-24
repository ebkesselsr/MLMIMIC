Supplementary material / reproducible files for the manuscript:
Title: "A comparison of the multilevel MIMIC model to the multilevel regression and mixed ANOVA model for the estimation and testing of a cross-level interaction effect: a simulation study".

Authors: Kessels, R., Moerbeek, M.
Code was written by Kessels, R.
In case of questions or comments please contact robkessels5@gmail.com

IMPORTANT: to obtain the results of the manuscript, ONLY the program "Simulation study code.R" needs to be ran.
(see below the details).

This folders contains the following data and files that can be used to reproduce all results of the manuscrtipt. 

./Mplus/  
	A subfolder containing the Mplus code to fit a ML-MIMIC model.  
	A working Mplus version is required. Make sure the R program "Simulation study code.R" is able to access this program. 
	The R-code creates a .dat Mplus data file and stores that in the same folder as the Mplus program. Mplus then reads
	this .dat file to estimate the parameters. These parameters are then re-loaded in R. This goes automatically when
	running the Simulation code below. 

./R-code/  
	-Simulation study code.R  
     	The simulation code to reproduce the results presented in the manuscript.
	To obtain the results of the manuscript, ONLY this program needs to be ran.
	The code takes a long time to run over all 1000 iterations. The results of each iteration is saved in three
	separate .txt files (for each model separately). These .txt files are saved in ./Output data/.   

	-Create tables.R  
 	R-program to reproduce the tables 1-7 of the supplementary material. This function reads the .txt tables in ./Output data/
	and creates LaTeX tables.   
	The tables were manually adjusted to fit the main manuscript.

	-Power analysis ML-MIMIC.R  
 	R-program to perform a power analysis with a ML-MIMIC model, in case a researcher is interested. THIS CODE WAS NOT USED FOR THIS MANUSCRIPT.   

./Output data/    
	-outputmlmimic.txt, outputml.txt, outputaov.txt  
	Simulation output for for the ML-MIMIC model, univariate multilevel model, and ANOVA model.  
	In this folder, the results of each iteration are saved during the simulation study.   
	Now, all output is already there. 
