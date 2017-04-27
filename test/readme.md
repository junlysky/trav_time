# A program to compute first arrival times in layered velocity model   

## usage
	1. the input file: model
           model format:  thickness vs vp rho
	2. get the binary program "trav" into the shell path in your computer
	3. perl trav.pl -Mmodel/depth distance
	4. the main out file is vps.dat, in whitch are: distance depth tp $tp ts $ts
