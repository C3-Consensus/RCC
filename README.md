# RCC

The RCC is a 26-dimensional vector, each value corresponding to a feature from the topology of the protein: the number of residues of a given cluster class.

Options to run code:

	-pdb			      <file to be processed>
	-chain			    <chain in the pdb file to be processed>
	-pathTomita		  <path to tomita executable file>
	-lateralChain 	<yes/no (no)>
	-minDist 		    <minimum contact distance (0)>
	-maxDist 		    <maximum contact distance (7)>
	-fileList 		  <csv file with filename, chain>
	-pathPdbFolder 	<path to pdb database folder>
  
The program will output the rcc on the terminal with the format:
<Protein> [<RCC separated by commas>] <Number of residues considered> <Number of residues missing>  <Ids of the missing residues>
e.g. 103lA00	[8, 10, 7, 0, 4, 2, 14, 5, 0, 0, 0, 1, 1, 11, 80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4]	157	5	35|36|37|38|39
For the protein 103l its rcc is [8, 10, 7, 0, 4, 2, 14, 5, 0, 0, 0, 1, 1, 11, 80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4], it has 157 and there are 5 residues missing from the pdb file (35, 36, 37, 38, 39). It is unexpected that these residues are missing so it is reported.
