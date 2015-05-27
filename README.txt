SWAlign.jar is a simple and robust alignment program for protein
sequence-to-sequence alignment based on the standard Smith-Waterman
dynamic programming algorithm. The mutation matrix is from BLOSUM62
with gap openning penaly=-11 and gap extension panalty=-1. 
The program is written by Renxiang Yan at Dr. Yang Zhang's laboratory 
in the Univerisity of Michigan.

Usage:
java -jar SWAlign.jar F1.fasta F2.fasta  (align two sequences in fasta file)
java -jar SWAlign.jar F1.pdb F2.pdb    1 (align two sequences in PDB file)
java -jar SWAlign.jar F1.fasta F2.pdb  2 (align 1 in fasta and 2 in pdb)
java -jar SWAlign.jar GKDGL EVADELVSE  3 (align sequences typed by keyboard)
java -jar SWAlign.jar GKDGL F.fasta    4 (align Seq-1 by keyboard and 2 in fasta)
java -jar SWAlign.jar GKDGL F.pdb      5 (align Seq-1 by keyboard and 2 in pdb)

Chanegs were made to display residue mappings and output in JSON format.