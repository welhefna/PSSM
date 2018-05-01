PSSM Parser Module Repository
=============================

PSSMParser Module is custom implementation of PSSM file reader for residue wise protein information.


What is a PSSM?
---------------

A PSSM, or Position-Specific Scoring Matrix, is a type of scoring matrix used in protein BLAST searches in which amino acid substitution scores are given separately 
for each position in a protein multiple sequence alignment. Thus, a Tyr-Trp substitution at position A of an alignment may receive a very different score than the same 
substitution at position B. This is in contrast to position-independent matrices such as the PAM and BLOSUM matrices, in which the Tyr-Trp substitution receives the same 
score no matter at what position it occurs.

PSSM scores are generally shown as positive or negative integers. Positive scores indicate that the given amino acid substitution occurs more frequently in the alignment 
than expected by chance, while negative scores indicate that the substitution occurs less frequently than expected. Large positive scores often indicate critical functional
residues, which may be active site residues or residues required for other intermolecular interactions.

PSSMs can be created using PSI-BLAST, which finds similar protein sequences to a query sequence, and then constructs a PSSM from the resulting alignment. Alternatively, 
PSSMs can be retrieved from the NCBI CDD database, since each CD is represented by a PSSM that encodes the observed substitutions in the seed alignments. These CD records 
can be found either by text searching in Entrez Conserved Domains or by using RPS-BLAST (Reverse Position-Specific BLAST), also known as CD-Search, to locate these domains on 
an input protein sequence

Usage
--------------
	
	from pssm import pssmparser as pssm_reader

Example
---------------

	try:
	
		p = PSSMParser("/home/welhefna/scope2.06New/fold2vector/database/PSSM/4INKA.pssm")

		
	except:
		#exception handler
		print "Exception is raised !", sys.exc_info()
	else:
		#print results of 
		print p.sequence_to_pssm('VVACQ')
		print len(p.sequence_to_pssm('VV'))

		
Setup
---------------
	
	python setup.y

