"""PSSM File Parser

PSSM Parser parse protein position specific matrix

Example:
	from pssm import pssmparser
	p = PSSMParser("/home/welhefna/scope2.06New/fold2vector/database/PSSM/4INKA.pssm")
	from pprint import pprint
	print (p.sequence_to_pssm('VV'))
	print len(p.sequence_to_pssm('VV'))

Note:
	None

TODO:
	*Improve parsing performace depend on the research goals.

"""

import os
import re

class PSSMParser(object):
	"""PSSM Parser
	
	Attributes:
		aa_sequence (str) :  protein AA sequence format. 
		lpssmc (list) : Last position-specific scoring matrix computed.
		woprd (list) : weighted observed percentages rounded down.
		ipp (float) : information per position.
		rwgrmtp (str) : relative weight of gapless real matches to pseudocounts.
		standard_ungapped  (float) : standard ungapped
		standard_gapped  (float) : standard gapped.
		psi_ungapped  (float) : psi ungapped.
		psi_gapped  (float) : psi gapped.
		file_path (str) : PSSM file path.
	
	"""
	
	__slots__ = [ '_aa_sequence', '_lpssmc', '_woprd', '_ipp', '_rwgrmtp', '_standard_ungapped', '_standard_gapped', '_psi_ungapped', '_psi_gapped', '_file_path']
	
	def __init__(self, file_path):
		"""Initialize PSSMParser Object
		
		Initialize PSSMParser object attributes from the given file path.
		
		Args:
			file_path (str) : target file full path.
		Kwargs :
			None.
		Returns:
			None.
		Rasis:
			IOError : file not available.
			
		"""
		
		self._file_path = file_path
		self._aa_sequence = [] #AA sequence format
		self._lpssmc = []  #Last position-specific scoring matrix computed Lx20 matrix
		self._woprd = [] #weighted observed percentages rounded down Lx20 matrix
		self._ipp = [] #information per position
		self._rwgrmtp = [] #relative weight of gapless real matches to pseudocounts
						 #K         Lambda
		self._standard_ungapped = None   #Standard Ungapped    0.1384     0.3187
		self._standard_gapped = None     #Standard Gapped      0.0410     0.2670
		self._psi_ungapped = None        #PSI Ungapped         0.1516     0.3179
		self._psi_gapped = None          #PSI Gapped           0.0464     0.2670
	
		
		if os.path.exists(file_path):
			with open(file_path,'r') as pssm_file:

				t_file_lines = pssm_file.readlines()
				t_aa_sequence = [] 
				for line in t_file_lines:
					line = self.clean(line)
					if len(line) >= 44:
						self._aa_sequence.append(line[2])
						self._lpssmc.append(map(float,line[3:3+20]))
						self._woprd.append(map(float,line[23:23+20]))
						self._ipp.append(float(line[-2]))
						self._rwgrmtp.append( float(line[-1]) if line[-1] != "inf" else "inf" )
				
				t_sgp = self.clean(t_file_lines[-4])
				self._standard_ungapped = { 'K': float(t_sgp[-2]), 'Lambda': float(t_sgp[-1]) }
				
				t_sgp = self.clean(t_file_lines[-3])
				self._standard_gapped = { 'K': float(t_sgp[-2]), 'Lambda': float(t_sgp[-1]) }
				
				t_sgp = self.clean(t_file_lines[-2])
				self._psi_ungapped = { 'K': float(t_sgp[-2]), 'Lambda': float(t_sgp[-1]) }
				
				t_sgp = self.clean(t_file_lines[-1])
				self._psi_gapped = { 'K': float(t_sgp[-2]), 'Lambda': float(t_sgp[-1]) }
				
		else:
			raise IOError( " File {} not valid path or not available ".format(self._file_path))
			
	def clean(self, string):
		"""Clean
		
		clean repeated spaces and carrage.
		
		Args:
			string (str) : input string to be cleaned.
		Kwargs:
			None.
		Returns:
			string (str) : cleaned string
		Raises:
			None.
			
		"""
		
		t_string = string.rstrip()
		t_string = re.sub(' +',' ',t_string)
		t_string = t_string.split(' ')
		return t_string
		
	def sequence_to_pssm(self, sequence):
		"""Sequence to PSSM
		
		return the PSSM information of the query sequence.
		
		Args: 
			sequence (str) :  AA format of protein amino-acid in string format.
		Kwargs:
			None.
		Returns:
			t_pssm (dict) :  PSSM information <key, value>, where key is PSSM attributte, and value is PSSM information.
		Raises:
			None.
			
		"""
		
		t_aa_sequence = "".join(self._aa_sequence)
		count = t_aa_sequence.count(sequence)
		
		t_pssm = []
		
		sequence_length = len(sequence)
		index = 0
		for c in range(count):
			index = t_aa_sequence.find(sequence, index)
			t ={
					"aa_sequence" : self._aa_sequence[index : index + sequence_length],
					"lpssm" : self._lpssmc[index : index + sequence_length],
					"woprd"	: self._woprd[index : index + sequence_length],
					"ipp"	: self._ipp[index : index + sequence_length],
					"rwgrmt" : self._rwgrmtp[index : index + sequence_length],
					"index" : range(index , index + sequence_length),
			}
	
			index += sequence_length
			t_pssm.append(t)
			
		return t_pssm

		
			
			
			
			
if __name__ == "__main__":
	pass