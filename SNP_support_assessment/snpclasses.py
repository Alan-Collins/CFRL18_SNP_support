#!/usr/bin/env python3

class SAM_data(object):
	"""stores columns of SAM entry as attributes"""
	def __init__(self, object):
		self.qname = object.split('\t')[0]
		self.flag = object.split('\t')[1]
		self.rname = object.split('\t')[2]
		self.pos = int(object.split('\t')[3])
		self.mapq = int(object.split('\t')[4])
		self.cigar = object.split('\t')[5]
		self.rnext = object.split('\t')[6]
		self.pnext = object.split('\t')[7]
		self.tlen = object.split('\t')[8]
		self.seq = object.split('\t')[9]
		self.qual = object.split('\t')[10]
		self.ln = len(self.seq)
		self.end = self.pos + self.ln
		self.mod_seq = ''
		self.mod_qual = ''
		self.refreshed = False

	def refresh(self):
		if not self.refreshed:
			self.seq = self.mod_seq
			self.qual = self.mod_qual
			self.ln = len(self.seq)
			self.end = self.pos + self.ln
			self.refreshed = True

class snp():
	"""Simple class to store contig ID and position info for SNPs"""
	def __init__(self, line):
		self.contig = line.split()[0].split(':')[0]
		self.position = line.split()[1]
		

class Read_padding():
	"""Store information about how much to pad reads when printing them in order to maintain register between lines"""
	def __init__(self):
		self.read_name = ''
		self.seq = ''
		self.qual = ''
		self.cigar = ''
		self.pad = 0 # 0: no pad, 1: left pad, 2: right pad, 3: both sides pad
		self.padl = 0
		self.padr = 0
		self.warning = ''
		
	def report_pad_seq(self):
		for i in self.qual:
			if ord(i) < 53:
				self.warning = '!!!WARNING: Low quality base call in sequence!!!'

		pseq = '-' * self.padl + self.seq + '-' * self.padr 
		pqual = '-' * self.padl + self.qual + '-' * self.padr
		return '\t'.join([self.read_name, pseq, pqual, self.cigar, self.warning])