import numpy as np
import itertools as it
import multiprocessing as mp

import os
import sys
import math
import random
import time
import gzip
import glob
import pickle

from scipy.linalg import eig
from data_processing import *


# --------------------------------------------------------------------------- #
# Basic timer
# --------------------------------------------------------------------------- #

def print_time_taken(start_time, task_name=""):
	t = time.time() - start_time
	st = "done " + task_name + " after = "
	if t < 60:
		st += (str(round(t,1)) + " sec")
	elif t < 3600:
		st += (str(round(t/60,1)) + " min")
	else:
		st += (str(round(t/3600,2)) + " h")
	print(st)
	return st
	



# --------------------------------------------------------------------------- #
# Read and Write matrix
# --------------------------------------------------------------------------- #

# basic read matrix from our format
def read_matrix(fname, symmetric = False):
	start_time = time.time()
	
	matrix = np.zeros((nbins,nbins))
	
	for line in open(fname):
		words = line.rstrip().split("\t")
		i = int(words[0])
		j = int(words[1])
		w = float(words[2])
		 
		matrix[i][j] = w
		
		if symmetric:
			matrix[j][i] = w
	
	print_time_taken(start_time, "reading input matrix")
	if symmetric:
		print("  Treated the matrix as SYMMETRIC!!")
	return matrix


# write down the matrix in basic format: bin_i bin_j val
def write_matrix(a, fname, write_vals_above = 0, upper_half = False):
	fo = open(fname, "w")
	
	for i in range(a.shape[0]):
		for j in range(a.shape[1]):
			if upper_half and j > i:
				continue
			if a[i][j] > write_vals_above:
				fo.write(str(i) + "\t" + str(j) + "\t" + str(a[i][j]) + "\n")
	fo.close


def pickle_write(fname, obj):
	with open(f"{fname}", 'wb') as f:
		    pickle.dump(obj, f)

def pickle_read(fname):
	with open(f"{fname}", 'rb') as f:
		    return pickle.load(f)


# --------------------------------------------------------------------------- #
# Read basic genome information
# --------------------------------------------------------------------------- #

def read_chr_lengths(fname):
	chr_length, chr_names = {}, []

	for line in open(fname):
		words = line.rstrip().split("\t")
		chr = words[0]
		length = int(words[1])
		
		chr_length[chr] = length
		chr_names.append(chr)
		
	return (chr_length, chr_names)


def read_chr_lengths_human(fname):
	chr_length, chr_names = {}, []
	for line in open(fname):
		words = line.rstrip().split("\t")
		chr = words[0]
		length = int(words[1])
		
		chr_length[chr] = length
		chr_names.append(chr)
		
	return (chr_length, chr_names)

# --------------------------------------------------------------------------- #
# Read the bin map
# --------------------------------------------------------------------------- #

# read bin map
def read_bin_map(fname):
	chr_map, bin_map = {}, {}
	
	for line in open(fname):
		words = line.rstrip().split("\t")
		chr   = words[0]
		start = int(words[1])
		end   = int(words[2])
		bin   = int(words[3])
		
		if chr not in chr_map:
			chr_map[chr] = []
			
		chr_map[chr].append([start, end, bin])
		bin_map[bin] = (chr, start, end)
	
	nbins = max(bin_map.keys()) + 1
	print("nbins = " + str(nbins))
	return (chr_map, bin_map, nbins)


def read_seed(fname):
	seed_bins = []
	for line in open(f"{fname}"):
		words = line.rstrip().split("\t")
		mid = int((int(words[1]) + int(words[2]))/2)
		
		seed_bins.append(get_bin_from_chr_loc(words[0], mid))
		
	return seed_bins
		

# --------------------------------------------------------------------------- #
# Find the bin for a given genomic location
# --------------------------------------------------------------------------- #

# returns the bin corresponding to the given genomic coordinates
def get_bin_from_chr_loc(chr, loc):
	x = int(loc/bin_size)
	
	if x >= len(chr_map[chr]):
		t = chr_map[chr][-1]
		if loc >= t[1] and loc < (t[1] + bin_size):
			return t[2]
		else:
			raise ValueError("location " + str(loc) + " is outside of chromosome " + str(chr))
	else:
		return chr_map[chr][x][2]


# --------------------------------------------------------------------------- #
# Generate mcool file from our npy matrix
# --------------------------------------------------------------------------- #

def generate_mcool(a, fname):
	# See: https://cooler.readthedocs.io/en/latest/cli.html#cooler-load
	
	# 1. write the matrix into COO format
	write_matrix(a, f"{T}{fname}.txt", 0, upper_half = True)

	# 2. call cooler
	os.system(f"cooler load -f coo {P}Data/genome/chr_lengths_{organism}.txt:{res} {T}{fname}.txt {T}{fname}.cool --field count:dtype=float")
	os.system(f"cooler zoomify {T}{fname}.cool")
	os.remove(f"{T}{fname}.txt")
	os.remove(f"{T}{fname}.cool")	


def write_array(array, fname):
	fo = open(fname, "w")
	for i in range(len(array)):
		fo.write(f"{i}\t{array[i]}\n")




# --------------------------------------------------------------------------- #
# Donut Filter
# --------------------------------------------------------------------------- #

# helper to calculate the mean of the donut-shaped neighborhood
def donut_mean(values, structure):
	return np.mean(values[structure])

def donut_structure(matrix, inner_radius=2, outer_radius=40):
	from scipy.ndimage import generic_filter
	
	# create a donut-shaped structuring element
	structure = np.zeros((2*outer_radius + 1, 2*outer_radius + 1), dtype=bool)
	for i in range(2*outer_radius + 1):
		for j in range(2*outer_radius + 1):
			distance = np.sqrt((i - outer_radius)**2 + (j - outer_radius)**2)
			if inner_radius <= distance <= outer_radius:
				structure[i, j] = True

	# apply the donut filter to calculate the mean of the donut neighborhood for each element
	mean_matrix = generic_filter(matrix, donut_mean, footprint=structure, mode='constant', cval=np.nan)

	# normalize the original matrix by dividing by the mean matrix, handling potential division by zero
	m = np.where(mean_matrix != 0, matrix / mean_matrix, np.nan)

	return m




# --------------------------------------------------------------------------- #
# Get the zero and non-zero rows, zero cis or special regions
# --------------------------------------------------------------------------- #

def find_zero_and_non_zero_rows(a):
	zero_bins, non_zero_bins = [], []
	
	zr = np.all((a == 0), axis = 1)
	for i in range(len(zr)):
		if zr[i]:
			zero_bins.append(i)
		else:
			non_zero_bins.append(i)

	print("found number of rows with only zeros = " + str(len(zero_bins)))
	return (zero_bins, non_zero_bins)

def zero_cis_contatcs(a):
	for (i,j) in it.combinations(list(range(nbins)), 2):
		if bin_map[i][0] == bin_map[j][0]:
			a[i][j] = 0
			a[j][i] = 0
			
	return a

def zero_special_regions(a):
	for i in range(len(a)):
		for j in range(nbinsh, len(a)):
			a[i][j] = 0
			a[j][i] = 0

	return a
	

# return the weight of a given clique
def find_clique_weight(clique, matrix):
	weight = 0
	for i in clique:
		for j in clique:
			if bin_map[i][0] == bin_map[j][0]:
				continue
							
			weight += matrix[i][j]
	
	return weight/2



# --------------------------------------------------------------------------- #
# report basic statistics for the Hi-C data
# --------------------------------------------------------------------------- #

def report_HiC_stats():
	cis_space, trans_space, num_cis, num_trans, empty_trans, empty_cis = 0, 0, 0, 0, 0, 0
	for i in range(nbins):
		for j in range(i, nbins):
			if bin_map[i][0] == bin_map[j][0]:
				num_cis += a[i][j]
				cis_space += 1
				if a[i][j] == 0:
					empty_cis += 1
			else:
				num_trans +=  a[i][j]
				trans_space += 1
				if a[i][j] == 0:
					empty_trans += 1
	
	total_space = cis_space + trans_space	
	total_reads = num_cis + num_trans
	
	print(f"total_space =\t{total_space}")
	print(f"cis_space = \t{cis_space}")
	print(f"trans_space = \t{trans_space}")
	print(f"trans_space_ratio = \t{trans_space/(trans_space + cis_space)}")
	print(f"num_cis = \t{num_cis}")
	print(f"num_trans = \t{num_trans}")
	print(f"total_reads = \t{total_reads}")
	print(f"trans_fract = \t{num_trans/total_reads}")
	print(f"empty_cis = \t{empty_cis}")
	print(f"empty_trans = \t{empty_trans}")
	print(f"frac_empty_trans = \t{empty_trans/trans_space}")




# --------------------------------------------------------------------------- #
# Cytoscape Visualization: write .gml file to use as input to Cytoscape 
# --------------------------------------------------------------------------- #

def create_gml_file_Cytoscape(fname, fn = "visual"):
	order = read_bins(f"{fname}")
	
	num_nodes = 40
	
	fo = open(f"{fout}clique_{fn}.gml", "w")
	fo.write(f"graph [\ncomment \"This is a sample graph\"\ndirected 0\nweighted 1\n")
	
	# we need to hande chr names nicely
	h_chr = {}
	for i in range(num_nodes):
		h_chr[i] = bin_map[order[i]][0].replace('chr', '')
	
	# write down the nodes
	for i in range(num_nodes):
		# label Mb position
		mb = round(float(bin_map[order[i]][1])/1000000,1)
		
		# order the nodes by bin_id with 0 in front to make it alphabetical
		lid = order[i]
		if lid < 10:
			lid = f"0000{lid}"
		elif lid < 100:
			lid = f"000{lid}"
		elif lid < 1000:
			lid = f"00{lid}"
		elif lid < 10000:
			lid = f"0{lid}"
		
		# chromosome order	
		ch = h_chr[i]
		if "X" in ch or "Y" in ch:
			ch = f"chr \"{ch}\""
		elif float(ch) < 10:
			ch = f"chr \"a0{ch}\""
		else:
			ch = f"chr \"a{ch}\""
		
		fo.write(f"node [\nid {i+1}\nlabel \"{mb}\"\nweight 7\norder_by \"{lid}\"\n{ch}\n]\n")
	
	# write down the edges
	for i in range(num_nodes):
		for j in range(i+1, num_nodes):
			chr1 = bin_map[order[i]][0]
			chr2 = bin_map[order[j]][0]
			
			if chr1 != chr2:
				x = a[order[i]][order[j]]
							
				fo.write(f"edge [\nsource {i+1}\ntarget {j+1}\nlabel \"{i}-{j}\"\nweight {x}\ncolor {x}\n]\n")
	fo.write(f"]\n")





# --------------------------------------------------------------------------- #
# read hashes
# --------------------------------------------------------------------------- #

def read_hashes(fname1, fname2):
	global hsh_entrz_ensmbl, hsh_ensmbl_entrz, hsh_entrz_human, hsh_human_entrz
	hsh_entrz_ensmbl, hsh_ensmbl_entrz, hsh_entrz_human, hsh_human_entrz = {}, {}, {}, {}
	
	for line in open(fname1):
		words = line.rstrip().split("\t")
		hsh_entrz_ensmbl[words[0]] = words[1]
		hsh_ensmbl_entrz[words[1]] = words[0]
	
	for line in open(fname2):
		words = line.rstrip().split("\t")
		hsh_human_entrz[words[0]] = words[1]
		hsh_entrz_human[words[1]] = words[0]
		
	# Probably don't return them but read them at the begining
	return (hsh_entrz_ensmbl, hsh_ensmbl_entrz, hsh_entrz_human, hsh_human_entrz)


def get_gene_name_from_ensmbl(e):
	if e in hsh_ensmbl_entrz and hsh_ensmbl_entrz[e] in hsh_entrz_human:
		return hsh_entrz_human[hsh_ensmbl_entrz[e]]
	else:
		return ""


# --------------------------------------------------------------------------- #
# read gff file
# --------------------------------------------------------------------------- #

def read_GFF_file(fname):
	global gff_dict_genes, gff_dict_bins
	gff_dict_genes, gff_dict_bins = {}, {}
	gene_types, feature_types = [], []
	start_time, num_lines = time.time(), 0
	
	for line in gzip.open(fname, "rt"):
		num_lines += 1
		words = line.rstrip().split("\t")
		
		# the genome sequence is below- we don't need it right now
		if words[0] == "##FASTA":
			print("encountered ##FASTA")
			break
			
		if words[0][0] == "#": continue 
		
		# feature_types.append(words[2]) # explore: what features are in this gff?
		
		if words[2] != "gene": continue
		
		gene      = words[8].split(";")[0][3:].split(".")[0]
		gene_type = words[8].split(";")[2][10:] # gene_types.append(gene_type) # explore: what types of genes are in this gff?
		tss       = int(words[3] if words[6] == "+" else words[4])
		
		
		chr = words[0]
		
		# if chrM or PF_M7661: skip
		if chr not in chr_map:
			continue
		

		if gene_type != "protein_coding":
			continue
		
		# finally, get the bin and populated the hashes
		bin  = get_bin_from_chr_loc(chr,tss)
		
		gff_dict_genes[gene] = (chr, tss, bin)
		
		if not bin in gff_dict_bins: gff_dict_bins[bin] = []
		gff_dict_bins[bin].append(gene)


	print_time_taken(start_time, "reading gff file with " + str(len(gff_dict_genes)) + " genes")
	return (gff_dict_genes, gff_dict_bins)




# --------------------------------------------------------------------------- #
# matrix normalization 
# --------------------------------------------------------------------------- #

# this functions carries the normalization steps from raw Hi-C matrix to one
# on which trans-C runs. 
# See im_matrix() which is now used instead to allow Sara to make rapid progress
def normalzie_matrix(a, fpath):
	pass


# Employ the ICE-normalizaation procedure
def ice_normalize_matrix(a, fpath, bias_name=""):
	from iced import normalization
	from iced import filter
	
	counts = filter.filter_low_counts(a, percentage=0.05)
	(normed, bias) = normalization.ICE_normalization(counts, output_bias=True)
	
	# the biases from ice
	if len(bias_name) > 1:
		#fo = open(bias_name, "w")
		for i in bias:
			i[0] = f"{i[0]}\n"
	
	return normed
	

# normalize the Hi-C matrix using binomial distribution to account for 
# chromosomal territories 
def normalize_chr_territories(a):	
	# number of possible bins per chromosome
	chr_loci = {}
	for c in chr_map.keys():
		chr_loci[c] = len(chr_map[c])
		
	from scipy.stats import binomtest
	
	# for each pair of chromosome
	for c1 in chr_map.keys():
		for c2 in chr_map.keys():
			if c1 == c2: continue
			num_loci_pairs, num_counts = 0, 0
			
			# obtain the necessary info for binomial function
			for e1 in chr_map[c1]:
				for e2 in chr_map[c2]:
					num_loci_pairs += 1
					num_counts += a[e1[2]][e2[2]]
					
			# use binomial distribution
			for e1 in chr_map[c1]:
				for e2 in chr_map[c2]:
						a[e1[2]][e2[2]] = -math.log(binomtest(int(a[e1[2]][e2[2]]), int(num_counts), float(1/num_loci_pairs), alternative='greater').pvalue + math.exp(-10))
	return a


# Apply donut filter normalization
def donut_filter(a, inner_radius=2, outer_radius=40):
	# loop over the Hi-C matrix
	for i in range(nbins):
		for j in range(nbins):
			
			# define boundaries that do not go outside the matrix
			x1 = max(0, i - outer_radius)
			y1 = max(0, j - outer_radius)
			x2 = min(nbins - 1, i + outer_radius)
			y2 = min(nbins - 1, j + outer_radius)
			
			# define donut stats
			num_loci, sum_counts = 0, 0
			
			# go over all points in the donut structure
			for i1 in range(x1, x2):
				for j1 in range(y1, y2):
					# do not count bins along the center axis of the donut
					if i1 == i or j1 == j:
						continue
					
					# if by any chance the donut spreads over chromosome boundaries
					if bin_map[i][0] != bin_map[i1][0] or bin_map[j][0] != bin_map[j1][0]:
						# skip those
						continue
					
					distance = np.sqrt((i - i1)**2 + (j - j1)**2)
					
					# points within the donut
					if inner_radius <= distance <= outer_radius:
						num_loci += 1
						sum_counts += a[i1][j1]
						
			# normalize the matrix entry
			if num_loci > 0 and sum_counts > 0:
				dn = sum_counts/num_loci
				a[i][j] = a[i][j]/dn
					
	return a


# --------------------------------------------------------------------------- #
# Matrix conversions
# --------------------------------------------------------------------------- #

def parse_matrix(fname, fpath, skip=False, store=False):	
	# if you want to load pre-processed matrix directly and skip the rest
	if skip:
		return np.load(fname)
	
	if "npy" in fname:
		return im_matrix(np.load(fname), fpath, store)
	elif "hic" in fname:
		print(f"please use hicConvertFormat to convert into mcool or npy file format\n")
		raise "not supporting .hic"
	elif "cool" in fname:
		m = convert_mcool_to_npy(fname, bin_size)
		return im_matrix(m, fpath, store)
	else:
		raise "unsupported format"



# wrapper function to potentially extend functionality if needed 
# it allows, given the flag store=True to save the normalized matrix
def im_matrix(a, fpath, store=False):
	m = ice_normalize_matrix(a, fpath)
	
	m = normalize_chr_territories(m)
	
	m = donut_filter(m)
	
	if store:
		np.save(f"{fpath}im_matrix.npy", m)
	return m
	

# --------------------------------------------------------------------------- #
# Read TF bedfile
# --------------------------------------------------------------------------- #

# read a TF input file (.bed.gz)
def read_tf_bed(tf_id, fpath):
	m, length = [0]*nbins, []
	
	with gzip.open(f"{fpath}/{tf_id}.bed.gz",'rb') as file_in:        
		for line in file_in:
			words    = line.decode('utf8').rstrip().split("\t")
			chr_name = words[0]
	
			if chr_name not in chr_map:
				continue
	
			bin = get_bin_from_chr_loc(chr_name, int(words[1]))
	
			m[bin] += 1
	
			length.append(int(words[2]) - int(words[1]))

	print(f"Peaks stats:")
	print(f"np.mean(length) = {np.mean(length)}")
	print(f"np.std(length) = {np.std(length)}")
	print(f"np.max(length) = {np.max(length)}")
	print(f"np.min(length) = {np.min(length)}\n")
	
	return m


# read our order
def read_pi(fname):
	order = []
	for line in open(fname):
		words = line.rstrip().split("\t")
		order.append([int(words[0]), float(words[1])])
	
	return order


# --------------------------------------------------------------------------- #
# Helper functions
# --------------------------------------------------------------------------- #

def read_tf_ids_to_gene_names(fname):
	hsh = {}
	for line in open(fname):
		words = line.rstrip().split("\t")
		hsh[words[0]] = words[1]
		
	return hsh


# get the file names TF idds for all TFs in the given cell line
def get_file_names_in_dir(dir):
	return [os.path.basename(fname).split(".")[0] for fname in glob.glob(f"{dir}")]


def clean_up():
	for used_dir in []:
		files = glob.glob(f"{used_dir}*")
		for f in files:
		    os.remove(f)

# --------------------------------------------------------------------------- #
# READ eCLIP data
# --------------------------------------------------------------------------- #

def read_eclip_file(fname):
	m = [0]*nbins
	with gzip.open(f"{fname}",'rb') as file_in:
		for line in file_in:
			words = line.decode('utf8').rstrip().split("\t")
			chr_name = words[0]

			if chr_name not in chr_map:
				continue

			bin = get_bin_from_chr_loc(chr_name, int(words[1]))
			m[bin] += 1
		
	return m

# get the ids of all eCLIP data tracks we downloaded from ENCODE
def get_eCLIP_ids(fname):
	eclip_ids = []
	for line in open(fname):
		words = line.rstrip()
		if "@@download/" in words:
			eclip_id = words.split("@@download/")[1].split(".")[0]
			eclip_ids.append(eclip_id)
	
	print(f"we read {len(eclip_ids)} eCLIP datasets")
	return eclip_ids


# --------------------------------------------------------------------------- #
#  Print genes in clique
# --------------------------------------------------------------------------- #

def print_genes_in_clique(order, fname, bname=""):
	(gff_dict_genes, gff_dict_bins) = read_GFF_file(fname)
	
	fo = open(f"{T}genes_{bname}.txt", "w")
	for i in order:
		if i not in gff_dict_bins:
			continue
		for gene in gff_dict_bins[i]:
			g = get_gene_name_from_ensmbl(gene)
			if len(g) > 2:
				fo.write(g + "\n")
	fo.close()


def write_genes_in_clique(clique, fname):	
	fo = open(f"{fname}", "w")
	for i in clique:
		if i not in gff_dict_bins:
			continue
		for gene in gff_dict_bins[i]:
			g = get_gene_name_from_ensmbl(gene)
			if len(g) > 1:
				fo.write(g + "\n")
	fo.close()


# --------------------------------------------------------------------------- #
# Process the user inputs and convert them into useful data objects
# --------------------------------------------------------------------------- #

def create_bin_map_user(chrom_sizes_file, fpath, res):
	# read the chromosome sizes supplied by the user and create a bin_map
	global chr_length, chr_names
	(chr_length, chr_names) = read_chr_lengths(chrom_sizes_file)
	
	for bin_size in [res]:

		bin, start, fo = 0, 0, open(f"{fpath}bin_map_{bin_size}.bed", "w")
	
		for chr in chr_names:
			while True:
				end = start + bin_size
			
				fo.write(f"{chr}\t{start}\t{end}\t{bin}\n")
			
				if end < chr_length[chr]:
					start += bin_size
					bin   += 1
				else:
					start = 0
					bin   += 1
					break


def process_inputs(chrom_sizes_file, fpath, res):
	global bin_size, chr_map, bin_map, nbins
	
	bin_size = res
	
	create_bin_map_user(chrom_sizes_file, fpath, res)

	(chr_map, bin_map, nbins) = read_bin_map(f"{fpath}bin_map_{bin_size}.bed")

	
	return (chr_map, bin_map, nbins)



if __name__ == "__main__":
	pass
