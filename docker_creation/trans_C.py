from base import *
from rand_walk import *
	
def read_user_inputs():
	print(f"Starting trans-C.\n")
	
	global fpath, a, bin_size, seed, chrom_sizes, bin_map, chr_map, nbins
	
	chrom_sizes_file = sys.argv[2]
	res = int(sys.argv[3])
	seed_file = sys.argv[4]
	fpath = sys.argv[5]
	
	# if the user specified an alpha, use it
	if len(sys.argv) == 7:
		alpha = float(sys.argv[6])
	# else use the default of 0.5
	else:
		alpha = 0.5
	
	process_inputs(chrom_sizes_file, fpath, res)
	seed = read_seed(seed_file)
	
	a = parse_matrix(sys.argv[1], fpath, skip=False, store=True)
	
	do_random_walk(a, alpha, seed, fpath)
	
	clean_up()
	
	print(f"\ntrans-C done.\n")

if __name__ == "__main__":
	read_user_inputs()

