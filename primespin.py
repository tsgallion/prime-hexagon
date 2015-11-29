"""  The Number Coil  """
#  List a prime; List mod6 of prime; Compare with another set (offset of mod6)
#  to see spin change, determine cell of each prime.
#
#

"""  This code plots the positons in the central triangle of the prime hexagon.

See accompanying notes and images for what these mean.

"""

import itertools
import numpy as np
import primesieve as ps
import logging
import re

# python 2/3 compatibility stuff
try:
    l = long(1)
except:
    long = int

try:
    zipper = itertools.izip
except:
    zipper = zip

logger = None
log_formatter = None

def setup_loggers():

    global logger, log_formatter
    if logger is not None:
        return
    
    # create logger
    logger = logging.getLogger('prime_hexagon')
    logger.setLevel(logging.INFO)

    # create formatter
    log_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')



setup_loggers()

def primes_from_a_to_b_primesieve(a,b):
    logger.info("starting generating primes from {} to {}".format(a,b))
    out = None
    try:
        print("trying direct numpy prime generation...")
        out = ps.generate_primes_numpy(a,b)
    except:
        pass
    # fall back to classic primesieve python 
    if out is None:
        print("fell back to classic primesieve")
        out = np.array(ps.generate_primes(a,b))        
    logger.info("done generating primes")
    return out


primes_from_a_to_b = primes_from_a_to_b_primesieve

def primes_from_2_to_n_np(n):
    # http://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python/3035188#3035188
    """ Input n>=6, Returns a array of primes, 2 <= p < n """
    sieve = np.ones(n/3 + (n%6==2), dtype=np.bool)
    sieve[0] = False
    for i in range(int(n**0.5)//3+1):
        if sieve[i]:
            k=3*i+1|1
            sieve[      ((k*k)//3)      ::2*k] = False
            sieve[(k*k+4*k-2*k*(i&1))//3::2*k] = False
    return np.r_[2,3,((3*np.nonzero(sieve)[0]+1)|1)]


def _compute_spins(primes, last_prime, last_spin):
    """Returns array of SPINS given an array of PRIMES (assumed to be 1d)

    SPINS result array is same shape as PRIMES, with values stored in SPINS[1:]; 
    SPINS[0] = 0, use that zero element to thread chunked computations together.
    """
    
    logger.info("compute_spins: generating mod6Values")
    m6val  = primes % 6
    m6_offset_sum     = np.empty_like(primes, dtype=np.int32)
    m6_offset_sum[0]  = m6val[0] + (last_prime % 6)            #  seed value from current val + prev m6val
    m6_offset_sum[1:] = m6val[1:] + m6val[:-1]                 #  cur m6val + prev m6val
    logger.info("compute_spins: done mod6Values")

    logger.info("compute_spins: starting to compute spins")

    z = np.zeros_like(primes, dtype=np.int32)
    z[ m6_offset_sum ==  6] =  1
    z[ m6_offset_sum == 10] = -1
    z[ m6_offset_sum ==  2] = -1
    z[0] *= last_spin
    spin = np.cumprod(z)

    logger.info("compute_spins: done computing spins")

    return spin

def _compute_positions(spin, seed_pos, seed_spin):
    """Given an array of SPINS and two SEED_POSITION and SEED_SPIN values, compute the positions along the prime hex 

    """

    logger.info("compute_positions: starting aux calculations")
    delta = np.zeros_like(spin)
    delta[0]  = spin[0] - seed_spin      # first delta is cur_spin - prev_spin from seed_spin
    delta[1:] = spin[1:] - spin[:-1]     # delta is cur_spin - prev_spin 

    increments = np.copy(spin)           # copy the spin array,
    increments[ delta != 0 ] = 0          # set any non-zero delta to zero in the increment array

    logger.info("compute_positions:\tdone with aux calculations")
    
    logger.info("compute_positions: starting primary calculation")

    # start at seed, cumulative add
    positions = np.copy(increments)
    #increments[0] = seed_pos
    outpositions = (seed_pos + np.cumsum(increments)) % 6
    logger.info("compute_positions:\tdone with primary calculation")
    return outpositions


def _compute_rotations(positions, rot_seed ):

    logger.info("compute_rotations: starting aux calculations")
    delta = np.diff(positions)

    z = np.zeros_like(delta)       # zero array like delta
    z[ delta == -5 ] =  1          # where delta = -5, set increment to 1
    z[ delta ==  5 ] = -1          # where delta is 5, set increment to -1
    logger.info("compute_rotations: done with aux calculations")

    logger.info("compute_rotations: starting primary calculations")
    r = np.cumsum( np.concatenate(([rot_seed],z)))           # cumulative sum the increment values
    logger.info("compute_rotations: done with primary calculations")

    return r

compute_spins = _compute_spins
compute_positions = _compute_positions
compute_rotations = _compute_rotations

def setup_computation_engine(use_cython):
    global compute_spins, compute_positions, compute_rotations

    has_primehex = False

    try:
        import primehexagon
        primehexagon.init()
        has_primehex = True
    except:
        pass

    if has_primehex and use_cython:
        print("using cython primehexagon implementation")
        compute_spins = primehexagon._compute_spins
        compute_positions = primehexagon._compute_positions
        compute_rotations = primehexagon._compute_rotations
    else:
        print("using numpy primehexagon implementation")
    
    
def write_collection_to_file(fobj, data):
    """Write to file object FOBJ each value in DATA collection
    """    
    w = fobj.write
    [ w( str(d) + '\n') for d in data]

def print_text_arrays_to_file( fobj, primes, spin, pos, rot, save_opts=None):
    skip_interval = 1
    if save_opts:
        skip_interval = save_opts.get('skip_interval',1)
    logger.info("start zipping and slicing data together")
    data  = zipper(primes, pos, spin, rot)
    d     = itertools.islice(data, None, None, skip_interval)
    logger.info("\tdone zipping and slicing data")

    write_collection_to_file(fobj, d)

def print_text_arrays( filename, primes, spin, pos, rot, save_opts=None):
    logger.info("saving text results to file {0}".format(filename))
    f = open(filename, "w") 
    print_text_arrays_to_file(f, primes, spin, pos, rot, save_opts=save_opts)
    f.close()
    logger.info("done saving text file")
    
def save_binary_arrays( filename, primes, spin, pos, rot, save_opts=None):

    skip_interval = 1
    if save_opts:
        skip_interval = save_opts.get('skip_interval',1)

    if skip_interval is not None and skip_interval > 1:
        s = slice(None, None, skip_interval)
    else:
        s = slice(None, None, 1)

    do_compress = save_opts and 'compress' in save_opts
        
    (saver,msg) = (np.savez,"uncompressed") if do_compress else (np.savez_compressed, "compressed")
        
    logger.info("start save binary arrays {} to file {} with {}".format(msg, filename, s))
    saver(filename, primes=primes[s], spin=spin[s], pos=pos[s], rot=rot[s])
    logger.info("\tdone save binary arrays")

def compute_chunks_from_files(start_file, start_val, nvals, nchunks = 10, save_opts=None):
    fname = start_file
    val = start_val
    for i in range(nchunks):
        res = compute_chunked_hex_positions(fname, val, nvals, save_opts=save_opts)
        val += nvals
        fname = res[0]
        del res

def compute_chunks_from_scratch(nvals, nchunks = 10, verbose=None, save_opts=None):
    res = compute_hex_positions(nvals, save_opts=save_opts)
    fname = res[0]
    del res
    start_val = nvals
    compute_chunks_from_files( fname, start_val, nvals, nchunks=nchunks-1, save_opts=save_opts)


def test_verbose(nvals=100, verbose=None):

    x = compute_hex_positions(nvals)
    fname,d0 = x[0],x[1:]

    x = compute_chunked_hex_positions('output.npz',nvals=nvals)
    (fname, d1) = (x[0], x[1:])
    newd1 = [ d[1:] for d in d1 ]
    for x in d0:
        print(x)
    for x in newd1:
        print(x)
    save = [d0,newd1]
    val = nvals
    for x in range(1,10):
        x = compute_chunked_hex_positions(fname, val, nvals)
        (fname, d1) = (x[0], x[1:])

        newd1 = [ d[1:] for d in d1 ]
        save.append( newd1 )
        for x in newd1:
            print(x)
        val += nvals
        
    primes =  np.concatenate([ x[0] for x in save])
    spins  =  np.concatenate([ x[1] for x in save])
    poss   =  np.concatenate([ x[2] for x in save])
    rots   =  np.concatenate([ x[3] for x in save])
    print_text_arrays  ( "chunked-output.txt", primes, spins, poss, rots, save_opts=save_opts)
    return save

class HexValues:
    def __init__(self, prime, spin, rot, pos):
        self.prime = prime
        self.spin = spin
        self.rot = rot
        self.pos = pos

def get_last_values_from_file(chunk_file):
    with np.load(chunk_file, mmap_mode='r') as data:
        old_primes = data['primes']
        spin = data['spin']
        rot = data['rot']
        pos = data['pos']
        
    (last_pos, last_spin, last_prime, last_rot) = (pos[-1], spin[-1], old_primes[-1], rot[-1])
    lastValues = HexValues(last_prime, last_spin, last_rot, last_pos)
    
    # we are done with these values, go ahead and close up shop
    del old_primes, spin, rot, pos
    return lastValues

# primes,spin,pos,rot = data['primes'],data['spin'],data['pos'],data['rot']
def compute_chunked_hex_positions(last_chunk_file, start_val, nvals, save_opts = None):

    last_vals = get_last_values_from_file(last_chunk_file)

    end_val   = start_val + nvals

    logger.info("compute_chunked_hex_positions: chunking values from {} to {}".format(start_val,end_val))
    
    working_primes = primes_from_a_to_b(start_val, end_val)
    logger.info("compute_chunked_hex_positions: found {} primes".format(len(working_primes)))

    newspin = compute_spins(working_primes, last_vals.prime, last_vals.spin)
    newpos = compute_positions(newspin, last_vals.pos, last_vals.spin)
    newrot = compute_rotations(newpos, last_vals.rot)

    basename  = "output-{:0>20d}-{:0>20d}".format(start_val, end_val)
    if save_opts is not None and save_opts.get('save_text',False):
        txtoutname = basename + '.txt'
        print_text_arrays( txtoutname, working_primes, newspin, newpos, newrot, save_opts=save_opts)
    outname = basename + '.npz'
    save_binary_arrays( outname, working_primes, newspin, newpos, newrot,save_opts=save_opts)
                      
    
    return (outname, working_primes, newspin, newpos, newrot)
    
def compute_hex_positions(end_num, save_opts=None):

    raw_primes = primes_from_a_to_b(1, end_num)

    # 2 and 3 are special, don't use them (they are the first two values, slice them out)
    working_primes = raw_primes[2:]

    spin = compute_spins(working_primes, 1, 1) # 1 as last prime is a special case to start off the computation
    pos = compute_positions(spin, 1, 1)
    rot = compute_rotations(pos, 0)

    basename = "output-{:0>20d}-{:0>20d}".format(0, end_num)
    if save_opts is not None and save_opts.get('save_text',False):
        txtname = basename + '.txt'
        print_text_arrays  ( txtname, working_primes, spin,  pos, rot, save_opts=save_opts)
    npzname = basename + '.npz'
    save_binary_arrays( npzname, working_primes, spin, pos, rot, save_opts=save_opts)
    
    return (npzname, raw_primes, working_primes, spin, pos, rot)

def countify2(ar):
    # http://stackoverflow.com/questions/4260645/how-to-get-running-counts-for-numpy-array-values
    # not working,
    ar2 = np.ravel(ar)
    ar3 = np.empty(ar2.shape, dtype=np.int32)
    uniques = np.unique(ar2)
    myarange = np.arange(ar2.shape[0])
    for u in uniques:
        ar3[ar2 == u] = myarange
    return ar3

def print_npz_vals(infile, slices):
    with np.load(infile, mmap_mode='r') as data:
        primes = data['primes']
        spin = data['spin']
        rot = data['rot']
        pos = data['pos']
    for s in slices:
        data  = zipper(primes[s], pos[s], spin[s], rot[s])
        write_collection_to_file(sys.stdout, data)
    
        
def get_slice_obj(slicearg):
    slicearg = re.sub(',','',slicearg)
    svals = [ int(n) if n else None for n in slicearg.split(':') ]
    svals = tuple(svals)
    s = slice(*svals)
    return s

def main(argv = None):
    import sys, argparse, os, re
    if argv is None:
        argv = sys.argv
            
    parser = argparse.ArgumentParser(description='Prime Spin Hexagons')
    parser.add_argument('--infile', help='Input file to start processing chunks',required=False)
    parser.add_argument('--viewvalues', help='Values to view in the file as a python slice, e.g. 1:100:', required=False, action='append',)
    parser.add_argument("-c", "--compress", help="flag to indicate whether or not to compress output files", action="store_true", default=True)
    parser.add_argument('--logfile', help='Save messagse to this log file',required=False)
    parser.add_argument('--verbose', help='Print messages to the terminal',required=False, default=1)
    parser.add_argument('--save-text', help='Flag to save text files',required=False, action="store_true", default=False)
    parser.add_argument('--save-binary', help='Flag to save binary files',required=False,action="store_true", default=True)
    parser.add_argument('--use-cython', help='Flag to use cython implementation if available',dest="use_cython",action="store_true")
    parser.add_argument('--no-use-cython', help='Flag to not use cython implementation',dest="use_cython",action="store_false")
    parser.add_argument('--nvalues', help="number of values to process in a chunk",
                        default=10**9, type=long)
    parser.add_argument('--skip', help="number of values to skip in printing output",
                        default=1, type=long)
    parser.add_argument("--chunks", help="number of chunks to process", default=10,  type=int)
    args = parser.parse_args()
    
    global logger
    
    if args.logfile:
        fh = logging.FileHandler(args.logfile)
        fh.setLevel(logging.INFO)
        fh.setFormatter(log_formatter)
        # add the handlers to the logger
        logger.addHandler(fh)

    if args.verbose:
        # create console handler and set level to debug
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        # add formatter to ch
        ch.setFormatter(log_formatter)
        
        # add ch to logger
        logger.addHandler(ch)

    setup_computation_engine(args.use_cython)
        
    save_opts = { 'compress' : args.compress,
                  'skip_interval' : args.skip,
                  'save_text' : args.save_text }
    print "Save options:" + str(save_opts)
    if args.infile is None:
        compute_chunks_from_scratch(args.nvalues, nchunks=args.chunks, save_opts=save_opts)

    else:
        # we have a starting file, figure some stuff out
        if not os.path.isfile(args.infile):
            sys.stderr.write("{} file does not exist".format(args.infile))
            sys.exit(1)
        startval, endval = [long(x) for x in re.findall('\d+',args.infile)]
        nvalues = endval - startval

        if args.viewvalues:
            slices = [ get_slice_obj(s) for s in args.viewvalues]
            print_npz_vals(args.infile, slices)
        else:
            compute_chunks_from_files(args.infile, endval, nvalues, nchunks=args.chunks, save_opts=save_opts)
        

if __name__ == "__main__":
    import sys
    main(sys.argv)



