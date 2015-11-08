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

logger = None
log_formatter = None

def setup_loggers():

    global logger, log_formatter
    if logger is not None:
        return
    
    # create logger
    logger = logging.getLogger('prime_hex')
    logger.setLevel(logging.INFO)

    # create formatter
    log_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')



setup_loggers()

def primes_from_a_to_b_primesieve(a,b):
    logger.info("starting generating primes from {} to {}".format(a,b))
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


def compute_spins(primes, last_spin):
    """Returns array of SPINS given an array of PRIMES (assumed to be 1d)

    SPINS result array is same shape as PRIMES, with values stored in SPINS[1:]; 
    SPINS[0] = 0, use that zero element to thread chunked computations together.
    """
    
    logger.info("compute_spins: generating mod6Values")
    m6val  = primes % 6
    m6_offset_sum     = np.zeros_like(primes)
    m6_offset_sum[1:] = m6val[:-1] + m6val[1:]
    logger.info("compute_spins: done mod6Values")

    logger.info("compute_spins: starting to compute spins")

    z = np.zeros_like(primes)
    z[ m6_offset_sum ==  6] =  1
    z[ m6_offset_sum == 10] = -1
    z[ m6_offset_sum ==  2] = -1
    z[0] = last_spin
    spin = np.cumprod(z)

    logger.info("compute_spins: done computing spins")

    return spin

def compute_positions(spin, seed_pos, seed_spin):
    """Given an array of SPINS and two SEED_POSITION and SEED_SPIN values, compute the positions along the prime hex 

    """

    logger.info("compute_positions: starting aux calculations")
    delta = np.zeros_like(spin)
    delta[0]  = spin[0] - seed_spin      # first delta is seed_spin from previous chunk to this first spin
    delta[1:] = spin[1:] - spin[0:-1]   # compute rest of deltas from input spins array

    increments = np.copy(spin)           # copy the spin array,
    increments[ delta != 0 ] = 0          # set any non-zero delta to zero in the increment array

    logger.info("compute_positions:\tdone with aux calculations")
    
    logger.info("compute_positions: starting primary calculation")

    # start at seed, cumulative add
    positions = np.copy(increments)
    positions[0] += seed_pos
    outpositions = np.cumsum(positions) % 6
    logger.info("compute_positions:\tdone with primary calculation")
    return outpositions


def compute_rotations(positions, rot_seed ):

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


def print_outputs(filename, data, skip=None):
    """Fold four major outputs into same result: 
    prime, position, spin, rotation
    """    
     
    f = open(filename, "w") 
 
    if skip is not None and skip > 1:
        logger.info("skipping values - using every {0}".format(skip))
        data = data[::skip]

    for out in data:
        s = str(out) + '\n'
        f.write(s)
    f.close()

def save_text_arrays( filename, primes, spin, pos, rot, skip_interval=1):
    logger.info("start zipping and slicing data together")
    data              = zip(primes, pos, spin, rot)
    d = itertools.islice(data, 0, len(primes), skip_interval)
    logger.info("\tdone zipping and slicing data")

    logger.info("saving text results to file {0}".format(filename))
    print_outputs(filename, d)
    logger.info("done saving text file")
    
def save_binary_arrays( filename, primes, spin, pos, rot, skip_interval=None, do_compress=None):

    if skip_interval is not None and skip_interval > 1:
        s = slice(None, None, skip_interval)
    else:
        s = slice(None, None, 1)

    saver = np.savez
    msg = "uncompressed"
    if do_compress:
        msg = "compressed"
        saver = np.savez_compressed
    logger.info("start save binary arrays {} to file {} with slice {}".format(msg, filename, s))
    saver(filename, primes=primes[s], spin=spin[s], pos=pos[s], rot=rot[s])
    logger.info("\tdone save binary arrays")

def blow_chunky_chunks(start_file, start_val, nvals, nchunks = 10, verbose=None, do_compress=None):
    fname = start_file
    val = start_val
    for i in range(nchunks-1):
        res = compute_chunked_hex_positions(fname, val, nvals, do_compress=do_compress)
        val += nvals
        fname = res[0]
        del res
    
def blow_chunks(nvals, nchunks = 10, verbose=None, do_compress=None):
    res = compute_hex_positions(nvals, do_compress=do_compress)
    fname = res[0]
    del res
    start_val = nvals
    blow_chunky_chunks( fname, start_val, nvals, nchunks=nchunks, do_compress=do_compress)


def test_verbose(nvals=100, verbose=None):

    x = compute_hex_positions(nvals)
    fname,d0 = x[0],x[1:]

    x = compute_chunked_hex_positions('output.npz',nvals=nvals)
    (fname, d1) = (x[0], x[1:])
    newd1 = [ d[1:] for d in d1 ]
    for x in d0:
        print x
    for x in newd1:
        print x
    save = [d0,newd1]
    val = nvals
    for x in range(1,10):
        x = compute_chunked_hex_positions(fname, val, nvals)
        (fname, d1) = (x[0], x[1:])

        newd1 = [ d[1:] for d in d1 ]
        save.append( newd1 )
        for x in newd1:
            print x
        val += nvals
        
    primes =  np.concatenate([ x[0] for x in save])
    spins  =  np.concatenate([ x[1] for x in save])
    poss   =  np.concatenate([ x[2] for x in save])
    rots   =  np.concatenate([ x[3] for x in save])
    save_text_arrays  ( "chunked-output.txt", primes, spins, poss, rots)
    return save
        
# primes,spin,pos,rot = data['primes'],data['spin'],data['pos'],data['rot']
def compute_chunked_hex_positions(last_chunk_file, start_val, nvals, skip_interval=1, do_compress=None):

    with np.load(last_chunk_file, mmap_mode='r') as data:
        old_primes = data['primes']
        spin = data['spin']
        rot = data['rot']
        pos = data['pos']

    (last_pos, last_spin, last_prime, last_rot) = (pos[-1], spin[-1], old_primes[-1], rot[-1])
    # we are done with these values, go ahead and close up shop
    del old_primes, spin, rot, pos

    end_val   = start_val + nvals

    logger.info("compute_chunked_hex_positions: chunking values from {} to {}".format(start_val,end_val))
    
    working_primes = primes_from_a_to_b(start_val, end_val)
    logger.info("compute_chunked_hex_positions: found {} primes".format(len(working_primes)))

    newspin = compute_spins(working_primes, last_spin)
    newspin[0] = last_spin     ## special condition of stiching start value to right 

    temppos = compute_positions(newspin[1:], last_pos, last_spin)
    newpos = np.concatenate( ([last_pos], temppos) )

    newrot = compute_rotations(newpos, last_rot)

    #save_text_arrays  ( "output.txt", working_primes, spin, pos, rot, skip_interval=skip_interval) 
    outname = "output-{:0>20d}-{:0>20d}.npz".format(start_val, end_val)
    save_binary_arrays( outname, working_primes, newspin, newpos, newrot,
                        skip_interval=skip_interval, do_compress=do_compress) 
    
    return (outname, working_primes, newspin, newpos, newrot)
    
def compute_hex_positions(end_num, skip_interval=1, do_compress=None):

    raw_primes = primes_from_a_to_b(1, end_num)

    # 2 and 3 are special, don't use them (they are the first two values, slice them out)
    working_primes = raw_primes[2:]

    spin = compute_spins(working_primes, 1)
    pos = compute_positions(spin, 1, 1)
    rot = compute_rotations(pos, 0)

    #save_text_arrays  ( "output.txt", working_primes, spin,  pos, rot, skip_interval=skip_interval) 
    outname = "output-{:0>20d}-{:0>20d}.npz".format(0, end_num)
    save_binary_arrays( outname, working_primes, spin, pos, rot,
                        skip_interval=skip_interval, do_compress=do_compress) 
    
    return (outname, raw_primes, working_primes, spin, pos, rot)

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
    
                                    
def main(argv = None):
    import sys, argparse, os, re
    if argv is None:
        argv = sys.argv
            
    parser = argparse.ArgumentParser(description='Prime Spin Hexagons')
    parser.add_argument('--startfile', help='Input file to start processing chunks',required=False)
    parser.add_argument('--startvalue', help='Starting value for resumed chunk computations',
                        required=False, type=long)
    parser.add_argument("-c", "--compress", action="store_true", default=False)
    parser.add_argument('--logfile', help='Save messagse to this log file',required=False)
    parser.add_argument('--verbose', help='Print messages to the terminal',required=False)
    parser.add_argument('--nvalues', help="number of values to process in a chunk",
                        default=10**9, type=long)
    parser.add_argument("--chunks", help="number of chunks to process", default=10,  type=int)
    args = parser.parse_args()
    
    if args.logfile:
        fh = logging.FileHandler(args.logfile)
        fh.setLevel(logging.INFO)
        fh.setFormatter(log_formatter)
        # add the handlers to the logger
        logger.addHandler(fh)

    if args.verbose:
        global logger
        # create console handler and set level to debug
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        # add formatter to ch
        ch.setFormatter(log_formatter)
        
        # add ch to logger
        logger.addHandler(ch)
        
    if args.startfile is None:
        print args
        blow_chunks(args.nvalues, nchunks=args.chunks, do_compress=args.compress)

    else:
        # we have a starting file, figure some stuff out
        if not os.path.isfile(args.startfile):
            sys.stderr.write("{} file does not exist".format(args.startfile))
            sys.exit(1)
        startval, endval = [long(x) for x in re.findall('\d+',args.startfile)]
        nvalues = endval - startval
        
        blow_chunky_chunks(args.startfile, startval, nvalues, nchunks=args.chunks)
        
    #compute_hex_positions(end_num, 1)

if __name__ == "__main__":
    import sys
    main(sys.argv)



