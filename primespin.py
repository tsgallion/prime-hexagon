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

logger = logging.getLogger(__name__)

def primes_from_a_to_b_primesieve(a,b):
    logger.info("starting generating primes: {} to {}".format(a,b))
    out = ps.generate_primes(a,b)
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


def compute_spins(primes):
    """Returns array of SPINS given an array of PRIMES (assumed to be 1d)

    SPINS result array is same shape as PRIMES, with values stored in SPINS[1:]; 
    SPINS[0] = 0, use that zero element to thread chunked computations together.
    """
    
    logger.info("compute_spins: generating mod6Values")
    m6val  = primes % 6
    m6_offset_sum = m6val[:-1] + m6val[1:]
    logger.info("compute_spins: done mod6Values")

    logger.info("compute_spins: starting to compute spins")

    z = np.copy(m6_offset_sum)
    z[ m6_offset_sum == 6]  = 1
    z[ m6_offset_sum == 10] = -1
    z[ m6_offset_sum == 2]  = -1
    spins = np.cumprod(z)

    logger.info("compute_spins: done computing spins")

    # assert that values all values are 1 or -1, i.e. there are no zeros
    assert np.count_nonzero(primes) == np.count_nonzero(spins) + 1

    out = np.append( [0], spins)
    #out = np.zeros_like(primes)
    #out[1:] = spins
    return out

def compute_positions(spins, seed_pos, seed_spin):
    """Given an array of SPINS and two SEED_POSITION and SEED_SPIN values, compute the positions along the prime hex 

    """

    logger.info("compute_positions: starting aux calculations")
    delta = np.zeros_like(spins)
    delta[0]  = spins[0] - seed_spin      # first delta is seed_spin from previous chunk to this first spin
    delta[1:] = spins[1:] - spins[0:-1]   # compute rest of deltas from input spins array

    increments = np.copy(spins)           # copy the spin array,
    increments[ delta != 0 ] = 0          # set any non-zero delta to zero in the increment array

    logger.info("compute_positions: done with aux calculations")
    
    logger.info("compute_spins: starting primary calculation")

    # start at seed, cumulative add
    positions = np.copy(increments)
    positions[0] = increments[0] +  seed_pos
    outpositions = np.cumsum(positions) % 6
    logger.info("compute_spins: done with primary calculation")
    return outpositions


def compute_rotations(positions, rot_seed):

    logger.info("compute_rotations: starting aux calculations")
    delta = np.diff(positions)

    z = np.zeros_like(delta)       # zero array like delta
    z[ delta == -5 ] = 1           # where delta = -5, set incrment to 1
    z[ delta == 5 ]  = -1          # where delta is 5, set increment to -1
    logger.info("compute_rotations: done with aux calculations")

    logger.info("compute_rotations: starting primary calculations")
    r = np.zeros_like(positions)   # construct output array shaped like input array
    r[0] = rot_seed                # set seed start value
    r[1:] = np.cumsum(z)           # cumulative sum the increment values
    logger.info("compute_rotations: done with primary calculations")

    return r


def print_outputs(filename, data, skip=None):
    """Fold four major outputs into same result: 
    prime, position, spin, rotation
    """    
     
    f = open(filename, "w") 
 
    if skip is not None:
        logger.info("skipping values - using every {0}".format(skip))
        data = data[::skip]

    for out in data:
        s = str(out) + '\n'
        f.write(s)
    f.close()


def compute_hex_positions(end_num):


    print("generating primes from 1 to {}...".format(end_num) )

    raw_primes = primes_from_a_to_b(1, end_num)

    # 2 and 3 are special, don't use them (they are the first two values, slice them out)
    working_primes = np.array(raw_primes[2:])

    spins = compute_spins(working_primes)
    spins[0] = 1     ## special condition of stiching start value to right 

    pos = compute_positions(spins, 1, 1)
    rot = compute_rotations(pos, 0)
    

    print("zipping all data together")
    data              = zip(working_primes, pos, spins, rot)
    print("\tzipped data.")
    print("slicing file")
    d = itertools.islice(data,0,end_num,100)

    f = "test.txt"
    
    """ print("convert via skip")
     #Prints every nth value set

    skip = data[1:end_num:100000]
     #print (skip)
    value = str(skip)
    print("\tskipped.")"""
    
    print("saving results to file {0}".format(f))
    print_outputs(f, d)
    print("\tprinted data")
    return (raw_primes, working_primes, pos, spins, rot)

if __name__ == '__main__':
    import sys

    if (len(sys.argv) > 1): 
        end_num = int(sys.argv[1])
    else:
        end_num = 100000

    print("using end: {0}".format(end_num))
    compute_hex_positions(end_num)




