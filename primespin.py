"""  The Number Coil  """
#  List a prime; List mod6 of prime; Compare with another set (offset of mod6)
#  to see spin change, determine cell of each prime.
#
#

"""  This code plots the positons in the central triangle of Tad Gallion's hexagon.
Initial prime within the triangle is at position 2 and is prime 5.  Numbers at 2,4,0
are dependent on status of those at 1,3, and 5 for spin determination for numbers
immediately following.  (Because they are n type primes and defined as residing
in the last hexagon.  The following number, n+1, spins with the n number, unless
it is followed by another prime, then spins with the second prime.
A sketch of the inner triangle:
.
           <- 3 ->
           
           2      4
           
       <-1    0     5->
         |          |
         v          v

"""

import itertools
import numpy as np
import primesieve as ps

def primes_from_a_to_b(a,b):
    return ps.generate_primes(a,b)
    
def primes_from_2_to(n):
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


def generate_primes_from_a_to_b_old(start,end):
    """Generates the list of Primes //
    else Converts primes by types n or n+2 types (mod 1 or 5) """
    p = []
    m6vals = []
    for i in range(start, end):
        for x in range(start, i):
            if( i%x == 0):
                break
        else:
            p.append(i)
            if(i%6==5):
                m6vals.append(5)
            else:
                m6vals.append(1)        

    """Get rid of artifacts cause by 2 and 3"""
    p.remove(2)
    p.remove(3)
    m6vals.remove(1)
    m6vals.remove(1)

    return (p,m6vals)



def compute_spins(primes):
    """Returns array of SPINS given an array of PRIMES (assumed to be 1d)

    SPINS result array is same shape as PRIMES, with values stored in SPINS[1:]; 
    SPINS[0] = 0, use that zero element to thread chunked computations together.
    """
    
    print("generating mod6Values...")
    m6val  = primes % 6
    print("generating offsets...")
    m6_offset_sum = m6val[:-1] + m6val[1:]
    z = np.copy(m6_offset_sum)
    z[ m6_offset_sum == 6]  = 1
    z[ m6_offset_sum == 10] = -1
    z[ m6_offset_sum == 2]  = -1
    spins = np.cumprod(z)

    print("\tgot spins")

    # TODO: assert that values all values are 1 or -1, i.e. there are no zeros
    assert np.count_nonzero(primes) == np.count_nonzero(spins) + 1

    out = np.zeros_like(primes)
    out[1:] = spins
    return out

def compute_positions(seed_pos, seed_spin, spins):
    delta = np.copy(spins)
    delta[0] = spins[0] - seed_spin     # first delta is seed_spin from previous chunk to this first spin
    delta[1:] = spins[1:] - spins[0:-1]  # compute rest of deltas from input spins array
    increments = np.copy(spins)
    increments[ delta != 0 ] = 0

    # start at seed, cumulative add
    positions = np.copy(increments)
    positions[0] = increments[0] +  seed_pos
    outpositions = np.cumsum(positions)
    return outpositions


def generate_spins_iterative(end, m6vals):
    """Gernerates spin //Gernerates Positions // mod6, to keep it base 6"""
    
    # Creates an offset list (offset by 1 item) to compare to the mod6
    # list - so if an n prime leads to another n prime, you see a turn,
    # etc.  That is - if 23 is mod 5 and 29 is mod 5
    # the total 10 says, turn around
    offSet = [1] + m6vals
    compare = [sum(x) for x in zip(m6vals, offSet)]
    
    spins = []
    positions = []
    pos = 1
    spin = 1
    for i in compare:
        if (i==10) or (i == 2):
            spin = spin* -1
        else:
            pos = (pos + spin)%6
        spins.append(spin)
        positions.append(pos)

    return (spins,positions) 



def generate_rotations(end, positions):
    """Determine Overall Rotation (from the centerpoint of the
    hexagon.  Positions:  another offset to Rotation to get central rotation   """

    rotOffset = [0] + positions
    posDelta = [i - j for i, j in zip(positions, rotOffset)]

    rotations = []
    rotCount = 0
    for i in posDelta:
        if (i==-5):
            rotCount=rotCount+1
        elif (i==5):
            rotCount=rotCount-1
        rotations.append(rotCount)
    
    return rotations



def print_outputs(filename, data, skip=None):
    """Fold four major outputs into same result: prime, position, spin
,rotation"""    
     
    f = open(filename, "w") 
 
    if skip is not None:
        print("skipping values - using every {0}".format(skip))
        data = data[::skip]
    #skip = datafullOutX4[::1]

    for out in data:
        s = str(out) + '\n'
        f.write(s)
    f.close()


def compute_hex_positions(end_num):


    print("generating primes...")
    #primes = primes_from_2_to(end_num)

    raw_primes = primes_from_a_to_b(1,1000)

    
    # 2 and 3 are special, don't use them (they are the first two values, slice them out)
    working_primes = np.array(raw_primes[2:])
    print("\tgot primes")

    spins = compute_spins(working_primes)
    spins[0] = 1     ## special condition of stiching start value to right 
    print("\tgot spins")

    pos = compute_positions(1, 1, spins)
    
    print("generating rotations...")
    rot            = generate_rotations(end_num, pos)
    print("\tgot rotations")

    print("zipping all data together")
    data              = zip(working_primes, pos, spins, rot)
    print("\tzipped data.")
    print("slicing file")
    d = itertools.islice(data,0,end_num,1)

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

if __name__ == '__main__':
    import sys

    if (len(sys.argv) > 1): 
        end_num = int(sys.argv[1])
    else:
        end_num = 100000

    print("using end: {0}".format(end_num))
    compute_hex_positions(end_num)




