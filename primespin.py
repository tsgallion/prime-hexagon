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
#YESTERDAY"S WORKED BUT THE TUPLE TO LIST STEP CAUSING TIME PROBLEMS AND ALSO SEEMS TO
#STOP AT 1B.  PYTHON 2 to 3 ISSUE _ WAS LIST RETURN FROM ZIP, NOW TUPLE
#LOW SWAP ISSUES AT 50Billion AND EVEN 25B AND 15B.  Despite 32GB, needs big swap.


# Generate a list of the primes below 40

import itertools
import numpy as np
import primesieve as ps

#input("To what number would you like to run the coil to?")

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






def generate_spins(end, m6vals):
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

def compute_hex_positions_np(end_num):
    """
    Try to do the same computations using array operations

    Status: incomplete, not known to work
    """

    #primes = primes_from_2_to(end_num)
    primes = primes_from_a_to_b(end_num)
    
    primes = primes[2:]

    m6vals = np.ones_like(primes)
    m6vals[primes % 6 == 5] = 5

    # 1 or -1s for each prime
    #spin_compares = np.insert(m6vals[:-1] + m6vals[1:],0,0)
    spin_compares = m6vals[:-1] + m6vals[1:]
    spin_signs = np.ones_like(primes)
    spin_signs[spin_compares == 10] = -1
    spin_signs[spin_compares == 2] = -1

    spin_sums = np.ones_like(primes)
    spin_sums[spin_compares == 10] = 0
    spin_sums[spin_compares == 2] = 0

    turns = np.insert(np.diff(spin_signs),0,0) // 2

    positions = np.ones_like(primes) + 1
    #positions = (positions + signs) % 6

def compute_hex_positions(end_num):

    primes_from_a_to_b(1,1000)

    print("generating primes...")
    primes = primes_from_2_to(end_num)
    primes = primes[2:]
    print("\tgot primes")

    print("generating mod6Values...")
    m6vals = np.ones_like(primes)
    m6vals[primes % 6 == 5] = 5
    print("\tgot mod6Values")

    # old way
    #(primes, mod6Val) = generate_primes_from_a_to_b_old(2,end_num)
    #(spinLib, posLib) = generate_spins(end_num, mod6Val)
    print("getting spins")
    (spinLib, posLib) = generate_spins(end_num, m6vals.tolist())
    print("\tgot spins")

    print("generating rotations...")
    rotLib            = generate_rotations(end_num, posLib)
    print("\tgot rotations")

    print("zipping all data together")
    data              = zip(primes, posLib,spinLib,rotLib)
    print("\tzipped data.")
    print("slicing file")
    d = itertools.islice(data,1,end_num,100)

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




