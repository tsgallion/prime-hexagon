# prime-hexagon

Prime numbers and their positions in the prime hexagon.

This code computes the positons in the central triangle of the prime hexagon. 

What is the prime hexagon, you may ask? TODO...

## Running

Compute the first set of output values

    python ./primespin.py --nvalues=10000000 --chunks 100

Resume a computation given a saved file:

    python ./primespin.py --startfile ./output-00000000000000107000-00000000000000108000.npz --chunks 100


Help with --help

```bash
python ./primespin.py --help
usage: primespin.py [-h] [--startfile STARTFILE] [--startvalue STARTVALUE]
                    [-c] [--logfile LOGFILE] [--nvalues NVALUES]
		                        [--chunks CHUNKS]

Prime Spin Hexagons

optional arguments:
  -h, --help            show this help message and exit
  --startfile STARTFILE
                  Input file to start processing chunks
  --startvalue STARTVALUE
                  Starting value for resumed chunk computations
  -c, --compress
  --logfile LOGFILE     File where to save the log
  --nvalues NVALUES     number of values to process in a chunk
  --chunks CHUNKS       number of chunks to process
```


output-00000000000000009000-00000000000000010000.npz

##  Computing a prime number's position on the prime hexagon

Initial prime within the triangle is at position 2 and is prime 5.  Numbers at 2,4,0
are dependent on status of those at 1,3, and 5 for spin determination for numbers
immediately following.  (Because they are n type primes and defined as residing
in the last hexagon.  The following number, n+1, spins with the n number, unless
it is followed by another prime, then spins with the second prime.
A sketch of the inner triangle:


           <- 3 ->
           
           2      4
           
       <-1    0     5->
         |          |
         v          v
         
