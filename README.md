# prime-hexagon

Prime numbers and their positions in the prime hexagon.

This code computes the positons in the central triangle of the prime hexagon. 

What is the prime hexagon, you may ask? TODO...

## Running

Compute the first set of output values

    python ./primespin.py --nvalues=10000000 --chunks 100

Resume a computation given a saved file:

    python ./primespin.py --startfile ./output-00000000000000107000-00000000000000108000.npz --chunks 100


##  Computing a prime number's position on the prime hexagon

TODO

## Command Line Help 

Run with the --help argument to see the help file.

```bash
python primespin.py --help
usage: primespin.py [-h] [--startfile STARTFILE] [--startvalue STARTVALUE]
                    [-c] [--logfile LOGFILE] [--verbose VERBOSE]
                    [--nvalues NVALUES] [--chunks CHUNKS]

Prime Spin Hexagons

optional arguments:
  -h, --help            show this help message and exit
  --startfile STARTFILE
                        Input file to start processing chunks
  --startvalue STARTVALUE
                        Starting value for resumed chunk computations
  -c, --compress
  --logfile LOGFILE     Save messagse to this log file
  --verbose VERBOSE     Print messages to the terminal
  --nvalues NVALUES     number of values to process in a chunk
  --chunks CHUNKS       number of chunks to process
```

Here are a few sample command lines:

Generate 1000 files of 1B values each

    primespin.py --nvalues=1000000000 --chunks=1000 

or

     python3 primespin.py --nvalues=$(echo 10,000,000,000 | tr -d ,) --chunks=10000 --skip=1000000 --logfile run.log

Same thing but write out only every 1Mth value:

    primespin.py --nvalues=10000000000 --chunks=10 --skip=1000000 

Now resume a computation starting at the specified file:

    primespin.py --chunks=10 --skip=1000000 --infile output-00000000009000000000-00000000010000000000.npz

Print some values from the binary .npz files

    primespin.py --viewvalues=:: --infile output-00000000009000000000-00000000010000000000.npz


Note viewvalues is evaluated as an [extended Python
slice](https://docs.python.org/2.3/whatsnew/section-slices.html). That is, the values are

   start:stop:step

Also any of those values can actually be negative values to indicate starting from the other end in the case of start and stop. If you make step negative it counts backward.

So let's print the last value in this file:

    primespin.py --viewvalues=-1:: --infile output-00000000009000000000-00000000010000000000.npz
