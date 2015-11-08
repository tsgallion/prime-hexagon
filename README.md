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
primespin.py:308: SyntaxWarning: name 'logger' is used prior to global declaration
  global logger
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
