"""  The Number Coil  """
#  List a prime; List mod6 of prime; Compare with another set (offset of mod6)
#  to see spin change, determine cell of each prime.
#
#

"""  This code plots the positons in the central triangle of the prime hexagon.

See accompanying notes and images for what these mean.

"""

import sys, os, glob
import logging
import re
import itertools

import numpy as np
import primesieve as ps

# python 2/3 compatibility stuff
try:
    l = long(1)
except:
    long = int

try:
    # python2 version of zip iterator is in itertools
    zipper = itertools.izip
except:
    # new python3 version is no longer in itertools
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
        print("fell back to classic primesieve prime generation")
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

# these functions are also implemented in cython versions
# cython versions are much faster.
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
        logger.info("using cython primehexagon implementation")
        compute_spins = primehexagon._compute_spins
        compute_positions = primehexagon._compute_positions
        compute_rotations = primehexagon._compute_rotations
    else:
        logger.info("using numpy primehexagon implementation")
    
    
def write_collection_to_file(fobj, data):
    """Write to file object FOBJ each value in DATA collection
    """
    w = fobj.write
    [ w( str(d) + '\n') for d in data]

def print_text_arrays_to_file( fobj, data, save_opts=None):
    skip_interval = 1
    if save_opts:
        skip_interval = save_opts.get('skip_interval',1)
    logger.info("start zipping and slicing data together")
    data  = zipper(data.prime, data.pos, data.spin, data.rot)
    d     = itertools.islice(data, None, None, skip_interval)
    logger.info("\tdone zipping and slicing data")

    write_collection_to_file(fobj, d)

def print_text_arrays( filename, data, save_opts):
    logger.info("saving text results to file {0}".format(filename))
    f = open(filename, "w") 
    print_text_arrays_to_file(f, data, save_opts=save_opts)
    f.close()
    logger.info("done saving text file")
    
def save_binary_arrays( filename, data, save_opts):

    skip_interval = save_opts.get('skip_interval', 1)
    s = slice(None, None, skip_interval)

    do_compress = save_opts.get('compress', 1)
        
    (saver,msg) = (np.savez,"uncompressed") if do_compress else (np.savez_compressed, "compressed")
        
    logger.info("start save binary arrays {} to file {} with {}".format(msg, filename, s))
    saver(filename, primes=data.prime[s], spin=data.spin[s], pos=data.pos[s], rot=data.rot[s])
    logger.info("\tdone save binary arrays")


def compute_hex_positions(start, end, boundaryVals, save_opts):

    primes = primes_from_a_to_b(start, end)

    # 2 and 3 are special, don't use them (they are the first two values, slice them out)
    if start < 2:
        primes = primes[2:]
    
    # 1 as last prime and last_spin and last_pos arespecial cases to start off the computation
    spin = compute_spins(primes, boundaryVals.prime, boundaryVals.spin)
    pos  = compute_positions(spin, boundaryVals.pos, boundaryVals.spin)
    rot  = compute_rotations(pos, boundaryVals.rot)

    data = HexValues(primes, pos, spin, rot)
    newHexFile = HexDataFile.save_new_file(start, end, data, save_opts)
    return newHexFile

def _require_infile(infile):
    # we have a starting file, figure some stuff out
    if not os.path.isfile(infile):
        raise IOError("required input file {} does not exist".format(infile))

_VERBOSE = 1
def dprint(level, *args):
    if _VERBOSE >= level:
        logger.info(*args)

_HEX_FILENAME_RE = re.compile(r"-(?P<start>\d+)-(?P<end>\d+)(-(?P<skip>\d+))?.npz")
_HEX_FILENAME_GLOB = '-[0-9]*-[0-9]*.npz'
        
class HexDataAssets:
    def __init__(self, opts):
        self.basedir = opts.get('basedir','output')
        self.basename = opts.get('basename','output')
        self._files = None

    @classmethod
    def get_hex_assets(klass, opts = None):
        if opts is None:
            basedir = 'data'
            basename = 'output'
        else:
            basedir = opts.dir or 'data'
            basename = opts.basename or 'output'
        h = klass(basedir, basename)
        return h

    def _clear_file_cache(self):
        self._files = None
    
    def list_files(self):
        # check cache
        if self._files:
            return files
        base_pattern = os.path.join(self.basedir, self.basename)
        full_pattern = base_pattern + _HEX_FILENAME_GLOB
        dprint(1, "looking for files matching pattern {}".format(full_pattern))
        filenames = sorted(glob.glob(full_pattern))
        dprint(1, "found {} files".format(len(filenames)))
        self._files = filenames
        return filenames
    
    def _find_files(self, vals, hex_predicate):
        """Return a SET of any files that contain any of the values specified"""
        s = set()
        for f in self.list_files():
            h = HexDataFile(f)
            dprint(2,"inspecting file {} for values {}".format(h, vals))
            if hex_predicate(h, vals):
                dprint(2,"\tsaving file for values {}".format(h, vals))
                s.add(h)
        return s

    def find_values_in_npz(self, vals):
        hexfiles = self.find_files_with_values(vals)
        if not hexfiles:
            print("no files found containing the value {} in directory {}".format(vals, basedir))
            return

        for hfile in hexfiles:
            # stupid check to see if any of the values are in this file.
            indexes = hfile.values_around(vals)
            if indexes:
                print("file {}".format(hfile.filename))
            else:
                next
            # now figure out which values are in this file...
            for v in vals:
                indexes = hfile.values_around((v,))
                if indexes:
                    dprint(2, "found values in file {}".format(hfile))
                for idx in indexes:
                    # TODO: get it if primes are in the previous or next file
                    #  i.e. if idx-2 < 0, get primes from previous files
                    #       if idx+2 > length, get primes from next file
                    low = max(0, idx-2)
                    hi  = min( hfile.length(), idx+2)
                    dprint(2,"for val={} found index={} with low={}, hi={}".format(v, idx, low, hi))

                    data = hfile.get_sliced_zipped_arrays_iterator(slice(low,idx))
                    write_collection_to_file(sys.stdout, data)

                    print('----> {} <----'.format(v))

                    zipped_arrays = hfile.get_zipped_arrays_iterator()
                    data = hfile.get_sliced_zipped_arrays_iterator(slice(idx, hi))
                    write_collection_to_file(sys.stdout, data)

                    print("")
    
    def find_files_with_values(self, vals):
        """Return a SET of any files that contain any of the values specified"""
        return self._find_files(vals, HexDataFile.contains_any_values)
        s = set()
        for f in self.list_files():
            h = HexDataFile(f)
            dprint(2,"inspecting file {} for values {}".format(h, vals))
            if h.contains_any_values(vals):
                dprint(2,"\tsaving file for values {}".format(h, vals))
                s.add(h)
        return s
    
    def find_file_with_value(self, val):
        """Return a SET of any files that contain the value VAL"""
        return self._find_files(vals, HexDataFile.contains)
        s = set()
        for f in self.list_files():
            h = HexDataFile(f)
            dprint(2,"inspecting file {} for val {}".format(h, val))
            if h.contains(val):
                dprint(2,"\tsaving file {} for value {}".format(h, val))
                s.add(h)
        return s

    def compute_next_chunks(self, startFile, nvalues, nchunks, save_opts):
        hexFile = startFile
        start_val = startFile.end
        for i in range(nchunks):
            end_val   = start_val + nvalues
            linkedValues = hexFile.get_last_values()
            new_hexfile = compute_hex_positions(start_val, end_val, linkedValues, save_opts)
            start_val = end_val
            hexFile = new_hexfile
        self._clear_file_cache()
        
    def compute_initial_chunks(self, nvalues, nchunks, save_opts):
        linkedValues = HexValues(1,1,1,0) # magic starting "previous" values
        hfile = compute_hex_positions(0, nvalues, linkedValues, save_opts)
        self.compute_next_chunks(hfile, nvalues, nchunks - 1, save_opts)


class HexValues:
    """Holds data for prime hexagon computations, including prime, position, spin, rotation data.

    This class treats the values as opaque, and we store either single
    values and also the arrays for many values
    """
    
    def __init__(self, prime, pos, spin, rot):
        self._prime = prime
        self._pos = pos
        self._spin = spin
        self._rot = rot

    @property
    def prime(self): return self._prime

    @property
    def pos(self): return self._pos

    @property
    def spin(self): return self._spin

    @property
    def rot(self): return self._rot

class HexDataFile:
    def __init__(self, filename, require_exists=True):
        self._filename = filename
        if require_exists:
            _require_infile( filename )
        startval,endval,skip = self.get_values(filename)
        self._start = startval
        self._end   = endval
        self._skip  = skip

    @property
    def filename(self): return self._filename
        
    @property
    def start(self): return self._start

    @property
    def end(self): return self._end

    @property
    def skip(self): return self._skip
    
    @staticmethod
    def get_values(filename):
        m = _HEX_FILENAME_RE.search(filename)
        if m is None or m.groups() < 2:
            raise ValueError("filename {} does not have expected format".format(filename))
        groups = m.groupdict()
        if 'skip' not in groups or groups['skip'] is None:
            groups['skip'] = 1
        startval = long(groups['start'])
        endval   = long(groups['end'])
        skip     = long(groups['skip'])
        dprint(2,"found values from filename {} = {},{},{}".format(filename,startval,endval,skip))
        return (startval, endval, skip)

    @staticmethod
    def get_filename( start, end, save_opts):
        skip = save_opts.get('skip_interval',1)
        outdir = save_opts.get('dir','output')
        basename = save_opts.get('basename','output')
        fname = "{}-{:0>20d}-{:0>20d}-{:0>6d}".format(basename,start, end, skip)
        name = os.path.join(outdir, fname)
        dprint(1,"new filename is {}".format(name))
        return name

    @classmethod
    def save_new_file( klass, start, end, data, save_opts):
        name = klass.get_filename(start, end, save_opts)
        outname = name + '.npz'
        save_binary_arrays( outname, data, save_opts)
        newFile = klass(outname)
        if save_opts.get('save_text', False):
            txtoutname = name + '.txt'
            print_text_arrays( txtoutname, data, save_opts=save_opts)
        return newFile
        
    def get_arrays(self):
        with np.load(self.filename, mmap_mode='r') as data:
            primes = data['primes']
            pos = data['pos']
            spin = data['spin']
            rot = data['rot']
        return (primes, pos, spin, rot)

    def get_zipped_arrays_iterator(self):
        primes, pos, spin, rot = self.get_arrays()
        return zipper(primes, pos, spin,rot)

    def get_sliced_zipped_arrays_iterator(self, s):
        primes, pos, spin, rot = self.get_arrays()
        # not using itertools because it does not support -1 index'ed slicing
        return zipper(primes[s], pos[s], spin[s],rot[s])
    
    def contains_value(self, n):
        return self.start <= n and n <= self.end

    def contains_any_values(self, vals):
        for val in vals:
            if self.contains_value(val):
                return True
        return False
    
    def length(self):
        return self.end - self.start

    def get_last_values(self):
        primes, pos, spin, rot = self.get_arrays()
        vals = HexValues(primes[-1], pos[-1], spin[-1], rot[-1])
        return vals
    
    def __eq__(self, other):
        return self.start == other.start and self.end == other.end

    def __cmp__(self, other):
        if other is None:
            return 1
        start_cmp = cmp(self.start, other.start)
        if 0 != start_cmp:
            return start_cmp
        else:
            return cmp(self.end, other.end)
                    
    def __hash__(self):
        "HexDataFile hash."
        return hash(self.start) ^ hash(self.end)
                    
    def __str__(self):
        return "HexDataFile(file={},start={}, end={}, skip={})".format(self.filename, self.start, self.end,self.skip)
    
    def values_around(self, vals):
        cur_vals = filter( self.contains_value, vals)
        if not cur_vals:
            IndexError("no values {} contained in HexDataFile {}".format(vals, self))
        primes, pos, spin, rot = self.get_arrays()
        idx = primes.searchsorted(cur_vals)
        for ii in idx:
            yield ii
        
    def merge(self, other, outfile):
        # assert skip values are consistent...
        # can you only merge contiguous chunks?
        # you could just insertion sort the values...
        pass

    def print_sliced_values(self, slices, ofile = sys.stdout):
        for s in slices:
            print("testing slice {} {} {}".format(s.start, s.stop, s.step))
            data = self.get_sliced_zipped_arrays_iterator(s)
            write_collection_to_file(ofile, data)

        
def get_slice_obj_from_str(slicearg):
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
    parser.add_argument('--infile', help='Input file to start processing chunks')
    parser.add_argument('--viewvalues', help='Values to view in the file as a python slice, e.g. 1:100:', action='append',)
    parser.add_argument('--find', help='Find values in files')
    parser.add_argument("-c", "--compress", help="flag to indicate whether or not to compress output files", action="store_true", default=True)
    parser.add_argument('--dir', help='Directory to read/write files to ',default='output')
    parser.add_argument('--basename', help='Basename of files to read and write', default='output')
    parser.add_argument('--logfile', help='Save messagse to this log file')
    parser.add_argument('--verbose', help='Print messages to the terminal', default=1)
    parser.add_argument('--save-text', help='Flag to save text files', action="store_true", default=False)
    parser.add_argument('--save-binary', help='Flag to save binary files',action="store_true", default=True)
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
                  'save_text' : args.save_text,
                  'dir' : args.dir,
                  'basename' : args.basename  }
    run_opts = { 'nvalues' : args.nvalues,
                 'nchunks' : args.chunks }
    
    logger.info("Save options:" + str(save_opts))

    hexAssets = HexDataAssets(save_opts)
    if args.find:
        vals = [long(v) for v in args.find.split(',')]
        vals.sort()
        hexAssets.find_values_in_npz(vals)
    elif args.viewvalues:
        _require_infile( args.infile )
        h = HexDataFile(args.infile)
        slices = [ get_slice_obj_from_str(s) for s in args.viewvalues]
        h.print_sliced_values(slices)
    else:
        if args.infile is not None:
            _require_infile( args.infile )
            start_file = HexDataFile(args.infile)
            hexAssets.compute_next_chunks(start_file, args.nvalues, args.chunks, save_opts)
        else:
            hexAssets.compute_initial_chunks(args.nvalues, args.chunks, save_opts)
        

if __name__ == "__main__":
    import sys
    main(sys.argv)



