# NAME

wisecondorx - an evolved WISECONDOR

# SYNOPSIS

wisecondorx

```
[--help|-h]
[--loglevel|-l]=[value]
[--version|-v]
```

**Usage**:

```
wisecondorx [global options] command [command options] [arguments...]
```

# GLOBAL OPTIONS

**--help, -h**: show help

**--loglevel, -l**="": Set the logging level (debug, info, warn, error) (default: info)

**--version, -v**: print the version


# COMMANDS

## convert

Convert and filter aligned reads to .npz format

>wisecondorx convert [options] <input.bam/cram> <prefix>

**--binsize, -b**="": Bin size (bp) (default: 5000)

**--exclude-contigs, -e**="": Glob pattern to exclude certain contigs from conversion (default: {*_alt,*_decoy,*_random,chrUn*,HLA*,chrM,chrEBV})

**--gonosomes, -g**="": Overwrite default gonosomes. (default: [chrX chrY])

**--help, -h**: show help

**--normdup, --no-remove-duplicates**: Do not remove duplicates

**--reference, -r**="": Reference genome file for cram conversion.

### help, h

Shows a list of commands or help for one command

## newref

Create a new reference using healthy reference samples

**--binsize**="": Scale samples to this bin size (bp), multiples of existing bin sizes only (default: 0)

**--help, -h**: show help

**--nipt**: Use NIPT reference presets

**--refsize**="": Number of reference locations per target (default: 300)

**--yfrac**="": Manually set the Y-fraction cutoff to determine sex (default: 0)

### help, h

Shows a list of commands or help for one command

## predict

Find copy number aberrations

**--alpha**="": p-value cut-off for calling a CBS breakpoint. (default: 0.0001)

**--beta**="": When beta is given, --zscore is ignored and a ratio cut-off is used to call aberrations. Beta is a number between 0 (liberal) and 1 (conservative) and is optimally close to the purity. (default: 0)

**--blacklist**="": Blacklist file for regions to ignore. Format: chr(	)start(	)end

**--help, -h**: show help

**--maskrepeats**="": Regions with distances > mean + sd * 3 will be masked. Number of masking cycles. (default: 5)

**--minrefbins**="": Minimum amount of sensible reference bins per target bin. (default: 150)

**--seed**="": Random seed for segmentation algoithm (default: 42)

**--sex**="": Force WiseCondorX to analyze the sample as the given sex. Options: M, F

**--zscore**="": z-score cut-off for aberration calling (default: 5)

### help, h

Shows a list of commands or help for one command

## docs, d

Generate CLI documentation

**--help, -h**: show help

### help, h

Shows a list of commands or help for one command

## help, h

Shows a list of commands or help for one command
