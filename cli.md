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

**--binsize, -b**="": Size per bin in bp. (default: 5000)

**--exclude-contigs, -e**="": Regex pattern to exclude certain contigs from conversion (default: ^(.*_alt|.*_decoy|.*_random|chrUn.*|HLA.*|chrM|chrEBV)$)

**--gonosomes, -g**="": Gonosome chromosomes to be used in the analysis. Should generally not be tweaked (default: [chrX chrY])

**--help, -h**: show help

**--normdup, --no-remove-duplicates**: Do not remove duplicates

**--reference, -r**="": Reference genome fasta for cram conversion.

### help, h

Shows a list of commands or help for one command

## newref

Create a new reference using healthy reference samples

**--binsize**="": Size per bin in bp, defines the resolution of the output. Multiples of existing bin sizes only. (default: 0)

**--help, -h**: show help

**--nipt**: Use NIPT reference presets

**--refsize**="": Number of reference locations per target. Should generally not be tweaked. (default: 300)

**--threads**="": Number of threads to use for processing. Defaults to the number of available CPU cores. (default: 0)

**--yfrac**="": Y read fraction cutoff, in order to manually define sex. Setting this to 1 will treat all samples as female (default: 0)

### help, h

Shows a list of commands or help for one command

## predict

Find copy number aberrations

**--alpha**="": p-value cut-off for calling a CBS breakpoint. (default: 0.0001)

**--beta**="": When beta is given, --zscore is ignored and a ratio cut-off is used to call aberrations. Beta is a number between 0 (liberal) and 1 (conservative) and is optimally close to the purity. (default: 0)

**--blacklist**="": Blacklist for masking additional regions; requires headerless .bed file. This is particularly useful when the reference set is too small to recognize some obvious loci, such as centromeres

**--help, -h**: show help

**--maskrepeats**="": Bins with distances > (mean + sd * 3) in the reference will be masked. This parameter represents the number of masking cycles and defines the stringency of the blacklist (default: 5)

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
