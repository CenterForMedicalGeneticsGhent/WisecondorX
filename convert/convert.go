package convert

import (
	"context"
	"fmt"
	"io"
	"log/slog"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/sbinet/npyio/npz"
	"github.com/urfave/cli/v3"

	"v.io/v23/glob"
)

var (
	infile string
	prefix string
)

var ConvertCmd cli.Command = cli.Command{
	Name:      "convert",
	Usage:     "Convert and filter aligned reads to .npz format",
	UsageText: "wisecondorx convert [options] <input.bam/cram> <prefix>",
	ArgsUsage: "<input.bam/cram> <prefix>",
	Arguments: []cli.Argument{
		&cli.StringArg{
			Name:        "bam/cram",
			UsageText:   "Input BAM or CRAM file to convert. If a CRAM file is provided, a reference genome must be specified using the --reference flag.",
			Destination: &infile,
		},
		&cli.StringArg{
			Name:        "prefix",
			UsageText:   "Prefix for output files. The output will be <prefix>.npz",
			Destination: &prefix,
		},
	},
	Flags: []cli.Flag{
		&cli.StringFlag{
			Name:      "reference",
			Aliases:   []string{"r"},
			Usage:     "Reference genome file for cram conversion.",
			TakesFile: true,
			Action: func(ctx context.Context, cmd *cli.Command, v string) error {
				if _, err := os.Stat(v); os.IsNotExist(err) {
					return cli.Exit("Error: Reference file does not exist", 1)
				}
				return nil
			},
		},
		&cli.IntFlag{
			Name:        "binsize",
			Aliases:     []string{"b"},
			Usage:       "Bin size (bp)",
			DefaultText: "5000",
			Value:       5000,
			Action: func(ctx context.Context, cmd *cli.Command, v int) error {
				if v <= 0 {
					return cli.Exit("Error: Binsize must be a positive integer", 1)
				}
				return nil
			},
		},
		&cli.BoolFlag{
			Name:    "normdup",
			Aliases: []string{"no-remove-duplicates"},
			Usage:   "Do not remove duplicates",
			Value:   false,
		},
		&cli.StringFlag{
			Name:        "exclude-contigs",
			Aliases:     []string{"e"},
			Usage:       "Glob pattern to exclude certain contigs from conversion",
			DefaultText: "{*_alt,*_decoy,*_random,chrUn*,HLA*,chrM,chrEBV}",
			Value:       "{*_alt,*_decoy,*_random,chrUn*,HLA*,chrM,chrEBV}",
		},
		&cli.StringSliceFlag{
			Name:        "gonosomes",
			Aliases:     []string{"g"},
			Usage:       "Overwrite default gonosomes.",
			DefaultText: "[chrX, chrY]",
			Value:       []string{"chrX", "chrY"},
		},
	},
	Before: func(ctx context.Context, cmd *cli.Command) (context.Context, error) {
		// Check if the correct number of arguments is provided
		if cmd.Args().Len() != 2 {
			cli.ShowSubcommandHelp(cmd)
			return nil, cli.Exit("Error: Incorrect number of arguments. Expected 2 arguments while "+strconv.Itoa(cmd.Args().Len())+" were given", 1)
		}

		// Check if the input file exists
		if _, err := os.Stat(infile); os.IsNotExist(err) {
			return nil, cli.Exit("Error: Input file does not exist", 1)
		}
		// Check if the input file has an index
		// if file extension is .bam, check for bai or csi index
		if ext := filepath.Ext(infile); ext == ".bam" {
			if _, err := os.Stat(infile + ".bai"); os.IsNotExist(err) {
				if _, err := os.Stat(infile + ".csi"); os.IsNotExist(err) {
					return nil, cli.Exit("Error: BAM file index does not exist. Please provide a .bai or .csi file", 1)
				}
			}
		}
		// if file extension is .cram, check if a reference is given and check for crai index
		if ext := filepath.Ext(infile); ext == ".cram" {
			if cmd.String("reference") == "" {
				return nil, cli.Exit("Error: Reference genome must be specified for CRAM files using the --reference flag", 1)
			}
			if _, err := os.Stat(infile + ".crai"); os.IsNotExist(err) {
				return nil, cli.Exit("Error: CRAM file index does not exist. Please provide a .crai file", 1)
			}
		}

		return nil, nil
	},
	Action: func(ctx context.Context, cmd *cli.Command) error {
		return nil
	},
}

type reader struct {
	io.ReadCloser
	cmd *exec.Cmd
}

func (r *reader) Close() error {
	if err := r.cmd.Wait(); err != nil {
		return err
	}
	return r.ReadCloser.Close()
}

// NewReader returns a bam.Reader from any path that samtools can read.
// This is a patch to enable cram compatibilty
func NewReader(path string, rd int, fasta string) (*bam.Reader, error) {
	cmd := exec.Command("samtools", "view", "-T", fasta, "-b", "-u", "-h", path)
	cmd.Stderr = os.Stderr
	pipe, err := cmd.StdoutPipe()
	if err != nil {
		return nil, err
	}
	if err = cmd.Start(); err != nil {
		pipe.Close()
		return nil, err
	}
	cr := &reader{ReadCloser: pipe, cmd: cmd}
	return bam.NewReader(cr, rd)
}

func WcxConvert(infile string, prefix string, binsize int, RemoveDup bool, reference string, excludeContigs string, gonosomes []string) {
	// Convert and filter aligned reads to .wcx format
	b, err := NewReader(infile, 0, reference)
	if err != nil {
		fmt.Printf("Error: Unable to read bam file %s", infile)
		return
	}
	defer b.Close()

	// Generate a map of chromosome names and the number of bins per chromosome.
	// Exclude chromosomes defined in the 'excludeContigs` glob pattern
	if excludeContigs == "" {
		excludeContigs = "{*_alt,*_decoy,*_random,chrUn*,HLA*,chrM,chrEBV}"
	}
	excludeContigsGlob, err := glob.Parse(excludeContigs)
	if err != nil {
		slog.Error("Unable to parse contig exclusion glob")
	}

	binsPerChromosome := make(map[string][]float64)
	for _, ref := range b.Header().Refs() {
		// skip the chromosome if not autosomes or gonosomes
		if excludeContigsGlob.Head().Match(ref.Name()) {
			continue
		}
		binsPerChromosome[ref.Name()] = make([]float64, ref.Len()/binsize+1)
	}

	currentChromosome := ""
	larp := -1
	larp2 := -1

	// Loop over the reads
	for {
		samRecord, err := b.Read()
		if err != nil {
			break
		}

		// Skip the read if contig is excluded
		if excludeContigsGlob.Head().Match(samRecord.Ref.Name()) {
			continue
		}

		if currentChromosome != samRecord.Ref.Name() {
			slog.Info("Processing chromosome: " + samRecord.Ref.Name())
			currentChromosome = samRecord.Ref.Name()
		}

		// Skip the read if duplicate && RemoveDup is true
		if samRecord.Flags&sam.Duplicate != 0 && RemoveDup == true {
			continue
		}

		// read is unpaired
		if samRecord.Flags&sam.Paired == 0 {
			if larp != samRecord.Pos {
				if samRecord.MapQ >= 1 {
					binsPerChromosome[samRecord.Ref.Name()][samRecord.Pos/binsize]++
				}
			}
			larp = samRecord.Pos
			// go to next read
			continue
		}

		// read is paired
		if samRecord.Flags&sam.Paired != 0 {
			// read is properly paired
			if samRecord.Flags&sam.ProperPair != 0 {
				if larp != samRecord.Pos && larp2 != samRecord.MatePos {
					if samRecord.MapQ >= 1 {
						binsPerChromosome[samRecord.Ref.Name()][samRecord.Pos/binsize]++
					}
				}
			}
			larp = samRecord.Pos
			larp2 = samRecord.MatePos
		}
	}

	// Write the output file
	npz, err := npz.Create(prefix + ".npz")
	if err != nil {
		fmt.Printf("Error: Unable to create file %s.npz: %s", prefix, err)
	}
	defer npz.Close()

	// Write bins per chromosome
	for chromosome, bins := range binsPerChromosome {
		err = npz.Write(chromosome, bins)
		if err != nil {
			fmt.Printf("Error: Unable to write to file %s.npz: %s", prefix, err)
		}
	}
	// Write binsize
	err = npz.Write("binsize", float64(binsize))
	if err != nil {
		fmt.Printf("Error: Unable to write to file %s.npz: %s", prefix, err)
	}

	err = npz.Close()
	if err != nil {
		fmt.Printf("Error: Unable to close file %s.npz: %s", prefix, err)
	}
}
