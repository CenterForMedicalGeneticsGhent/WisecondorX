package newref

import (
	"context"
	"fmt"
	"os"
	"runtime"

	"github.com/matthdsm/wisecondorx/utils"
	"github.com/urfave/cli/v3"
	"gonum.org/v1/gonum/floats"
)

var (
	infiles []string // Input files for the new reference
	prefix  string   // Prefix for the new reference output

)

var NewRefCmd cli.Command = cli.Command{
	Name:      "newref",
	Usage:     "Create a new reference using healthy reference samples",
	ArgsUsage: "<input1.npz> <input2.npz> ... <prefix>",
	Arguments: []cli.Argument{
		&cli.StringArgs{
			Name:        "input.npz",
			UsageText:   "Input .npz files containing healthy reference samples. Multiple files can be provided.",
			Destination: &infiles,
			Min:         0,
			Max:         -1,
		},
		&cli.StringArg{
			Name:        "prefix",
			UsageText:   "Prefix for output files. The output will be <prefix>.npz.",
			Destination: &prefix,
		},
	},
	Flags: []cli.Flag{
		&cli.BoolFlag{
			Name:        "nipt",
			Usage:       "Use NIPT reference presets",
			Value:       false,
			DefaultText: "false",
		},
		&cli.Float64Flag{
			Name:        "yfrac",
			Usage:       "Y read fraction cutoff, in order to manually define sex. Setting this to 1 will treat all samples as female",
			DefaultText: "Automatically determined",
			Action: func(ctx context.Context, cmd *cli.Command, v float64) error {
				if v < 0 || v > 1 {
					return cli.Exit("Error: Y fraction must be between 0 and 1", 1)
				}
				return nil
			},
		},
		&cli.IntFlag{
			Name:        "refsize",
			Usage:       "Number of reference locations per target. Should generally not be tweaked.",
			Value:       300,
			DefaultText: "300",
			Action: func(ctx context.Context, cmd *cli.Command, v int) error {
				if v <= 0 {
					return cli.Exit("Error: Reference size must be a positive integer", 1)
				}
				return nil
			},
		},
		&cli.IntFlag{
			Name:  "binsize",
			Usage: "Size per bin in bp, defines the resolution of the output. Multiples of existing bin sizes only.",
			Action: func(ctx context.Context, cmd *cli.Command, v int) error {
				if v <= 0 {
					return cli.Exit("Error: Binsize must be a positive integer", 1)
				}
				return nil
			},
		},
		&cli.IntFlag{
			Name:  "threads",
			Usage: "Number of threads to use for processing. Defaults to the number of available CPU cores.",
			Value: 0, // 0 means use all available cores
			Action: func(ctx context.Context, cmd *cli.Command, v int) error {
				if v < 0 {
					return cli.Exit("Error: Number of threads must be a positive integer", 1)
				} else if v > 0 {
					runtime.GOMAXPROCS(v) // Set the number of threads to use
				}
				return nil
			},
		},
	},
	Before: func(ctx context.Context, cmd *cli.Command) (context.Context, error) {
		// Check if all provided input files exist
		for _, i := range infiles {
			if _, err := os.Stat(i); os.IsNotExist(err) {
				return nil, cli.Exit("Error: Input file does not exist", 1)
			}
		}
		return ctx, nil
	},
	Action: func(ctx context.Context, cmd *cli.Command) error {
		return nil
	},
}

func WcxNewRef(infiles []string, prefix string, nipt bool, yfrac float64, refsize int, binsize int) {
	// Create a new reference using healthy reference samples

	samples := make([]map[string][]float64, len(infiles))
	// Load the input files
	for _, infile := range infiles {
		infileData, infileBinsize, err := utils.LoadSampleNpzFile(infile)
		if err != nil {
			fmt.Println("Error: Unable to load file: ", infile)
			return
		}

		// Scale the sample to the new bin size
		scaledSample, err := utils.ScaleSample(infileData, int(infileBinsize), binsize)
		if err != nil {
			fmt.Println("Error: Unable to scale sample: ", err)
			return
		}
		samples = append(samples, scaledSample)
	}

	sexes, _ := trainSexDeterminationModel(samples, yfrac)

	// Do some sanity checks if a NIPT reference is requested
	if nipt == true {
		femaleCount := 0
		for _, sex := range sexes {
			if sex == "F" {
				femaleCount++
			}
		}
		if femaleCount < 5 {
			fmt.Println("A NIPT reference should have at least 5 female feti samples")
			return
		}
	}
	if nipt == false {
		for i, sample := range samples {
			samples[i] = correctSex(sample, sexes[i])
		}
	}
}

func trainSexDeterminationModel(samples []map[string][]float64, yfrac float64) (sexes []string, cutoff float64) {
	// Train the sex determination model
	sexes = make([]string, len(samples))
	// Create the training data
	//sexes := make([]string, len(samples))
	yfracs := make([]float64, len(samples))

	// Calculate the Y fraction for each sample
	for _, sample := range samples {
		sliceSumY := floats.Sum(sample["chrY"])
		sliceSumAll := 0.0
		for _, bins := range sample {
			sliceSumAll += floats.Sum(bins)
		}
		yfrac := float64(sliceSumY) / sliceSumAll
		yfracs = append(yfracs, yfrac)
	}

	// Train the Guassian Mixture Model
	return sexes, cutoff
}

func correctSex(sample map[string][]float64, sex string) map[string][]float64 {
	//Levels gonosomal reads with the one at the autosomes.
	if sex == "M" {
		correctedChrX := make([]float64, len(sample["chrX"]))
		for _, v := range sample["chrX"] {
			correctedChrX = append(correctedChrX, v*2)
		}
		sample["chrX"] = correctedChrX
		correctedChrY := make([]float64, len(sample["chrY"]))
		for _, v := range sample["chrY"] {
			correctedChrY = append(correctedChrY, v*2)
		}
		sample["chrY"] = correctedChrY
	}
	return sample
}

func calculateMask(samples []map[string][]int32) (totalMask int, binsPerChromosome map[string][]int32) {
	// Calculate the mask for the samples
	totalMask = 0
	binsPerChromosome = make(map[string][]int32)

	// Get all chromosomes
	chromosomes := make([]string, 0)
	for chr, _ := range samples[0] {
		chromosomes = append(chromosomes, chr)
	}
	// Get max length of the bin slice among all chromosomes
	maxLen := 0
	for _, sample := range samples {
		for _, bins := range sample {
			if len(bins) > maxLen {
				maxLen = len(bins)
			}
		}
	}

	return totalMask, binsPerChromosome
}
