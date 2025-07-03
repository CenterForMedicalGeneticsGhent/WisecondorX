package predict

import (
	"context"
	"fmt"
	"os"
	"strconv"

	"github.com/sbinet/npyio/npz"
	"github.com/urfave/cli/v3"
)

var (
	infile    string // Input file for prediction
	reference string // Reference file for prediction
	prefix    string // Prefix for output files
)

var PredictCmd = &cli.Command{
	Name:      "predict",
	Usage:     "Find copy number aberrations",
	ArgsUsage: "<input.npz> <reference.npz> <prefix>",
	Arguments: []cli.Argument{
		&cli.StringArg{
			Name:        "input.npz",
			UsageText:   "Input .npz file containing the sample to analyze.",
			Destination: &infile,
		},
		&cli.StringArg{
			Name:        "reference.npz",
			UsageText:   "Reference .npz file containing healthy reference samples.",
			Destination: &reference,
		},
		&cli.StringArg{
			Name:        "prefix",
			UsageText:   "Prefix for output files.",
			Destination: &prefix,
		},
	},
	Flags: []cli.Flag{
		&cli.IntFlag{
			Name:        "minrefbins",
			Usage:       "Minimum amount of sensible reference bins per target bin.",
			Value:       150,
			DefaultText: "150",
			Action: func(ctx context.Context, cmd *cli.Command, v int) error {
				if v <= 0 {
					return cli.Exit("Error: Minimum reference bins must be a positive integer", 1)
				}
				return nil
			},
		},
		&cli.IntFlag{
			Name:        "maskrepeats",
			Usage:       "Regions with distances > mean + sd * 3 will be masked. Number of masking cycles.",
			Value:       5,
			DefaultText: "5",
			Action: func(ctx context.Context, cmd *cli.Command, v int) error {
				if v < 0 {
					return cli.Exit("Error: Mask repeats must be a positive integer", 1)
				}
				return nil
			},
		},
		&cli.Float64Flag{
			Name:        "alpha",
			Usage:       "p-value cut-off for calling a CBS breakpoint.",
			Value:       0.0001,
			DefaultText: "0.0001",
			Action: func(ctx context.Context, cmd *cli.Command, v float64) error {
				if v < 0 || v > 1 {
					return cli.Exit("Error: Alpha must be between 0 and 1", 1)
				}
				return nil
			},
		},
		&cli.Float64Flag{
			Name:  "beta",
			Usage: "When beta is given, --zscore is ignored and a ratio cut-off is used to call aberrations. Beta is a number between 0 (liberal) and 1 (conservative) and is optimally close to the purity.",
			Action: func(ctx context.Context, cmd *cli.Command, v float64) error {
				if v < 0 || v > 1 {
					return cli.Exit("Error: Beta must be between 0 and 1", 1)
				}
				return nil
			},
		},
		&cli.Float64Flag{
			Name:        "zscore",
			Usage:       "z-score cut-off for aberration calling",
			Value:       5.0,
			DefaultText: "5",
			Action: func(ctx context.Context, cmd *cli.Command, v float64) error {
				if v <= 0 {
					return cli.Exit("Error: Z-score cutoff must a positive integer", 1)
				}
				return nil
			},
		},
		&cli.StringFlag{
			Name:      "blacklist",
			Usage:     "Blacklist file for regions to ignore. Format: chr(\t)start(\t)end",
			Value:     "",
			TakesFile: true,
			Action: func(ctx context.Context, cmd *cli.Command, v string) error {
				if v != "" {
					if _, err := os.Stat(v); os.IsNotExist(err) {
						return cli.Exit("Error: Blacklist file does not exist", 1)
					}
				}
				return nil
			},
		},
		&cli.StringFlag{
			Name:  "sex",
			Usage: "Force WiseCondorX to analyze the sample as the given sex. Options: M, F",
			Action: func(ctx context.Context, cmd *cli.Command, v string) error {
				if v != "M" && v != "F" {
					return cli.Exit("Error, unknown sex: "+v+". Options are: M, F", 1)
				}
				return nil
			},
		},
		&cli.IntFlag{
			Name:        "seed",
			Usage:       "Random seed for segmentation algoithm",
			Value:       42,
			DefaultText: "42",
			Action: func(ctx context.Context, cmd *cli.Command, v int) error {
				if v < 0 {
					return cli.Exit("Error: Seed must be a non-negative integer", 1)
				}
				return nil
			},
		},
	},
	Before: func(ctx context.Context, cmd *cli.Command) (context.Context, error) {
		// Check if the correct number of arguments is provided
		if cmd.Args().Len() != 3 {
			cli.ShowSubcommandHelp(cmd)
			return nil, cli.Exit("Error: Incorrect number of arguments. Expected 3 arguments while "+strconv.Itoa(cmd.Args().Len())+" were given", 1)
		}

		// Check if the input file exists
		if _, err := os.Stat(infile); os.IsNotExist(err) {
			return nil, cli.Exit("Error: Input file does not exist", 1)
		}

		// Check if the reference file exists
		if _, err := os.Stat(reference); os.IsNotExist(err) {
			return nil, cli.Exit("Error: Reference file does not exist", 1)
		}
		return ctx, nil
	},
	Action: func(ctx context.Context, cmd *cli.Command) error {
		return nil
	},
}

// SexEnum represents the available sexes for analysis.
type SexEnum string

const (
	Male   SexEnum = "M"
	Female SexEnum = "F"
)

func WcxPredict(infile string, reference string, prefix string, minrefbins int, maskrepeats int, alpha float64, beta float64, zscore float64, blacklist string, sex string, seed int) {
	// Predict chromosomal aberrations

	// Load the input npz file
	inFileData, inFileBinsize, err := LoadSampleNpzFile(infile)
	if err != nil {
		fmt.Println("Error: Unable to load file: ", infile)
		return
	}

	// Load the reference npz file
	_, refFileBinsize, err := LoadReferenceNpzFile(reference)
	if err != nil {
		fmt.Println("Error: Unable to load file: ", reference)
		return
	}

	// Scale the sample to the reference bin size
	_, err = ScaleSample(inFileData, int(inFileBinsize), int(refFileBinsize))
	if err != nil {
		fmt.Println("Error: Unable to scale sample: ", err)
		return
	}

	// Load the reference npz file
	ref, err := npz.Open(reference)
	if err != nil {
		fmt.Println("Error: Unable to open file", reference)
		return
	}
	defer ref.Close()
}
