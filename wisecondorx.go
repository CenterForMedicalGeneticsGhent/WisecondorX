package main

import (
	"context"
	"net/mail"
	"os"
	"strconv"
	"time"

	"github.com/urfave/cli/v3"
)

func main() {
	Cmd := &cli.Command{
		Name:    "wisecondorx",
		Version: "2.0.0",
		Authors: []any{
			&mail.Address{
				Name:    "Matthias De Smet",
				Address: "matthias.desmet@ugent.be",
			},
			&mail.Address{
				Name:    "CMGG ICT Team",
				Address: "ict.cmgg@uzgent.be",
			},
		},
		Copyright: "Copyright (c) " + time.Now().Format("2006") + " Center for Medical Genetics Ghent, Ghent University Hospital",
		Usage:     "an evolved WISECONDOR",
		UsageText: "wisecondorx [global options] command [command options] [arguments...]",
		ArgsUsage: "",
		Commands: []*cli.Command{
			{
				Name:      "convert",
				Usage:     "Convert and filter aligned reads to .npz format",
				UsageText: "wisecondorx convert [options] <input.bam/cram> <prefix>",
				ArgsUsage: "<input.bam/cram> <prefix>",
				Flags: []cli.Flag{
					&cli.StringFlag{
						Name:    "reference",
						Aliases: []string{"r"},
						Usage:   "Reference genome file for cram conversion.",
					},
					&cli.IntFlag{
						Name:        "binsize",
						Aliases:     []string{"b"},
						Usage:       "Bin size (bp)",
						DefaultText: "5000",
						Value:       5000,
					},
					&cli.BoolFlag{
						Name:    "normdup",
						Aliases: []string{"no-remove-duplicates"},
						Usage:   "Do not remove duplicates",
						Value:   false,
					},
				},
				Before: func(ctx context.Context, cmd *cli.Command) (context.Context, error) {
					// Check if the correct number of arguments is provided
					if cmd.Args().Len() != 2 {
						cli.ShowSubcommandHelp(cmd)
						return nil, cli.Exit("Error: Incorrect number of arguments. Expected 2 arguments while "+strconv.Itoa(cmd.Args().Len())+" were given", 1)
					}

					// Check if the input file exists
					if _, err := os.Stat(cmd.Args().Get(0)); os.IsNotExist(err) {
						return nil, cli.Exit("Error: Input file does not exist", 1)
					}

					// Check if the reference file exists if provided
					if cmd.String("reference") != "" {
						if _, err := os.Stat(cmd.String("reference")); os.IsNotExist(err) {
							return nil, cli.Exit("Error: Reference file does not exist", 1)
						}
					}
					// Check if the binsize is a positive integer
					if cmd.Int("binsize") <= 0 {
						return nil, cli.Exit("Error: Binsize must be a positive integer", 1)
					}
					return nil, nil
				},
				Action: func(ctx context.Context, cmd *cli.Command) error {
					return nil
				},
			},
			{
				Name:      "newref",
				Usage:     "Create a new reference using healthy reference samples",
				ArgsUsage: "<input1.npz> <input2.npz> ... <prefix>",
				Flags: []cli.Flag{
					&cli.BoolFlag{
						Name:        "nipt",
						Usage:       "Use NIPT reference presets",
						Value:       false,
						DefaultText: "false",
					},
					&cli.Float64Flag{
						Name:        "yfrac",
						Usage:       "Manually set the Y-fraction cutoff to determine sex",
						DefaultText: "Automatically determined",
					},
					&cli.IntFlag{
						Name:        "refsize",
						Usage:       "Number of reference locations per target",
						Value:       300,
						DefaultText: "300",
					},
					&cli.IntFlag{
						Name:  "binsize",
						Usage: "Scale samples to this bin size (bp), multiples of existing bin sizes only",
					},
				},
				Before: func(ctx context.Context, cmd *cli.Command) (context.Context, error) {
					// Check if all provided input files exist
					for i := 0; i < cmd.Args().Len()-1; i++ {
						if _, err := os.Stat(cmd.Args().Get(i)); os.IsNotExist(err) {
							return nil, cli.Exit("Error: Input file does not exist", 1)
						}
					}

					// Check if fraction is between 0 and 1
					if cmd.Float64("yfrac") < 0 || cmd.Float64("yfrac") > 1 {
						return nil, cli.Exit("Error: Y fraction must be between 0 and 1", 1)
					}
					return ctx, nil
				},
				Action: func(ctx context.Context, cmd *cli.Command) error {
					return nil
				},
			},
			{
				Name:      "predict",
				Usage:     "Find copy number aberrations",
				ArgsUsage: "<input.npz> <reference.npz> <prefix>",
				Flags: []cli.Flag{
					&cli.IntFlag{
						Name:        "minrefbins",
						Usage:       "Minimum amount of sensible reference bins per target bin.",
						Value:       150,
						DefaultText: "150",
					},
					&cli.IntFlag{
						Name:        "maskrepeats",
						Usage:       "Regions with distances > mean + sd * 3 will be masked. Number of masking cycles.",
						Value:       5,
						DefaultText: "5",
					},
					&cli.Float64Flag{
						Name:        "alpha",
						Usage:       "p-value cut-off for calling a CBS breakpoint.",
						Value:       0.0001,
						DefaultText: "0.0001",
					},
					&cli.Float64Flag{
						Name:  "beta",
						Usage: "When beta is given, --zscore is ignored and a ratio cut-off is used to call aberrations. Beta is a number between 0 (liberal) and 1 (conservative) and is optimally close to the purity.",
					},
					&cli.Float64Flag{
						Name:        "zscore",
						Usage:       "z-score cut-off for aberration calling",
						Value:       5.0,
						DefaultText: "5",
					},
					&cli.StringFlag{
						Name:  "blacklist",
						Usage: "Blacklist file for regions to ignore. Format: chr(\t)start(\t)end",
					},
					&cli.StringFlag{
						Name:  "sex",
						Usage: "Force WiseCondorX to analyze the sample as the given sex. Options: M, F",
					},
					&cli.IntFlag{
						Name:        "seed",
						Usage:       "Random seed for segmentation algoithm",
						Value:       0,
						DefaultText: "0",
					},
				},
				Before: func(ctx context.Context, cmd *cli.Command) (context.Context, error) {
					// Check if the correct number of arguments is provided
					if cmd.Args().Len() != 3 {
						cli.ShowSubcommandHelp(cmd)
						return nil, cli.Exit("Error: Incorrect number of arguments. Expected 3 arguments while "+strconv.Itoa(cmd.Args().Len())+" were given", 1)
					}

					// Check if the input file exists
					if _, err := os.Stat(cmd.Args().Get(0)); os.IsNotExist(err) {
						return nil, cli.Exit("Error: Input file does not exist", 1)
					}

					// Check if the reference file exists
					if _, err := os.Stat(cmd.Args().Get(1)); os.IsNotExist(err) {
						return nil, cli.Exit("Error: Reference file does not exist", 1)
					}

					// Check if the blacklist file exists
					if cmd.String("blacklist") != "" {
						if _, err := os.Stat(cmd.String("blacklist")); os.IsNotExist(err) {
							return nil, cli.Exit("Error: Blacklist file does not exist", 1)
						}
					}

					// Check if zscore is bigger than 0
					if cmd.Float64("zscore") <= 0 {
						return nil, cli.Exit("Error: Z-score must be bigger than 0", 1)
					}

					// Check if alpha is between 0 and 1
					if cmd.Float64("alpha") < 0 || cmd.Float64("alpha") > 1 {
						return nil, cli.Exit("Error: Alpha must be between 0 and 1", 1)
					}

					// Check if beta is between 0 and 1
					if cmd.Float64("beta") != 0 && (cmd.Float64("beta") < 0 || cmd.Float64("beta") > 1) {
						return nil, cli.Exit("Error: Beta must be between 0 and 1", 1)
					}
					return ctx, nil
				},
				Action: func(ctx context.Context, cmd *cli.Command) error {
					return nil
				},
			},
		},
		EnableShellCompletion: true,
		HideHelp:              false,
		HideVersion:           false,
		ShellComplete: func(ctx context.Context, cmd *cli.Command) {
			// Custom shell completion logic can be added here if needed
			// For now, we just call the default completion
			cli.DefaultAppComplete(ctx, cmd)
		},
		Action: func(ctx context.Context, cmd *cli.Command) error {
			cli.ShowAppHelp(cmd)
			return nil
		},
	}

	Cmd.Run(context.Background(), os.Args)
}
