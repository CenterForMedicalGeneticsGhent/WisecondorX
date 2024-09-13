package main

import (
	"log"
	"os"

	"github.com/CenterForMedicalGeneticsGhent/WisecondorX/internal/wisecondorx"
	"github.com/urfave/cli/v2"
)

func main() {
	app := &cli.App{
		Name:    "wisecondorx",
		Version: "2.0.0",
		Usage:   "an evolved WISECONDOR",
	}

	app.Commands = []*cli.Command{
		{
			Name:      "convert",
			Aliases:   []string{"c"},
			Usage:     "Convert and filter a aligned reads to .wcx format",
			ArgsUsage: "<input.bam> <prefix>",
			Flags: []cli.Flag{
				// &cli.StringFlag{
				// 	Name:    "reference",
				// 	Aliases: []string{"r"},
				// 	Usage:   "Reference genome file for cram conversion.",
				// },
				&cli.IntFlag{
					Name:        "binsize",
					Aliases:     []string{"b"},
					Usage:       "Bin size (bp)",
					DefaultText: "5000",
					Value:       5000,
				},
				&cli.BoolFlag{
					Name:  "normdup",
					Usage: "Do not remove duplicates",
					Value: false,
				},
			},
			Before: func(cCtx *cli.Context) error {
				// Check if the correct number of arguments is provided
				if cCtx.Args().Len() != 2 {
					return cli.ShowSubcommandHelp(cCtx)
				}

				// Check if the input file exists
				if _, err := os.Stat(cCtx.Args().Get(0)); os.IsNotExist(err) {
					return cli.Exit("Error: Input file does not exist", 1)
				}

				// Check if the input file ends with .bam
				if cCtx.Args().Get(0)[len(cCtx.Args().Get(0))-4:] != ".bam" {
					return cli.Exit("Error: Input file be in bam format", 1)
				}

				// Check if the reference file exists
				// if cCtx.String("reference") != "" {
				// 	if _, err := os.Stat(cCtx.String("reference")); os.IsNotExist(err) {
				// 		return cli.Exit("Reference file does not exist", 1)
				// 	}
				// }
				return nil
			},
			Action: func(cCtx *cli.Context) error {
				infile := cCtx.Args().Get(0)
				prefix := cCtx.Args().Get(1)
				binsize := cCtx.Int("binsize")
				normdup := cCtx.Bool("normdup")
				reference := cCtx.String("reference")
				wisecondorx.WcxConvert(infile, prefix, binsize, !normdup, reference)
				return nil
			},
		},
		{
			Name:      "newref",
			Aliases:   []string{"n"},
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
			Before: func(cCtx *cli.Context) error {
				// Check if all provided input files exist
				for i := 0; i < cCtx.Args().Len()-1; i++ {
					if _, err := os.Stat(cCtx.Args().Get(i)); os.IsNotExist(err) {
						return cli.Exit("Error: Input file does not exist", 1)
					}
				}

				// Check if fraction is between 0 and 1
				if cCtx.Float64("yfrac") < 0 || cCtx.Float64("yfrac") > 1 {
					return cli.Exit("Error: Y fraction must be between 0 and 1", 1)
				}
				return nil
			},
			Action: func(cCtx *cli.Context) error {
				infiles := cCtx.Args().Slice()[0 : cCtx.Args().Len()-1]
				prefix := cCtx.Args().Get(cCtx.Args().Len() - 1)
				nipt := cCtx.Bool("nipt")
				yfrac := cCtx.Float64("yfrac")
				refsize := cCtx.Int("refsize")
				binsize := cCtx.Int("binsize")
				wisecondorx.WcxNewRef(infiles, prefix, nipt, yfrac, refsize, binsize)
				return nil
			},
		},
		{
			Name:      "predict",
			Aliases:   []string{"p"},
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
			Action: func(cCtx *cli.Context) error {
				infile := cCtx.Args().Get(0)
				ref := cCtx.Args().Get(1)
				prefix := cCtx.Args().Get(2)
				minrefbins := cCtx.Int("minrefbins")
				maskrepeats := cCtx.Int("maskrepeats")
				alpha := cCtx.Float64("alpha")
				beta := cCtx.Float64("beta")
				zscore := cCtx.Float64("zscore")
				blacklist := cCtx.String("blacklist")
				sex := cCtx.String("sex")
				seed := cCtx.Int("seed")

				wisecondorx.WcxPredict(infile, ref, prefix, minrefbins, maskrepeats, alpha, beta, zscore, blacklist, sex, seed)
				return nil
			},
		},
	}

	if err := app.Run(os.Args); err != nil {
		log.Fatal(err)
	}
}
