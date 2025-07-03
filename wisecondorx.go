package main

import (
	"context"
	"log/slog"
	"net/mail"
	"os"
	"time"

	"github.com/matthdsm/wisecondorx/convert"
	"github.com/matthdsm/wisecondorx/docs"
	"github.com/matthdsm/wisecondorx/newref"
	"github.com/matthdsm/wisecondorx/predict"

	"github.com/urfave/cli/v3"
)

func main() {
	// Initialize the logger
	loglevel := slog.LevelVar{}
	logger := slog.New(slog.NewTextHandler(os.Stdout, &slog.HandlerOptions{
		Level: &loglevel,
	}))
	slog.SetDefault(logger)

	// Initialize the CLI command with its metadata
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
		Arguments: []cli.Argument{},
		ArgsUsage: "",
		Flags: []cli.Flag{
			&cli.StringFlag{
				Name:        "loglevel",
				Aliases:     []string{"l"},
				Usage:       "Set the logging level (debug, info, warn, error)",
				Value:       "info",
				DefaultText: "info",
				Action: func(ctx context.Context, cmd *cli.Command, v string) error {
					switch v {
					case "debug":
						logger.Debug("Setting log level to debug")
						loglevel.Set(slog.LevelDebug)
					case "info":
						loglevel.Set(slog.LevelInfo)
					case "warn":
						logger.Info("Setting log level to warn")
						loglevel.Set(slog.LevelWarn)
					case "error":
						logger.Info("Setting log level to error")
						loglevel.Set(slog.LevelError)
					default:
						logger.Error("Invalid log level, defaulting to info")
						loglevel.Set(slog.LevelInfo)
					}
					return nil
				},
			},
		},
		Commands: []*cli.Command{
			&convert.ConvertCmd,
			&newref.NewRefCmd,
			&predict.PredictCmd,
			&docs.BuildCmd,
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

	// Run the CLI command with the provided arguments
	Cmd.Run(context.Background(), os.Args)
}
