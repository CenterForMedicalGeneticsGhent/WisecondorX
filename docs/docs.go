package docs

import (
	"context"
	"os"

	docs "github.com/urfave/cli-docs/v3"
	"github.com/urfave/cli/v3"
)

var BuildCmd = cli.Command{
	Name:    "docs",
	Aliases: []string{"d"},
	Usage:   "Generate CLI documentation",
	Action: func(ctx context.Context, cmd *cli.Command) error {
		md, err := docs.ToMarkdown(cmd.Root())
		if err != nil {
			panic(err)
		}
		fi, err := os.Create("cli.md")
		if err != nil {
			panic(err)
		}
		defer fi.Close()
		if _, err := fi.WriteString(md); err != nil {
			panic(err)
		}
		return nil
	},
}
