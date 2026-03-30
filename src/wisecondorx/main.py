#!/usr/bin/env python3

import logging
import warnings
import typer
from wisecondorx.convert import wcx_convert
from wisecondorx.newref import wcx_newref
from wisecondorx.predict import wcx_sex, wcx_predict
from wisecondorx import __version__

VERSION: str = __version__


def setup_logging(loglevel: str = "INFO") -> None:
    logging.basicConfig(
        format="[%(levelname)s - %(asctime)s]: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=getattr(logging, loglevel.upper(), None),
    )


app = typer.Typer(
    name="wisecondorx",
    help="WisecondorX: Copy Number Aberration detection from Whole Genome Sequencing data.",
    add_completion=False,
)
app.command(name="convert")(wcx_convert)
app.command(name="newref")(wcx_newref)
app.command(name="gender")(wcx_sex)
app.command(name="predict")(wcx_predict)


@app.callback()
def main_callback(
    loglevel: str = typer.Option(
        "INFO",
        "--loglevel",
        help="Logging level (info, warning, debug, error, critical)",
    ),
) -> None:
    warnings.filterwarnings("ignore")
    setup_logging(loglevel=loglevel)


if __name__ == "__main__":
    app()
