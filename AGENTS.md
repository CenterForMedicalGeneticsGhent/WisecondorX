# WisecondorX Agent Guidelines

This document provides project-specific instructions and architectural context for agents when working on the WisecondorX codebase.

## Project Overview

WisecondorX is a bioinformatics tool for Copy Number Variation (CNV) analysis in whole-genome sequencing data. It is an evolution of the original WISECONDOR algorithm, optimized for speed and accuracy.

### Core Stack

- **Language:** Python 3.12+
- **Data Science:** NumPy, SciPy, Polars, Scikit-learn, PyArrow, Ruptures
- **Bioinformatics:** Pysam
- **Environment & Dependency Management:** [Pixi](https://pixi.sh)
- **Linting**: Prek  & Ruff
- **Testing:** Pytest (using `unittest` style tests)

## Architectural Context

### Pipeline Stages

1. **Convert:** BAM/CRAM files are converted into binned read counts stored in or parquet (`.parquet`) file.
2. **Newref:** Multiple `.parquet` files from control samples are combined to create a reference parquet file. This involves PCA and gender modeling.
3. **Predict:** A test sample `.parquet` is compared against a reference `.parquet` to detect CNVs.

### CLI Structure

- `src/wisecondorx/cli.py`: The main entry point using `typer` for subcommands.
- `src/wisecondorx/predict.py`: Functional logic for the predict stage.
- `src/wisecondorx/convert.py`: Functional logic for the convert stage.
- `src/wisecondorx/reference.py`: Functional logic for creating reference files.

## Development Workflow

### Pixi Tasks

Use `pixi run` for all standard development tasks to ensure environment consistency:

- `pixi run lint`: Check code quality with Ruff.
- `pixi run format`: Format code with Ruff.
- `pixi run test`: Run the test suite.
- `pixi run hooks`: Run pre-commit hooks (managed by `prek`).

### Code Style

- Adhere to the Ruff configuration in `pyproject.toml`.
- **Import Sorting:** Ruff handles this; ensure it's run before committing.
- **Typing:** Use maximal type annotations

### Testing Guidelines

- **Contract Tests:** When adding features that change how data is stored or processed, update the contract tests in `tests/`.
- **Functional Tests:** Aim for high coverage of the computational logic.

## Agent Specific Instructions

- **Environment Awareness:** Always prefer `pixi run <command>` over direct execution if a pixi task exists.
- **Surgical Updates:** When fixing bugs, prioritize small, targeted changes and adding a corresponding test case.
- **Documentation:** Update this document and README if you make changes that affect the overall architecture or development workflow.