# WisecondorX Gemini CLI Guidelines

This document provides project-specific instructions and architectural context for Gemini CLI when working on the WisecondorX codebase.

## Project Overview

WisecondorX is a bioinformatics tool for Copy Number Variation (CNV) analysis in whole-genome sequencing data. It is an evolution of the original WISECONDOR algorithm, optimized for speed and accuracy.

### Core Stack
- **Language:** Python 3.12+
- **Data Science:** NumPy, SciPy, Pandas, Scikit-learn, Matplotlib
- **Bioinformatics:** Pysam
- **CLI:** Typer
- **Supporting Logic:** R (CBS analysis via `DNAcopy` from Bioconductor, `jsonlite`, `png`)
- **Environment & Dependency Management:** [Pixi](https://pixi.sh)
- **Linting & Formatting:** Ruff (Line length: 79)
- **Testing:** Pytest (using `unittest` style tests)

## Architectural Context

### Pipeline Stages
1. **Convert:** BAM/CRAM files are converted into binned read counts stored in compressed NumPy (`.npz`) files.
2. **Newref:** Multiple `.npz` files from control samples are combined to create a reference `.npz` file. This involves PCA and gender modeling.
3. **Gender:** (Optional) Predicts the sex of a sample based on the reference model.
4. **Predict:** A test sample `.npz` is compared against a reference `.npz` to detect CNVs using Circular Binary Segmentation (CBS).
5. **RefQC:** (Supplementary) Quality control for a generated reference.

### Data Flow & "Contracts"
The project relies heavily on `.npz` files for intermediate storage. These files have strict "contracts" regarding their keys and data structures.
- **Contract Validation:** When modifying data processing logic, ensure that the `.npz` output remains compatible with subsequent stages. Refer to `tests/test_newref_contract.py` and `tests/test_predict_contract.py` for examples of these invariants.
- **Intermediate Files:** Creation of references often involves temporary "prep" and "part" files.

### CLI Structure
- `src/wisecondorx/main.py`: The main entry point using `typer` for subcommands.
- `src/wisecondorx/convert.py`: `convert` subcommand.
- `src/wisecondorx/newref.py`: `newref` subcommand.
- `src/wisecondorx/predict.py`: `predict` and `gender` subcommands.
- `src/wisecondorx/refqc.py`: `refqc` subcommand.
- `src/wisecondorx/utils.py`: Helper functions used by the pipeline stages.
- `src/wisecondorx/plotter.py`: Plotting functions for the pipeline stages.
- `src/wisecondorx/include/`: Contains R script `CBS.R` (specifically for Circular Binary Segmentation (CBS)) invoked via Python's `subprocess`.

## Development Workflow

### Pixi Tasks
Use `pixi run` for all standard development tasks to ensure environment consistency:
- `pixi run lint`: Check code quality with Ruff.
- `pixi run format`: Format code with Ruff.
- `pixi run test`: Run the test suite.
- `pixi run docker-build`: Build the Docker image.

### Code Style
- Adhere to the Ruff configuration in `pyproject.toml`.
- **Line Length:** 79 characters.
- **Import Sorting:** Ruff handles this; ensure it's run before committing.
- **String Formatting:** Use f-strings instead of `.format()` and limit outputs to 3 decimals where applicable.

### Testing Guidelines
- **Contract Tests:** When adding features that change how data is stored or processed, update the contract tests in `tests/`.
- **Functional Tests:** Aim for high coverage of the computational logic.
- **R Integration:** Be mindful that changes to R scripts may require updates to the Python wrapper logic that calls them.

## Gemini CLI Specific Instructions

- **Environment Awareness:** Always prefer `pixi run <command>` over direct execution if a pixi task exists.
- **Binary Files:** When investigating issues related to `.npz` files, use `numpy.load` within a script or via `run_shell_command` to inspect the contents.
- **Surgical Updates:** When fixing bugs, prioritize fixing the underlying logic lower level functions with minimal code changes and adding a corresponding test case.