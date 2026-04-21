# Sandbox for Outer Loop Analysis (Sola)

This MATLAB library solves equality constrained optimization problems in reduced space with derivative (1st and 2nd) based optimization computed via adjoints.
In addition, it automatically enables hyper-differential sensitivity analysis with respect to model discrepancy.

## File Structure

- [**src/**](./src/): Source code defining class interfaces and codes for outer loop analysis.
- [**examples/**](./examples/): Applications / examples with real physics.
- [**tests/**](./tests/): Unit tests and toy examples for verification.
- [**docs/**](./docs/): Documentation, powered by `sphinx`.

## Formatting MATLAB Code (Style Enforcement)

Code in this repository is formatted using the `mh_style` command line tool provided in the [MISS_HIT](https://misshit.org/) package by Florian Schanda.

### Installation

```bash
pip install miss_hit
```

### Usage

Use the following commands to format the MATLAB source code.

```bash
mh_style            # List style problems in detail without changing any files.
mh_style --brief    # List style problems briefly without changing any files.
mh_style --fix      # Fix style problems, overwriting offending files.
```

The following will format MATLAB code blocks in the documentation.

```bash
cd docs
python3 format_code_blocks.py
```

To ensure proper formatting before committing changes, install a git `pre-commit` hook with the following command.

```bash
cp .pre-commit-formatting-hook .git/hooks/pre-commit
```

Now `git commit` automatically runs the formatter and cleans the code.
A commit is **aborted** if the formatter makes any changes; in this case, review / add the changes and commit again.

To run the formatter without a commit:

```bash
./.pre-commit-formatting-hook
```
