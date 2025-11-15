# Pixi Development Guide for HippUnfold

This guide walks you through setting up and using Pixi for HippUnfold development.

## What is Pixi?

Pixi is a fast package manager for Python projects that uses `conda` under the hood but with improved speed and reproducibility. It manages dependencies defined in `pyproject.toml` and locks them in `pixi.lock`.

## Initial Setup (One-time)

### 1. Install Pixi

If you don't have Pixi installed, run:

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

Then restart your terminal or source your shell configuration:

```bash
source ~/.bashrc  # for bash
# or
source ~/.zshrc   # for zsh
```

Verify installation:

```bash
pixi --version
```

### 2. Navigate to the HippUnfold directory

```bash
cd /nfs/khan/trainees/msalma29/hippunfold_dev/hippunfold
```

### 3. Install the project environment

```bash
pixi install
```

This command will:
- Read `pyproject.toml` to understand your project's dependencies
- Download all required packages from conda-forge and bioconda
- Create a locked `pixi.lock` file (reproducible across machines)
- Set up the development environment locally in `.pixi/envs/`

**This only needs to be done once!**

## Daily Development Workflow

### Using the Development Environment

HippUnfold has three environments defined in `pyproject.toml`:

- **`default`**: Runtime dependencies only (minimal, for production use)
- **`dev`**: Runtime + development dependencies (recommended for development)
- **`dev-only`**: Development tools only (useful for CI/linting)

### Running Commands in the Dev Environment

Use the `-e dev` flag to specify the development environment:

```bash
# Run HippUnfold with help
pixi run -e dev hippunfold --help

# Run HippUnfold quick version
pixi run -e dev hippunfold-quick --help

# Run Python directly
pixi run -e dev python script.py

# Run tests
pixi run -e dev pytest
```

### Interactive Shell (Recommended)

Activate the environment interactively for a more seamless workflow:

```bash
pixi shell -e dev
```

Now you're inside the environment and can run commands directly without `pixi run`:

```bash
# Inside the pixi shell
hippunfold --help
hippunfold-quick --help
python script.py
pytest

# Exit the environment
exit
```

## Quality Assurance

HippUnfold includes pre-configured quality tasks in `pyproject.toml`:

### Check Code Quality

```bash
pixi run -e dev quality_check
```

This runs:
- `isort` - Import sorting
- `black` - Code formatting
- `snakefmt` - Snakefile formatting

All in check mode (won't modify files).

### Fix Code Quality Issues Automatically

```bash
pixi run -e dev quality_fix
```

This runs the same tools but **modifies files** to fix formatting and import issues.

## Adding New Dependencies

### Add to Runtime Dependencies

```bash
pixi add package-name
```

### Add to Development Dependencies Only

```bash
pixi add --feature dev package-name
```

### Update All Dependencies

```bash
pixi update
```

## Updating the Lock File

After modifying `pyproject.toml` manually, regenerate the lock file:

```bash
pixi lock
```

## Common Development Tasks

### Run Tests

```bash
pixi run -e dev pytest
```

### Check Code Style

```bash
pixi run -e dev black --check hippunfold/
pixi run -e dev isort --check hippunfold/
```

### Format Code

```bash
pixi run -e dev black hippunfold/
pixi run -e dev isort hippunfold/
```

### Use Python Interactively

```bash
pixi run -e dev python
```

## Troubleshooting

### Environment is stale or corrupted

Remove and reinstall:

```bash
rm -rf .pixi/
pixi install
```

### A package won't install

Check if it's available on conda-forge or bioconda:

```bash
pixi search package-name
```

If not available, you may need to use a PyPI dependency instead. Add it to `[tool.pixi.pypi-dependencies]` in `pyproject.toml`.

### Lock file conflicts

If `pixi.lock` has merge conflicts, regenerate it:

```bash
git checkout HEAD -- pixi.lock  # discard local changes
pixi lock                        # regenerate from pyproject.toml
```

## Key Files

- **`pyproject.toml`**: Defines project metadata, dependencies, and tasks
- **`pixi.lock`**: Auto-generated lock file with exact package versions (commit to git)
- **`.pixi/`**: Local directory where environments are installed (ignore from git)

## Summary

```bash
# First time setup
pixi install

# Every day - interactive development
pixi shell -e dev
# ... run commands ...
exit

# Or run individual commands
pixi run -e dev command-name

# Check and fix code quality
pixi run -e dev quality_check
pixi run -e dev quality_fix
```

Happy developing! 🦆
