# Contributing to HippUnfold

To contribute code or documentation to HippUnfold:

```bash
git clone https://github.com/khanlab/hippunfold.git
cd hippunfold
pixi install
```

If you want the development tools as well, install the development environment:

```bash
pixi install --environment dev
```

## Deep learning / nnU-net model files

HippUnfold downloads nnU-net model files on demand and stores them in the cache directory by default:

```bash
~/.cache/hippunfold/
```

You can override this location with:

```bash
export HIPPUNFOLD_CACHE_DIR=/path/to/custom/cache
```
