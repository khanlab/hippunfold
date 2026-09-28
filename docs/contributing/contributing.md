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

## Release process

- Create and publish a GitHub Release tag in the format `vX.Y.Z`, with optional PEP 440 suffixes (for example `vX.Y.Zrc1`).
- Publishing the release triggers both publish workflows:
  - Conda package publication to prefix.dev
  - Docker image publication to both dockerhub and ghcr.io
- Package versioning is now tag-driven:
  - Python package version is derived dynamically from git tags at build time.
  - Conda package version is injected from the release tag in CI before publish.
- Docker images are tagged with both `vX.Y.Z` and `X.Y.Z`; `latest` is only pushed for non-prereleases.


