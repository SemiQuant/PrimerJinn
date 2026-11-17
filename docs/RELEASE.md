# Releasing primerJinn

Everything below is prepared; nothing runs until a tag is pushed. Order: GitHub, then PyPI, then conda.

## 0. Before tagging

1. Set the version in `primerJinn/__init__.py` (`__version__`) and in `conda/meta.yaml` (`version`). It is the single source for the wheel, the GUI footer and the release check.
2. Move the "unreleased" heading in `CHANGELOG.md` to the version and date.
3. Run the suite and build locally:

```bash
python -m pytest tests -q
python -m build                      # dist/primerjinn-<version>.tar.gz and .whl
python -m twine check dist/*
```

4. Try the wheel in a fresh environment:

```bash
python -m venv /tmp/pj && /tmp/pj/bin/pip install "dist/primerjinn-<version>-py3-none-any.whl[gui]"
/tmp/pj/bin/primerjinn-multiplex --genome example/ref.fasta --positions example/mtb_positions.tsv --out /tmp/pj/mtb --check
```

5. Commit on `main` (merge the feature branch first).

## 1. GitHub

```bash
git tag -a v2.1.0 -m "primerJinn 2.1.0"
git push origin main --tags
```

The `release` workflow (`.github/workflows/release.yml`) builds the sdist and wheel, checks that the tag equals `__version__`, creates the GitHub release with the files attached and generated notes, and publishes to PyPI.

## 2. PyPI

The workflow publishes with trusted publishing, so no token is kept anywhere. One-time setup before the first tag, on https://pypi.org/manage/account/publishing/ (the `primerJinn` project already exists there, version 2.0.3):

- Owner `SemiQuant`, repository `PrimerJinn`, workflow `release.yml`, environment `pypi`.
- In the GitHub repository, Settings > Environments, create `pypi` (optionally with a required reviewer, which turns the publish job into a manual approval).

Manual fallback, from the same `dist/` as above:

```bash
python -m twine upload dist/*
```

## 3. Conda (bioconda)

primerJinn is not on conda yet. The recipe is `conda/meta.yaml`. After the PyPI release:

1. Fill in the sdist checksum: `sha256sum dist/primerjinn-<version>.tar.gz` (or from the PyPI file page).
2. Fork https://github.com/bioconda/bioconda-recipes, copy the recipe to `recipes/primerjinn/meta.yaml`, and open a pull request. Bioconda's CI builds it (`noarch: python`); every run dependency (`viennarna`, `pysam`, `primer3-py`, `numba`, `blast`, `streamlit`, `altair`) is already on bioconda or conda-forge.
3. Local check, optional: `conda build conda/ -c conda-forge -c bioconda` (needs `conda-build`), or `bioconda-utils lint` from the recipes checkout.

Once merged, `conda install -c conda-forge -c bioconda primerjinn` works and `environment.yml` can be reduced to that line.

## Later versions

Bump `__version__` and `conda/meta.yaml`, update the changelog, tag `vX.Y.Z`, push. For conda, bioconda's bot opens the version bump automatically once the recipe exists; only the sha256 and the dependency list need checking.
