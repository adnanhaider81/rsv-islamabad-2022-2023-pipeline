# Contributing

Thanks for considering a contribution.

- Open an issue to discuss substantial changes before a PR.
- Follow the style and structure already used in the Snakefile and config.
- Keep test data small. If you need larger data use Git LFS or generate on the fly.
- Add or update tests under `tests/` and make sure CI passes.

## Development setup
```bash
conda env create -f env/environment.yml
conda activate rsv-env
snakemake -s workflow/Snakefile -n --configfile tests/config_test.yaml
```
