## Preparing environment for building documentation

The detailed instructions to install Conda, which is an open-source package management and environment management system, for Linux environment can be found in [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).

The following commands can be used to create Conda environment and build documentation.

```
# updates conda
conda update -n base -c defaults conda
# creates new environment to build documentation
conda create --name docs
# activates newly created environment for documentation
conda activate docs
# deactivates environment
conda deactivate
# install required modules such as sphinx and bitex extension
conda install sphinx
conda install -c conda-forge sphinxcontrib-bibtex
```

## Building documentation

The documentation can be found under **docs/** directory. To build the documentation, `make html` command needs to be run under this directory. If all the required [Sphinx](https://www.sphinx-doc.org/en/master/) pacages found and installed correctly, the compiled documentation can be found under `build/html/` directory and you could just open `index.html` with browser.
