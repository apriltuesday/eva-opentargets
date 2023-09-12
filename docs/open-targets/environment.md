# Setting up the common environment

1. Log in to the LSF cluster (currently `codon`), where all data processing must take place.
1. Using a `become` command, switch to a common EVA production user instead of your personal account.
1. Adjust and execute the commands below to set up the environment. Notes:
    - The first five variables are installation-specific and are blanked in this repository. You can get the values for the EVA installation from the [private repository](https://github.com/EBIvariation/configuration/blob/master/open-targets-configuration.md).
    - These instructions use the Gitlab deployment of the pipeline onto the cluster. To use another installation just modify `CODE_ROOT` and `PYTHON_BIN` accordingly.

```bash
# This variable should point to the directory where the clone of this repository is located on the cluster
export CODE_ROOT=

# Location of Python installation which you configured using build instructions
export PYTHON_INSTALL_PATH=

# Location of Nextflow installation path
export NEXTFLOW_INSTALL_PATH=

# The directory where subdirectories for each batch will be created
export BATCH_ROOT_BASE=

# Base path of FTP directory on the cluster
export FTP_PATH_BASE=

# Setting up paths
export PATH=${PYTHON_INSTALL_PATH}:${NEXTFLOW_INSTALL_PATH}:$PATH
export PYTHONPATH=${PYTHON_INSTALL_PATH}

# Location of Python executable, pointing to the virtualenv
export PYTHON_BIN=${CODE_ROOT}/env/bin/python
```
