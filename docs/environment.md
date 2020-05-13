# Setting up the common environment

Log in to the LSF cluster, where all data processing must take place. Switch to a common EVA production user instead of your personal account. Then adjust and execute the commands below. They will set up the environment, fetch and build the code.

Notes:
* The first five variables are installation-specific and are blanked in this repository. You can get the values for the EVA use case from the [private repository](https://github.com/EBIvariation/configuration/blob/master/open-targets-configuration.md).
* By modifying the `*REMOTE` and `*BRANCH` variables, you can run arbitrary versions of both the main and the VEP pipeline. This is highly useful for development and debugging. By default it fetches master branches of both repositories.
* Running these commands will overwrite any local changes you had in the repository copy on the cluster, including any changes to the Java ClinVar XML parser. Be sure to commit and push those before re-running this block of commands.

```bash
# This variable should point to the directory where the clone of this repository is located on the cluster
export CODE_ROOT=

# Location of Python installation which you configured using build instructions
export PYTHON_INSTALL_PATH=

# Location of bcftools installation path
export BCFTOOLS_INSTALL_PATH=

# The directory where subdirectories for each batch will be created
export BATCH_ROOT_BASE=

# Base path of FTP directory on the cluster
export FTP_PATH_BASE=

# Base bsub command line for all commands.
export BSUB_CMDLINE="bsub"

# Setting up Python paths
export PATH=${PYTHON_INSTALL_PATH}:${PYTHON_INSTALL_PATH}/bin:${BCFTOOLS_INSTALL_PATH}:$PATH
export PYTHONPATH=${PYTHON_INSTALL_PATH}

# External service paths
CLINVAR_PATH_BASE="ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar"

export MAIN_REMOTE=origin
export MAIN_BRANCH=master
export VEP_REMOTE=origin
export VEP_BRANCH=master

cd ${CODE_ROOT}
git fetch ${MAIN_REMOTE}
git checkout ${MAIN_BRANCH}
git reset --hard ${MAIN_REMOTE}/${MAIN_BRANCH}

git submodule update --init --recursive
cd vep-mapping-pipeline
git fetch ${VEP_REMOTE}
git checkout ${VEP_BRANCH}
git reset --hard ${VEP_REMOTE}/${VEP_BRANCH}
cd ..

python3 -m pip -q install -r requirements.txt
python3 -m pip -q install -r vep-mapping-pipeline/requirements.txt
python3 setup.py install
```
