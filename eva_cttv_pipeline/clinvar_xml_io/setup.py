from setuptools import setup, find_packages

classifiers = [
    "Natural Language :: English",
    "License :: OSI Approved :: Apache Software License",
    "Programming Language :: Python :: 3"
]


def get_requires(filename):
    """ Function used to read the required dependencies (e.g. in 'requirements.txt')    
    """
    requires = []
    with open(filename, "rt") as req_file:
        for line in req_file:
            requires.append(line.rstrip())
    return requires


# Extract the markdown description. Supported by PyPi in its native format
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='clinvar_xml_io',
      version='0.0.1',
      author_email="opentargets-clinvar@ebi.ac.uk",
      url="https://github.com/EBIvariation/eva-opentargets/eva_cttv_pipeline/clinvar_xml_io",
      packages=find_packages(),
      install_requires=get_requires("requirements.txt"),
      # List additional files that are included in the package (e.g. tests, datafiles...)
      package_data={
          'clinvar_xml_io': ["data/*"]
      },
      package_dir={'': '.'},
      description="Utilities and classes to parse the ClinVar XML",
      long_description=long_description,
      long_description_content_type="text/markdown",
      tests_require=get_requires("test_requirements.txt"),  # ! Review the required packages for the tests
      setup_requires=get_requires("test_requirements.txt"),
      test_suite='tests',
      # More classifiers to be added (https://pypi.org/classifiers/ - e.g. Operating System)
      classifiers=classifiers
      )
