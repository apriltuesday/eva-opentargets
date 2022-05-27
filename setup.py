import os
import os.path
from setuptools import setup, find_packages


def get_package_data():
    package_data = []
    for root, dirnames, filenames in os.walk('./eva_cttv_pipeline/evidence_string_generation/resources'):
        root = root.replace("./eva_cttv_pipeline/", "")
        for filename in filenames:
            new_fn = os.path.join(root, filename)
            package_data.append(new_fn)
    return package_data


def get_requires():
    requires = []
    with open("requirements.txt", "rt") as req_file:
        for line in req_file:
            requires.append(line.rstrip())
    return requires


setup(name='eva_cttv_pipeline',
      version='2.5.4',
      packages=find_packages(),
      install_requires=get_requires(),
      #! TBD: list as a dependency subpackage 'clinvar_xml_utils.clinvar_xml_utils.clinvar_xml_utils'
      package_data={
          'eva_cttv_pipeline': get_package_data()
      },
      tests_require=get_requires(),
      setup_requires=get_requires(),
      test_suite='tests'
      )
