from setuptools import setup, find_packages


def get_requires():
    requires = []
    with open("requirements.txt", "rt") as req_file:
        for line in req_file:
            requires.append(line.rstrip())
    return requires


setup(name='eva_cttv_pipeline',
      version='2.7.2',
      packages=find_packages(),
      install_requires=get_requires(),
      #! TBD: list as a dependency subpackage 'clinvar_xml_utils.clinvar_xml_utils.clinvar_xml_utils'
      package_data={
          'eva_cttv_pipeline': ['OT_SCHEMA_VERSION']
      },
      tests_require=get_requires(),
      setup_requires=get_requires(),
      test_suite='tests'
      )
