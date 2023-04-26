from setuptools import setup, find_packages


def get_requires():
    requires = []
    with open("requirements.txt", "rt") as req_file:
        for line in req_file:
            requires.append(line.rstrip())
    return requires


setup(name='cmat',
      version='3.0.0-dev',
      packages=find_packages(),
      install_requires=get_requires(),
      package_data={
          'cmat': ['OT_SCHEMA_VERSION']
      },
      tests_require=get_requires(),
      setup_requires=get_requires(),
      test_suite='tests'
      )
