from setuptools import setup


DESCRIPTION = "Umap and Bismap: tools for genome " +\
    "and methylome mappability"


def readme():
    with open('README.rst') as f:
        return f.read()


setup(name='umap',
      version='0.1',
      description=DESCRIPTION,
      long_description=readme(),
      classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics"],
      url='https://bitbucket.org/hoffmanlab/ubismap',
      author='Mehran Karimzadeh, Anshul Kundaje',
      author_email='mehran.karimzadeh@uhnresearch.ca',
      license='MIT',
      packages=['umap'],
      install_requires=[
          "argparse",
          "numpy",
          "pandas"],
      include_package_data=True,
      zip_safe=False)
