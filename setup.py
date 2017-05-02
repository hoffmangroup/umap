from setuptools import setup


DESCRIPTION = "Umap and Bismap: tools for genome " +\
    "and methylome mappability"


def readme():
    with open('README.rst') as f:
        return f.read()


setup(name='umap',
      version='1.1.0',
      description=DESCRIPTION,
      long_description=readme(),
      classifiers=[
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: GNU General Public License v3 ",
        "(GPLv3)"],
      url='https://bitbucket.org/hoffmanlab/proj/bismap',
      author="Mehran Karimzadeh, Carl Ernst, " +
      "Anshul Kundaje, Michael M. Hoffman",
      author_email='mehran.karimzadeh@uhnresearch.ca',
      license='GPLv3',
      packages=['umap'],
      install_requires=[
          "argparse",
          "numpy",
          "pandas"],
      include_package_data=True,
      zip_safe=False)
