from setuptools import setup, find_packages

setup(name='variantdetective',
      version='1.0.0',

description='variantdetective',
      url='https://github.com/OLF-Bioinformatics/VariantDetective',
      author='Phil Charron', 
      author_email='phil.charron@inspection.gc.ca', 
      license='MIT', 
      packages=find_packages(), 
      install_requires=['numpy','pandas'], 
      zip_safe=False) 
