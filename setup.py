from setuptools import setup, find_packages

setup(name='variantdetective',
      version='1.0.0',
      description='variantdetective',
      url='https://github.com/OLF-Bioinformatics/VariantDetective',
      author='Phil Charron', 
      author_email='phil.charron@inspection.gc.ca', 
      license='MIT', 
      install_requires=['numpy','pandas'], 
      package_dir={"": "variantdetective"},
      packages=find_packages(where="variantdetective"),
      zip_safe=False) 
