from setuptools import setup, find_packages

setup(name='variantdetective',
      version='1.0.1',
      description='VariantDetective: an accurate all-in-one pipeline for detecting consensus bacterial SNPs and SVs',
      url='https://github.com/OLF-Bioinformatics/VariantDetective',
      author='Phil Charron', 
      author_email='phil.charron@inspection.gc.ca', 
      license='MIT', 
      #install_requires=['pandas>=2.1.3', 'numpy>=1.26', 'tensorflow>=2.8.0'], 
      packages=find_packages(),
      include_package_data=True,
      package_data={
          'variantdetective': ['clair3_models/*/*'],
      },
      entry_points={
          'console_scripts': [
              'variantdetective=variantdetective.main:main',
          ],
      }
)
