from setuptools import setup, find_packages

setup(name='NISEpy',
      version='1.0',
      description='a library to enable the quick calculation of patch potential forces on curved surfaces from (KPFM) images of the surface potential',
      url='https://github.com/NTNU-SmallSat-Lab/Imaging-Pipeline',
      author='Joe Garrett',
      author_email='jgarret2@umd.edu',
      license='BSD License 2.0',
      packages=find_packages(),
      install_requires = ['numpy ', 'scipy', 'igor',
            'mpmath', 'tqdm', 'tables','matplotlib',
            'scikit-image'],
      zip_safe=False)
