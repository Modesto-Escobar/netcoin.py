from setuptools import setup

setup(name='netcoin',
      version='1.0.0',
      description='Interactive Analytic Networks',
      long_description='Create interactive analytic networks.',
      classifiers=[
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Visualization',
      ],
      keywords='network javascript d3',
      url='https://sociocav.usal.es/blog/nca/',
      author='Modesto Escobar, David Barrios, Carlos Prieto, Luis Martinez-Uribe',
      author_email='metal@usal.es',
      license='GPLv2+',
      packages=['netcoin'],
      install_requires=[
          'pandas',
          'numpy',
          'scipy',
      ],
      python_requires='>=3.8',
      include_package_data=True,
      zip_safe=False)
