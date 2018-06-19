from setuptools import setup


setup(name='bioinfo_package',
      version='1.0',
      py_modules=['cli.biocli_module'],
      install_requires=['Click'],
      entry_points={'console_scripts': ['bioinfo=cli.biocli_module:biocli']},
      )
