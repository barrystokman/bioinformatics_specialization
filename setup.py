from setuptools import setup


def parse_reqs(req_path='./requirements.txt'):
    """Recursively parse requirements from nested pip files."""
    install_requires = []
    with open(req_path, 'r') as handle:
        # remove comments and empty lines
        lines = (line.strip() for line in handle
                 if line.strip() and not line.startswith('#'))
        for line in lines:
            # check for nested requirements files
            if line.startswith('-r'):
                # recursively call this function
                install_requires += parse_reqs(req_path=line[3:])
            else:
                # add the line as a new requirement
                install_requires.append(line)
    return install_requires


setup(name='bioinfo_package',
      version='1.0',
      py_modules=['cli.biocli_module'],
      install_requires=parse_reqs(),
      entry_points={'console_scripts': ['bioinfo=cli.biocli_module:biocli']},
      )
