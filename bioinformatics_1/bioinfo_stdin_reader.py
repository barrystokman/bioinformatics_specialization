"""
This module reads files from STDIN and returns a read file object
"""
# TODO: DON'T READ FROM STDIN, READ FROM CLI DATASET ARGUMENT


import sys


class ReadFile:
    """
    Creates a class instance of a file read with STDIN
    The type of file is derived from the number of line that the file contains
    The second element in the value tuple in the class variable 'file_types' indicates which line
    contains the input data
    """
    def __init__(self):
        """
        initialize Readfile object by reading from STDIN
        """
        self.lines = sys.stdin.read().splitlines()

    def __str__(self):
        """
        returns a readable print output
        """
        return_str = ""
        for line in self.lines:
            return_str += line + "\n"
        return return_str

    def _count_lines(self):
        """
        returns the number of lines in the read file. Used to determine the type of file
        """
        return len(self.lines)

    def get_full_text(self):
        """
        return the full text of a file
        """
        return self.lines


def main():
    read_file = ReadFile()
    print(read_file)
    pass


if __name__ == '__main__':
    main()
