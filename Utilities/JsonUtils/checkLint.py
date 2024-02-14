import json
import sys


def check(filename):
    try:
        with open(filename) as fid:
            json.load(fid)
            return True
    except:
        return False


if __name__ == "__main__":
    # Allow for command line call
    check(sys.argv[1])
