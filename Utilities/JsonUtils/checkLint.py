import json


def check(filename):
    try:
        with open(filename) as f:
            lf = json.load(f)
            return True
    except:
        return False


if __name__ == "__main__":
    # Allow for command line call
    check(sys.argv[1])
