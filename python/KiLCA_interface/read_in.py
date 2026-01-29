def read_in(path):
    """Reads in text file and returns list of the lines."""

    raw = {}
    with open(path) as f:
        content = f.readlines()

    return content
