def save_file(raw, to_path):
    """Write (blueprint) data to file."""

    file = open(to_path, 'w')

    for l in raw:
        file.writelines([l])
    file.close()