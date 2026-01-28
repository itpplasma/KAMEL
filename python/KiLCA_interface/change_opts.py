def change_opts(raw, ind, dat, sep):
    """Change options of blueprint file contained in raw."""

    if len(ind) != len(dat):
        raise ValueError("Unequal length of index and new data.")

    count = 0
    for i in ind:
        c = raw[i].split(sep, 1)  # 1 to split only after first occurrence of "#"
        c[0] = opt2str(dat[count]) + " "
        raw[i] = c[0] + "     # " + c[1]
        count = count + 1

    return raw


def change_opts_fieldpy(raw, ind, dat, sep):
    """Change options of blueprint file contained in raw."""

    if len(ind) != len(dat):
        raise ValueError("Unequal length of index and new data.")

    count = 0
    for i in ind:
        c = raw[i].split(sep, 1)  # 1 to split only after first occurrence of "#"
        c[0] = opt2str(dat[count]) + " "
        raw[i] = c[0] + "       " + c[1]
        count = count + 1

    return raw


def opt2str(opt):
    """Change option to string."""

    if isinstance(opt, str):
        return opt
    if isinstance(opt, float):
        return "{:.2e}".format(opt)
    if isinstance(opt, int):
        return str(opt)

    if isinstance(opt, list):
        # for complex numbers:
        if len(opt) == 2:
            # if len(opt[0]) == 1:
            return "(" + "{:.2e}".format(opt[0]) + ", " + "{:.2e}".format(opt[1]) + ")"
            if len(opt[0]) != 1:
                l = []
                for i in range(0, len(opt[0])):
                    pass
