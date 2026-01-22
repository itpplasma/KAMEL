import h5py
import numpy as np

mrks = ["s", "d", "x", "."]
cols = [
    "tab:blue",
    "tab:orange",
    "tab:green",
    "tab:red",
    "tab:purple",
    "tab:brown",
    "tab:pink",
    "tab:gray",
    "tab:olive",
    "tab:cyan",
]


def determine_phase_of_phi_WKB_struct(dat: dict):
    integral = np.zeros(len(dat["r_found_only_one"]), dtype=complex)
    intrange = np.where(dat["r_found_only_one"] > dat["r_res"])[0][0]
    print(intrange)
    for i, r in enumerate(dat['r_found_only_one']):
        if r > dat['r_res']:
            integral[i] = -np.trapezoid(dat['k_r1_only_one'][intrange:i], dat['r_found_only_one'][intrange:i])
            #integral[i] = np.trapezoid(dat['k_r1_only_one'][intrange:i], dat['r_found_only_one'][intrange:i])
        if r < dat['r_res']:
            integral[i] = -np.trapezoid(dat['k_r1_only_one'][i:intrange], dat['r_found_only_one'][i:intrange])
    phase = np.exp(1j * integral)
    return phase


def determine_phase_of_phi_WKB(r, kr, r_res):
    integral = np.zeros(len(r), dtype=complex)
    intrange = np.where(r > r_res)[0]
    for i, r_val in enumerate(r):
        if r_val > r_res:
            integral[i] = np.trapezoid(kr[intrange:i], r[intrange:i])
        if r_val < r_res:
            integral[i] = -np.trapezoid(kr[i:intrange], r[i:intrange])
    phase = np.exp(1j * integral)
    return phase


def only_one_root(dat: dict):
    dat["k_r1_only_one"] = []
    dat["k_r2_only_one"] = []
    dat["r_found_only_one"] = []
    for i, r in enumerate(dat["r_found"]):
        if not r == dat["r_found"][i - 1]:
            dat["r_found_only_one"].append(r)
            dat["k_r1_only_one"].append(
                np.min(np.abs(np.real(dat["k_r1"][np.where(dat["r_found"] == dat["r_found"][i])])))
                - 1j
                * np.min(
                    np.abs(np.imag(dat["k_r1"][np.where(dat["r_found"] == dat["r_found"][i])]))
                )
                # +1j*np.min(np.imag(dat['k_r1'][np.where(dat['r_found'] == dat['r_found'][i])]))
            )
            dat["k_r2_only_one"].append(
                np.min(np.abs(np.real(dat["k_r2"][np.where(dat["r_found"] == dat["r_found"][i])])))
                + 1j
                * np.min(
                    np.abs(np.imag(dat["k_r2"][np.where(dat["r_found"] == dat["r_found"][i])]))
                )
            )
    dat["r_found_only_one"] = np.array(dat["r_found_only_one"])
    dat["k_r1_only_one"] = np.array(dat["k_r1_only_one"], dtype=complex)
    dat["k_r2_only_one"] = np.array(dat["k_r2_only_one"], dtype=complex)
    return dat


def read_dat(file: str):
    h5f = h5py.File(file, "r")
    dat = {}
    dat["r_found"] = np.array(h5f["r_found"][:])
    dat["k_r1"] = np.array(h5f["k_r1"][:])
    dat["k_r2"] = np.array(h5f["k_r2"][:])
    dat["options"] = {}
    for option in h5f["options"].keys():
        dat["options"][option] = h5f["options/" + option][()]
    dat["r_res"] = np.array(h5f["r_res"])
    h5f.close()
    return dat
