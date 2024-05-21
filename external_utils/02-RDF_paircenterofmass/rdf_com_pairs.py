#!/usr/bin/python
import argparse
import os
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt

"""
 rdf_com_pairs.py
 
 This program calculates the RDF (Radial Distribution Function) of a pair of molecules. This pair is generated 
 in Chiripa. The program needs MDAnalysis library (https://www.mdanalysis.org)
 
 The program has been only checked with python3.
 
 usage: rdf_com_pairs -f gaussian_inputs_1_1.xyz -set1 0:19 -set2 20:39

"""


def print_header():

    """
    Print a header of the program
    """

    print("#############################################")
    print("#              RDF COM Pairs                #")
    print("#############################################")


def check_arguments():

    """
    Check CLI arguments

    Returns:
        A tuple:
            Name of the xyz file containing the trajectory (filename)
            String with indexes for the first COM (set1)
            String with indexes for the xeconf COM (set2)

    """

    isok = True
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", "-f", help="XYZ trajectory file")
    parser.add_argument("--set1", "-s1", help="String for indexes, i.e: \"1:3\" or \"1,2,3\"")
    parser.add_argument("--set2", "-s2", help="String for indexes, i.e: \"1:3\" or \"1,2,3\"")

    args = parser.parse_args()

    print_header()

    if args.file:
        filename = args.file
    else:
        print("ERROR --> Please give a file containing a XYZ trajectory.\n")
        isok = False

    if args.set1:
        set1 = args.set1
    else:
        print("ERROR --> Please give a valid index-range i.e: 0:19 or 1,2,3,4 for set1.")
        isok = False

    if args.set2:
        set2 = args.set2
    else:
        print("ERROR --> Please give a valid index-range i.e: 0:19 or 1,2,3,4 for set2.")
        isok = False

    if not isok:
        parser.print_help()
        exit(0)

    return filename, set1, set2


def get_all_indexes_mol(inputstring):

    """
    Get an inputstring from the command line and returns a list with all atom indices.

    Example:
        "1:4, 15, 20:22" --> [1,2,3,4,15,20,21,22]

    Args:
        inputstring: A command input string. Examples =  "1:4, 15, 20:22", "1:4", "1, 2, 3 ,8"

    Returns:
        A list of atoms

    """

    latoms = []
    ss = inputstring.split(",")
    for item in ss:
        try:
            at0, at1 = item.split(":")
            for i in range(int(at0), int(at1)+1):
                latoms.append(i)
        except ValueError:
            try:
                latoms.append(int(item))
            except ValueError:
                return None

    return latoms


def xyz_to_gro(fname, atmol, atch):

    """
    XYZ trajectory file to GRO trajectory file
    Args:
        fname: Name of the XYZ file
        atmol: A list os list containing the index atom for each molecule,
            i.e [[0,1,2],[3,4,5,6]] represents two molecules with three and four atoms, respectively.

    Returns:

    """

    with open(fname) as f:
        fnamegro = os.path.basename(fname).split(".")[0]+".gro"
        fout = open(fnamegro, 'w')
        iframe = 0

        while True:

            # Read number of atoms from xyz
            natoms = f.readline()

            # End condition
            if not natoms:
                break
            natoms = int(natoms)

            # Read the title in the xyz trajectory
            title = "Frame {}\n".format(iframe)

            # Write title and number of atoms in the gro file
            fout.writelines(title)
            fout.writelines("{0:10d}\n".format(natoms))

            # Skip title line in the XYZ trajectory
            nolines = f.readline()

            # Read all the coordinates for the current frame
            lines_coords = list()
            for iatom in range(natoms):
                lines_coords.append(f.readline()[:-1])

            # Write coordinates in the gro file
            iat = 0
            for item in lines_coords:
                atom_name, x, y, z = item.split()
                vx = 0.0
                resname="R{0:03d}".format(atch[iat])
                l = "{0:>5d}{1:>5s}{2:>5s}{3:>5d}{4:>8.3f}{5:>8.3f}{6:>8.3f}\n".format(atch[iat]+1, resname, atom_name, iat+1,
                                                            float(x)/10, float(y)/10, float(z)/10)
                fout.writelines(l)
                iat += 1
            l = "{0:>8.5f} {1:>8.5f} {2:>8.5f}\n".format(3.0, 3.0, 3.0)
            fout.writelines(l)
            iframe += 1


def run_program():

    """
    Docs
    Returns:

    """

    # Filename
    filename, range1, range2 = check_arguments()
    print("Filename         : {}".format(filename))
    print("Atoms to COM1    : {}".format(range1))
    print("Atoms to COM2    : {}".format(range2))

    # From xyz to gro
    # Getting the chain identification for each atom
    atmol = []
    lat = get_all_indexes_mol(range1)
    atmol.append(lat)
    lat = get_all_indexes_mol(range2)
    atmol.append(lat)
    ich = 0
    atch = []
    for imol_list in atmol:
        for iat in imol_list:
            atch.append(ich)
        ich += 1

    xyz_to_gro(filename, atmol, atch)

    # Load MDA trajectory
    distances = list()
    u = mda.Universe(filename)
    nframes = u.trajectory.n_frames
    print("Number of frames : {}".format(nframes))
    set1 = u.select_atoms('index {}'.format(range1))
    set2 = u.select_atoms('index {}'.format(range2))
    for ts in u.trajectory:
        com1 = set1.center_of_mass()
        com2 = set2.center_of_mass()
        d = np.linalg.norm(com2 - com1)
        distances.append(d)

    print(distances)
    maxdist = np.ceil(np.max(distances))
    mindist = np.floor(np.min(distances))
    delta = 0.2 #angstroms
    bins = int((maxdist-0)/delta) + 1
    # # Take a histogram
    counts, lenghts = np.histogram(distances, bins=bins)

    # print(counts)
    # print(lenghts)
    # # Normalize counts by averaging across the Nframes
    # #counts = counts / nframes
    # # calculate the volume of each spherical shell:
    # shell_volumes = 4/3*np.pi*(lenghts[1:]**3-lenghts[:-1]**3)
    # # normalize by volume of each shell
    # counts = counts / shell_volumes
    #
    # # Plot
    plt.plot((lenghts[:-1]+lenghts[1:])/2.0, counts, '-o')

    plt.show()

    print("Job Done!!!!!")


if __name__ == '__main__':

    run_program()

