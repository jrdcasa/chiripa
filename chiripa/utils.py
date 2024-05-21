from __future__ import print_function
from sys import getsizeof, stderr
from itertools import chain
from collections import deque, defaultdict
from reprlib import repr
import numpy as np
import os
import shutil
from zipfile import ZipFile
import logging
import colorlog
import sqlite3
import pandas as pd
import chiripa as chi


# ########################################################################################
def total_size(o, handlers={}, verbose=False):
    """ Returns the approximate memory footprint an object and all of its contents.

    Automatically finds the contents of the following builtin containers and
    their subclasses:  tuple, list, deque, dict, set and frozenset.
    To search other containers, add handlers to iterate over their contents:

        handlers = {SomeContainerClass: iter,
                    OtherContainerClass: OtherContainerClass.get_elements}
                    
    https://code.activestate.com/recipes/577504/

    """
    dict_handler = lambda d: chain.from_iterable(d.items())
    all_handlers = {tuple: iter,
                    list: iter,
                    deque: iter,
                    dict: dict_handler,
                    set: iter,
                    frozenset: iter,
                   }
    all_handlers.update(handlers)     # user handlers take precedence
    seen = set()                      # track which object id's have already been seen
    default_size = getsizeof(0)       # estimate sizeof object without __sizeof__

    def sizeof(obj):
        if id(obj) in seen:       # do not double count the same object
            return 0
        seen.add(id(obj))
        s = getsizeof(obj, default_size)

        if verbose:
            print(s, type(obj), repr(obj), file=stderr)

        for typ, handler in all_handlers.items():
            if isinstance(obj, typ):
                s += sum(map(sizeof, handler(obj)))
                break
        return s

    return sizeof(o)

# ########################################################################################
def test_total_size():
    """
    Test the total_size function with some data
     
    """

    print ("=========== Dictionary ==============")
    d = dict(a=1, b=2, c=3, d=[4,5,6,7], e='a string of chars')
    # print(total_size(d, verbose=False))
    # print(getsizeof(d))
    print("Size calculated with total size function      : {} bytes".format(total_size(d, verbose=False)))
    print("Size calculated with getsizeof native function: {} bytes".format(getsizeof(d)))

    print ("=========== List ==============")
    l = [10,100, 200, 300]
    print("Size calculated with total size function      : {} bytes".format(total_size(l, verbose=False)))
    print("Size calculated with getsizeof native function: {} bytes".format(getsizeof(l)))

    print ("=========== Set ==============")
    l = (10,100, 200, 300)
    print("Size calculated with total size function      : {} bytes".format(total_size(l, verbose=False)))
    print("Size calculated with getsizeof native function: {} bytes".format(getsizeof(l)))

    print ("=========== File ==============")
    f = open("utils.py", 'r')
    print("Size calculated with total size function      : {} bytes".format(total_size(f, verbose=False)))
    print("Size calculated with getsizeof native function: {} bytes".format(getsizeof(f)))

    print ("=========== Numpy ==============")
    arr = np.zeros((100,100), dtype=float)
    print("Size calculated with total size function      : {} bytes".format(total_size(arr, verbose=False)))
    print("Size calculated with getsizeof native function: {} bytes".format(getsizeof(arr)))

    print ("=========== String ==============")
    s = "El Quijote de la Mancha"
    print("Size calculated with total size function      : {} bytes".format(total_size(s, verbose=False)))
    print("Size calculated with getsizeof native function: {} bytes".format(getsizeof(s)))

# ########################################################################################
def compress_files(localdir="./", dirbase="./", zipname="gaussian_inputs.zip"):

    # Empty path file paths list
    file_paths = []
    prev_dir = os.getcwd()+"/"
    #os.chdir(dirbase)
    #dirbase=prev_dir

    if dirbase[-1] != "/":
        dirbase += "/"

    # Walking through directory and subidrectories
    for root, directories, files in os.walk(prev_dir):
        for filename in files:
            # join the two strings in order to form the full filepath.
            if filename.find("zip") != -1 or filename.find("full_send") !=-1: continue
            filepath = os.path.join(root, filename)
            file_paths.append("./"+filepath.replace(prev_dir,''))

    # writing files to a zipfile
    with ZipFile(zipname, 'w') as zip:
        #writing each file
        for file in file_paths:
            zip.write(file)

    # Remove the folder of the input files
    delete_folder(dirbase)

# ########################################################################################
def delete_folder(folder):

    shutil.rmtree(folder, ignore_errors=False)

# #######################################################################################
def energy_convergence_MC(dfl, boxl, T=300, seed=9837724721, nseeds=200):


    energy_1 = dfl[0]['energy'].values[0]  # energy_monomer_1 (hartrees)
    energy_2 = dfl[1]['energy'].values[0]  # energy_monomer_2 (hartrees)

    # Energy monomer1, momomer2, monomer1+monomer1(1_1), monomer1+monomer2(1_2), monomer2+monomer1(2_1), monomer2+monomer2(2_2)
    energy_ref = [energy_1, energy_2, energy_1 * 2, energy_1 + energy_2, energy_2 + energy_1, energy_2 * 2]

     # MONTECARLO STATISTICS =========================================
    iseed = seed
    eij_mc = {2:[], 3:[], 4:[], 5:[]}
    samp_mc = {2:[], 3:[], 4:[], 5:[]}
    sij_mc = {2:[], 3:[], 4:[], 5:[]}
    wij = []
    wij_12 = []

    for i in range(0, nseeds):

        print ("{} of {} seeds".format(i, nseeds))
        E_dictaccepted = {2: [], 3: [], 4: [], 5: []}

        for ipair in range(2, 6):
            E_current = 0.0; E_sum = 0.0; accepted = 0; E_listgenerated = [];  cog_list = []
            for index, row in dfl[ipair].iterrows():
                E_listgenerated.append(row['delta_energy'])
                isaccepted, RT = chi.metropolis_MC(E_current, row['delta_energy'], T, iseed=iseed)
                iseed += 1
                if isaccepted:
                    E_current = row['delta_energy']
                    E_dictaccepted[ipair].append(E_current)
                    cog_list.append([float(s) for s in row.cog.split()])

            eij = np.mean(E_dictaccepted[ipair])
            eij_mc[ipair].append(eij)
            samp_mc[ipair].append(len(E_dictaccepted[ipair]))
            sij, sij_std = chi.testPointInsidePyramid(boxl, cog_list)
            sij_mc[ipair].append(sij)

    return eij_mc, samp_mc, sij_mc




