import numpy as np
import math
import random
import datetime
import sys, os
from copy import copy
import multiprocessing as mp
import chiripa as chi

# #########################################################################
def Zwork_per_core(list_of_parameters):

    """
    Calculate a single Z using the list of parameters. This is used in parallel calculation of Z.
    This is called by ``_Z_calc_group`` function


    Args:
        list_of_parameters (list): List of the parameters.


    .. note:: Items in list of parameters (by position):

        * ``filenamecoords``: (list, string, 2-element) A 2-element lists containing the path to the coordinate files for each segment
        * ``filenametop``: (list, string, 2-element) A 2-element lists containing the path to the topology files for each segment
        * ``dirZ``: (string) Directory to store the pairwise structures generated during the process if ``Z_debug`` is ``True``
        * ``Z_samples`` (integer): Number of samples to average Z
        * ``Z_putTrialsMonomer`` (integer): Maximum number of trials to put a Segment
        * ``Z_nonbondedMethod`` (string): Non-bonded method
        * ``Z_Debug`` (Boolean). If ``True`` saves information in ``dirZ``
        * ``seed`` (integer): Seed for the random number generator.


    Returns:
        (tuple): tuple containing:

            * ``filenamecoords[0]`` (str) Name of the file containing the coordinatates of the Base segment.
            * ``filenamecoords[1]`` (str) Name of the file containing the coordinatates of the Screen segment.
            * ``Z_avg`` (float) Averaged Z-value.
            * ``Z_std`` (float) Standard deviation of the Z-value.
            * ``list_of_parameters[2]`` (str) Label "11", "12", "21" or "22".
            * ``Z_hist`` (list) Returns the values of all calculated Z which can be uesed to build an histogram.

    """

    name = list_of_parameters[0]
    filenamecoords = list_of_parameters[1]
    filenametop = list_of_parameters[2]
    filetypeatoms = list_of_parameters[3]
    dirZ = list_of_parameters[5]
    Z_samples = list_of_parameters[6]
    Z_putTrialsMonomer = list_of_parameters[7]
    Z_nonbondedMethod = list_of_parameters[8]
    Z_debug = list_of_parameters[9]
    iseed = list_of_parameters[10]

    # Segments to calculate the Z-number
    s0 = chi.Segment(name[0] ,filecoord=filenamecoords[0], filetop=filenametop[0], filetypeatoms=filetypeatoms[0])
    s1 = chi.Segment(name[1] ,filecoord=filenamecoords[1], filetop=filenametop[1], filetypeatoms=filetypeatoms[1])

    # args = task_list
    Z_avg, Z_std, Z_hist = calculate_average_z_number(s0, s1, samples_Z=Z_samples,
                                putTrialsMonomer=Z_putTrialsMonomer, idir=dirZ,
                                box=chi.CubicBox(40.0, 40.0, 40.0),
                                debug=Z_debug,
                                iseed = iseed,
                                nonbondedEvaluation=Z_nonbondedMethod)

    del s0
    del s1
    return filenamecoords[0], filenamecoords[1], Z_avg, Z_std, list_of_parameters[4], Z_hist

# #########################################################################
def Z_calc_group(filename_list, filecoords_list, filetop_list, filetypes_list, logger, Z_samples, Z_putTrialsMonomer,
                 Z_nonbondedMethod, Z_debug, workers=4, iseed=None):

    """
    Group data to run in parallel the values of Z for each pair ("11", "12", "21", "22").


    Args:
        filecoords_list (list): A 2-element lists containing the path to the coordinate files for each segment.
        filetop_list (list): A 2-element lists containing the path to the topology files for each segment.
        logger (logger): Log results.
        Z_samples (int): int Number of samples to average Z.
        Z_putTrialsMonomer (int): Maximum number of trials to put a Segment.
        Z_nonbondedMethod (str): Non-bonded method.
        Z_Debug (bool): if ``True`` saves information in ``dirZ``.
        workers (int): Number of processes.
        iseed (int): Seed for the random number generator.


    Returns:
        results (list): A list containing the results of the Z calculation in parallel


    **Example of results**

    .. code-block::

        0 = {tuple: 6} ('../data/n-hexane.pdb', '../data/n-hexane.pdb', 10.6, 0.7999999999999999, '11', [10, 12, 10, 11, 10])
        1 = {tuple: 6} ('../data/n-hexane.pdb', '../data/nitrobenzene.pdb', 9.2, 0.9797958971132713, '12', [8, 9, 9, 11, 9])
        2 = {tuple: 6} ('../data/nitrobenzene.pdb', '../data/n-hexane.pdb', 9.8, 0.7483314773547883, '21', [9, 9, 10, 10, 11])
        3 = {tuple: 6} ('../data/nitrobenzene.pdb', '../data/nitrobenzene.pdb', 9.4, 1.0198039027185568, '22', [11, 9, 9, 10, 8])

    """

    # *****************************
    PROCESSES_AVAILABLE = int(mp.cpu_count())
    NUMBER_OF_TASKS = workers

    start_time = datetime.datetime.now()
    fmt = "%a %d/%m/%Y, %H:%M:%S"

    with open("Z_running_tmp.txt",'w') as f:
        fmt = "%a %d/%m/%Y, %H:%M:%S"
        l = "Z calculation starting at {} using {:d} cores.\n".format(start_time.strftime(fmt), workers)
        f.writelines(l)

    # Write parameters  ===================================
    message  = ""
    message += "\t Starting time                                         = {} \n".format(start_time.strftime(fmt))
    message += "\t Number of samples to calculate coordination number(Z) = {} \n".format(Z_samples)
    message += "\t Number of trials to put a Segment                     = {} \n".format(Z_putTrialsMonomer)
    message += "\t Non-bonded method (VdW radii)                         = {} \n".format(Z_nonbondedMethod)
    message += "\t Debug flag                                            = {} \n".format(Z_debug)
    message += "\t Number of processes available                         = {} \n".format(PROCESSES_AVAILABLE)
    message += "\t Number of tasks in parallel                           = {} \n".format(NUMBER_OF_TASKS)
    print(message) if logger is None else logger.info(message)

    # Starting mp.pool ====================================
    pool_tasks = mp.Pool(processes=NUMBER_OF_TASKS)

    # Create data =========================================
    list_name = []
    for i in filename_list:
        for j in filename_list:
            list_name.append([i, j])

    list_coords = []
    for i in filecoords_list:
        for j in filecoords_list:
            list_coords.append([i, j])

    list_top = []
    for i in filecoords_list:
        for j in filecoords_list:
            list_top.append([i, j])

    list_types = []
    for i in filetypes_list:
        for j in filetypes_list:
            list_types.append([i, j])

    data = [ [list_name[0], list_coords[0], list_top[0], list_types[0], "11", "./aZ11-coordination", Z_samples,
              Z_putTrialsMonomer, Z_nonbondedMethod, Z_debug, iseed] ,
             [list_name[1], list_coords[1], list_top[1], list_types[1], "12", "./aZ12-coordination", Z_samples,
              Z_putTrialsMonomer, Z_nonbondedMethod, Z_debug, iseed] ,
             [list_name[2], list_coords[2], list_top[2], list_types[2], "21", "./aZ21-coordination", Z_samples,
              Z_putTrialsMonomer, Z_nonbondedMethod, Z_debug, iseed] ,
             [list_name[3], list_coords[3], list_top[3], list_types[3], "22", "./aZ22-coordination", Z_samples,
              Z_putTrialsMonomer, Z_nonbondedMethod, Z_debug, iseed] ]

    results = pool_tasks.map(Zwork_per_core, (data[i] for i in range(4)))

    pool_tasks.close()
    pool_tasks.join()

    message = "\n"
    for i in range(workers):
        name1 = results[i][0].split("/")[-1].split(".")[0]
        name2 = results[i][1].split("/")[-1].split(".")[0]
        label = results[i][4]
        line = "\t#Z# Z_{0:s} ({1:15s} / {2:15s}) = {3:.2f} +- {4:.2f}\n".format(label, name1, name2,
                                                           results[i][2],results[i][3])
        message += line

    end_time = datetime.datetime.now()
    elapsed_time = end_time - start_time
    message += "\n\t#Z# TIME: Z calculation in {:.1f} seconds using {:d} cores.\n".format(elapsed_time.total_seconds(), workers)

    with open("Z_done.txt",'w') as f:
        os.remove("Z_running_tmp.txt")
        l = "Z calculation done in {:.1f} seconds using {:d} cores.\n".format(elapsed_time.total_seconds(), workers)
        f.writelines(l)

    print(message) if logger is None else logger.info(message)

    return results

# #########################################################################
def write_Z_results(Z_results):

    with open("Z_results.log", 'w') as f:
        for item in Z_results:
            name1 = item[0].split("/")[-1].split(".")[0]
            name2 = item[1].split("/")[-1].split(".")[0]
            label = item[4]
            line = "\t Z_{0:s} ({1:15s} / {2:15s}) = {3:.2f} +- {4:.2f}\n".format(label, name1, name2,
                                                                                  item[2],item[3])
            f.writelines(line)

    with open("Z_histogram.log",'w') as f:
        max_freq = -np.infty
        min_freq = np.infty
        for item in Z_results:
            label = item[4]
            zlist = item[5]
            end = max(zlist)
            start = min(zlist)
            hist = np.histogram(zlist, bins =end-start+1)
            if max(hist[0]) > max_freq: max_freq = max(hist[0])
            if min(hist[0]) < min_freq: min_freq = min(hist[0])
            i = 0
            f.writelines("# Z_{0:s}\n".format(label))
            for ibin in range(start, end+1):
                f.writelines("{0:d} {1:d}\n".format(ibin, hist[0][i]))
                i+= 1
            f.writelines("\n\n")
            #print(hist[0], start, end, label)

    create_gnu_histogram(max_freq, min_freq)

# #########################################################################
def create_gnu_histogram(max_freq, min_freq):

    lines= """
set term wxt 1 enhanced dashed size 800,800 font "Arial,9"
set multiplot layout 2,2
set encoding iso_8859_1
set style line 1 lt 1 ps 1.4 lc rgb "black"  pt 6 lw 2.0
set style line 2 lt 1 ps 1.4 lc rgb "blue"  pt 8 lw 2.0
set style line 3 lt 1 ps 1.4 lc rgb "red"  pt 10 lw 2.0

f1="./Z_histogram.log"

set style data histogram
set style histogram cluster gap 2
set style fill solid
set boxwidth 1.0
set xtics rotate by 0
set yrange [%d:%d]
set grid

p "Z_histogram.log" i 0 u 2:xtic(1) title "Pair 1-1" linecolor "red"

p "Z_histogram.log" i 1 u 2:xtic(1) title "Pair 1-2" linecolor "blue"

p "Z_histogram.log" i 2 u 2:xtic(1) title "Pair 2-1" linecolor "green"

p "Z_histogram.log" i 3 u 2:xtic(1) title "Pair 2-2" linecolor "brown"

unset multiplot
"""%(min_freq, max_freq)

    with open("histogram.gnu", "w") as f:
        f.writelines(lines)

# #########################################################################
def atoms_with_minimum_distances(s):
    """
    Calculate the minimun distance between two set of atoms and return a list
    with the following format [iat_set1, iat_set2, distance]

   Args:
        s (Simulator): Simulator to perform the calculations

    """

    # Check if the instance is a simulator
    if not isinstance(s, Simulator):
        print("Warning: {0:} is not a Simulator object".format(type(s)))
        print("Warning: gen_conf_fan_algorithm function needs a Simulator object")
        print("Warning: Nothing is done, returning None")
        return None


    imol1 = s.get_ith_molecule(0)
    imol2 = s.get_ith_molecule(1)

    ref = imol1.get_coords()
    conf = imol2.get_coords()
    dij, rijx, rijy, rijz = distance_array(ref, conf)

    min_dist = np.amin(dij)
    print("min_dist, np.where (dij == min_dist)")

    print(min_dist, np.where (dij == min_dist))
    print("====")

# #########################################################################
def calc_single_anisotropy(central_mon, neigh_mon, nrotations=100, nconf=2000, box= chi.CubicBox(30.0, 30.0, 30.0),
                           iseed=None, debug=False, nonbondedEvaluation="truhlar", idir=".", itry=0 ):


    if debug:
        if os.path.isfile("{:s}/anisotropy_full.pdb".format(idir)):
            os.remove("{:s}/anisotropy_full.pdb".format(idir))

    # Seed for the random generator ========
    if iseed is None:
        iseed = int(982628732*random.random())

    sim = chi.Simulator(mol=(central_mon, neigh_mon), box=box)

    cog_list = list()
    for iconf in range(nconf):
        gen_conf_fan_algorithm_c(sim, iseed=iseed, nonbondedEvaluation=nonbondedEvaluation)
        cog_list.append(sim.get_ith_molecule(-1).center_of_geom())
        iseed += 1
        if debug:
            sim.write_PDB_simulator(filenamePDB="{:s}/anisotropy_full.pdb".format(idir),box=sim.box,
                                    with_vmd_tcl=True, appendpdb=True)

    S_avg, S_std = chi.testPointInsidePyramid(sim.box.Lx, cog_list, debug=debug, iseed=iseed, nrotations=nrotations)

    if debug:
        p = "{:s}/acog_debug.xyz".format(idir, itry)
        with open(p, 'w') as f:
            l = len(cog_list)
            f.writelines("{}\n".format(str(l)))
            f.writelines("Mol\n")
            for item in cog_list:
                l = "C {0:f} {1:f} {2:f}\n".format(item[0], item[1], item[2])
                f.writelines(l)

    del sim
    return S_avg, S_std

# #########################################################################
def gen_conf_fan_algorithm(s, iseed=None, chosen_solution="MAX"):

    """
    The function accepts a simulator object and then change the coordinates of the last molecule
    in the molecule list in accordance
    with the previous algorithm.

    Args:
        s (Simulator): Simulator to perform the calculations
        iseed (int): Seed for the random generator number
        chosen_solution (str): default="MAX", allowed values=("MAX", "MIN", "RANDOM"). Blanco et al. shown that\
        the maximum and minimum solution can avoid the overlap between molecules. Fan et al. used the MAX
        mol_move (int). Number of the molecule to move. The followed order is that\
        in the mol_list of the simulator object

    Returns:
        r_max: (float)

    """

    # Check if the instance is a simulator
    if not isinstance(s, chi.Simulator):
        print("Warning: {0:} is not a Simulator object".format(type(s)))
        print("Warning: gen_conf_fan_algorithm function needs a Simulator object")
        print("Warning: Nothing is done, returning None")
        return None

        # Generate euler angles ========
    if iseed is None:
        random.seed()
    else:
        random.seed(iseed)

    # Translate the center of geometry of the first and last monomers to the origin (0,0,0)
    # s.write_PDB_simulator(filenamePDB="./sim01.pdb",
    #                       title="Frame 01. Original",with_vmd_tcl=False, appendpdb=False)
    nmols = s.get_nmolecules()
    for i in range(0, nmols, nmols-1):
        imol = s.get_ith_molecule(i)
        cog = imol.center_of_geom()
        imol.translate_vector(-cog)

    # s.write_PDB_simulator(filenamePDB="./sim01.pdb",
    #                       title="Frame 02. Translate",with_vmd_tcl=True, appendpdb=True)

    # New orientation for the molecule 2 in the simulator object
    imol1 = s.get_ith_molecule(0)
    imol2 = s.get_ith_molecule(-1)
    imol2.euler_orientation(iseed=iseed)
    # s.write_PDB_simulator(filenamePDB="./sim01.pdb",
    #                       title="Frame 03. Orientation",with_vmd_tcl=True, appendpdb=True)

    # Generate a vector on a sphere
    n = point_on_unit_sphere(iseed=iseed)

    # Calculate all distances between the atoms of mol1(i) and atoms of mol2(j)
    # Equation (10) --> Fan et al. Macromolecules, 25(14), 1992, 3667-3676
    ref = imol1.get_coords()
    conf = imol2.get_coords()

    elements1 = imol1._elements
    elements2 = imol2._elements
    dij, rijx, rijy, rijz = chi.distance_array(ref, conf)
    sol_max = list()
    sol_min = list()
    index_mask = list()
    # pij, rijx, rijy, rijz = distance_array_purepython(ref, conf)
    # print(np.testing.assert_almost_equal(dij, pij))
    #for i in range(imol1._natoms):
    for i in range(imol1._natoms):
        if i in [imol1._dummy_head_atom, imol1._dummy_tail_atom]:
            Rvdw_i = chi.element_vdw_truhlar_radius['C']+0.4
        else:
            Rvdw_i = chi.element_vdw_truhlar_radius[elements1[i]]
        for j in range(imol2._natoms):
            rij = np.array([rijx[i,j], rijy[i,j], rijz[i,j]])
            #np.testing.assert_almost_equal(dij[i,j], np.linalg.norm(rij))
            b =  np.dot(n,rij)
            Rvdw_j = chi.element_vdw_truhlar_radius[elements2[j]]
            discriminat = b*b - dij[i,j]*dij[i,j] + (Rvdw_i + Rvdw_j)**2
            #print(i, j, elements1[i], elements2[j], Rvdw_i, Rvdw_j)
            # For real solutions get both solutions and store an global index in
            # index_mask. index_mask[0] = [0,2] : Global index 0 for the i-th 0 (local number in mol1)
            #                                     and    index 2 for the j-th 2 (local number in mol2)
            if discriminat >= 0.0:
                s1 = -b + math.sqrt(discriminat)
                s2 = -b - math.sqrt(discriminat)
                sol_max.append(s1)
                sol_min.append(s2)
                index_mask.append([i, j])
                #print(i, j, dij[i,j], Rvdw_i, Rvdw_j, imol1.get_coords()[i], imol2.get_coords()[j],   discriminat, -b + math.sqrt(discriminat), -b - math.sqrt(discriminat))
                #print(i, j, discriminat, -b + math.sqrt(discriminat), -b - math.sqrt(discriminat))

    # Get the maximum (positive direction) and minimum (negative direction)
    # distance as well as the indexes of the atoms.
    r_max = max(sol_max)
    i_max = index_mask[sol_max.index(max(sol_max))]
    r_min = min(sol_min)
    i_min = index_mask[sol_min.index(min(sol_min))]

    # New coordinates -->
    if chosen_solution.upper() == "MAX":
        r = r_max
    elif chosen_solution.upper() == "MIN":
        r = r_min
    elif chosen_solution.upper() == "RANDOM":
        p = random.random()
        if p < 0.5:
            r = r_min
        else:
            r = r_max
    else:
        print("Warning {:s} is not allowed in generate_pair_conformations function. "
              "Using MAX as default.".format(chosen_solution))
        r = r_max

    # Move molecule 2 alogn the vector r*n
    imol2.translate_vector(np.dot(r,n))
    # s.write_PDB_simulator(filenamePDB="./sim01.pdb",
    #                       title="Frame 04. Move m2",with_vmd_tcl=True, appendpdb=True)

    return r_max, r_min, i_max, i_min


# #########################################################################
def gen_conf_fan_algorithm_c(s, nonbondedEvaluation, iseed=None,
                             chosen_solution="MAX"):

    """
    Solving the overlap problem using the algortihms exposed in :

    Molecular silverware. I. General solutions to excluded volume constrained problems
    (Blanco, M. Journal of Computational Chemistry, 1991, 12(2), 237-247

    Application of Molecular Simulation To Derive Phase Diagrams of Binary Mixtures
    (Fan C. et al. Macromolecules, 1992, 25(14), 3667-3676)

    The function accepts a simulator object and then change the coordinates of the last molecule
    in the molecule list in accordance
    with the previous algorithm.

    ``Parameters``:
        * **s** (type: Simulator): Simulator to perform the calculations
        * **iseed** (type int, default=None): Seed for the random generator number
        * **chosen_solution** (type: string, default="MAX", allowed values=("MAX", "MIN", "RANDOM"). Blanco et al. shown that\
        the maximum and minimum solution can avoid the overlap between molecules. Fan et al. used the MAX
        * **mol_move** (type int, default=-1). Number of the molecule to move. The followed order is that\
        in the mol_list of the simulator object

    ``Returns``:
        * **r_max**

    """

    # Check if the instance is a simulator
    if not isinstance(s, chi.Simulator):
        print("Warning: {0:} is not a Simulator object".format(type(s)))
        print("Warning: gen_conf_fan_algorithm function needs a Simulator object")
        print("Warning: Nothing is done, returning None")
        return None

        # Generate euler angles ========
    if iseed is None:
        random.seed()
    else:
        random.seed(iseed)

    # Translate the center of geometry of the first and last monomers to the origin (0,0,0)
    # s.write_PDB_simulator(filenamePDB="./sim01.pdb",
    #                       title="Frame 01. Original",with_vmd_tcl=False, appendpdb=False)
    nmols = s.get_nmolecules()
    for i in range(0, nmols, nmols-1):
        imol = s.get_ith_molecule(i)
        cog = imol.center_of_geom()
        imol.translate_vector(-cog)

    # s.write_PDB_simulator(filenamePDB="./sim01.pdb",
    #                       title="Frame 02. Translate",with_vmd_tcl=True, appendpdb=True)

    # New orientation for the molecule 2 in the simulator object
    imol1 = s.get_ith_molecule(0)
    imol2 = s.get_ith_molecule(-1)
    imol2.euler_orientation(iseed=iseed)
    # s.write_PDB_simulator(filenamePDB="./sim01.pdb",
    #                       title="Frame 03. Orientation",with_vmd_tcl=True, appendpdb=True)

    # Generate a vector on a sphere
    n = point_on_unit_sphere(iseed=iseed)

    # Calculate all distances between the atoms of mol1(i) and atoms of mol2(j)
    # Equation (10) --> Fan et al. Macromolecules, 25(14), 1992, 3667-3676
    ref = imol1.get_coords()
    conf = imol2.get_coords()
    elements1 = imol1._elements
    elements2 = imol2._elements

    Rvdw_i = []
    Rvdw_j = []
    if nonbondedEvaluation.upper() == "TRUHLAR":
        for i in range(imol1._natoms):
            if i in [imol1._dummy_head_atom, imol1._dummy_tail_atom]:
                Rvdw_i.append(chi.element_vdw_truhlar_radius['C']+0.4)
            else:
                Rvdw_i.append(chi.element_vdw_truhlar_radius[elements1[i]])
        for j in range(imol2._natoms):
             Rvdw_j.append(chi.element_vdw_truhlar_radius[elements2[j]])

        r_max, r_min, i_max, i_min = chi.calc_discriminant(ref, conf, n, imol1._natoms, imol2._natoms, Rvdw_i, Rvdw_j)

    elif nonbondedEvaluation.upper() == "OKUWAKI_CORRECTION":
        ref_typeelements = imol1._typeelements
        conf_typeelements = imol2._typeelements
        dij = np.zeros((imol1._natoms, imol2._natoms))

        for i in range(imol1._natoms):
            itype = ref_typeelements[i]
            for j in range(imol2._natoms):
                jtype = conf_typeelements[j]
                str1 = itype+"-"+jtype
                str2 = jtype+"-"+itype
                try:
                    dij[i,j] = chi.distances_revaluated_Okuwaki[str1]
                except KeyError:
                    try:
                        dij[i,j] = chi.distances_revaluated_Okuwaki[str2]
                    except KeyError:
                        print("========= ERROR ============")
                        print("Gen_conf_fan_algorithm_c cannot be used.")
                        print("Method: {:} is not defined".format(nonbondedEvaluation))
                        print("Available methods are: truhlar and okuwaki_correction")
                        print("========= ERROR ============")
                        sys.exit()

        r_max, r_min, i_max, i_min = chi.calc_discriminant_okuwaki(ref, conf, n, imol1._natoms, imol2._natoms, dij)

    else:
        print("========= ERROR ============")
        print("Gen_conf_fan_algorithm_c cannot be used.")
        print("Method: {:} is not defined".format(nonbondedEvaluation))
        print("Available methods are: truhlar and okuwaki_correction")
        print("========= ERROR ============")
        sys.exit()

    # Get the maximum (positive direction) and minimum (negative direction)
    # distance as well as the indexes of the atoms.

    # New coordinates -->
    if chosen_solution.upper() == "MAX":
        r = r_max
    elif chosen_solution.upper() == "MIN":
        r = r_min
    elif chosen_solution.upper() == "RANDOM":
        p = random.random()
        if p < 0.5:
            r = r_min
        else:
            r = r_max
    else:
        print("Warning {:s} is not allowed in generate_pair_conformations function. "
              "Using MAX as default.".format(chosen_solution))
        r = r_max

    # Move molecule 2 along the vector r*n
    imol2.translate_vector(np.dot(r,n))
    # s.write_PDB_simulator(filenamePDB="./sim01.pdb",
    #                       title="Frame 04. Move m2",with_vmd_tcl=True, appendpdb=True)
    return r_max, r_min, i_max, i_min


# #########################################################################
def point_on_unit_sphere(iseed=None):

    """
    "Understanding Molecular Simulation", Frenkel-Smith, 2nEdition
    Algorithm 42 (Random Vector on a Unit Sphere), Pag 578

    A random vector on a sphere of diameter 1 (origin 0,0,0)

    :param iseed
    :return: Point on a unit sphere from xcenter
    """

    if iseed is None:
        random.seed()
    else:
        random.seed(iseed)

    ransq = 2.0
    ran1 = 1. - 2. * random.random()
    ran2 = 1. - 2. * random.random()
    while ransq >= 1:
        ran1 = 1. - 2. * random.random()
        ran2 = 1. - 2. * random.random()
        ransq = ran1*ran1 + ran2*ran2

    ranh=2.*math.sqrt(1.-ransq)

    bx=ran1*ranh
    by=ran2*ranh
    bz=(1-2.*ransq)

    point = np.array([bx, by, bz])

    return point


# #########################################################################
def calculate_single_z_number(central_mon, neigh_mon, ntry=2000,
                              iseed=None, box=chi.CubicBox(30.0, 30.0, 30.0),
                              nonbondedEvaluation="truhlar"):

    """
    This function calculate the coordination number Z using the algorithm
    proposed by Fan et al. Macromolecules, 25(14), 1992, 3667-3676.

    The Z number is calculated using a nearest neighouurs around a center molecule.
    The algorithm is as follows:

        1. Generate  the Simulator object with the central monomer and the
        first position for the other monomer
        2. Try to put the next neighbour around the central monomer
        2.1. Get the neighbour position around the central monomer (gen_conf_fan_algorithm)
        2.2. Check overlap
        2.3. if overlap = False. z=z+1 and go to step 2
        2.4. if overlap = True & itry < ntry --> Go to step 2.1
        2.5. if overlap = True & itry > ntry --> Go to step 3
        3. There is no more room for a new position

    Parameters:

        central_mon: (type: Segment): Central monomer
        neigh_mon: (type: Segment): Neighbour monomer
        ntry: (type: integer, default: 2000): Number of trials to put the neighbour monomer around the central monomer
        iseed: (type: integer, default:None): Seed for the random number generator (RNG). The value\
         ``None`` is given to get a seed value from random.seed() which uses the current\
        time to generate it. Otherwise, use iseed as intial value for the seed RNG.

    Returns:

        z: (type: float) Number of coordination
        sim: (type: Simulator) Simulator containing the coordination configuration

    """

    if iseed is None:
        iseed = int(982628732*random.random())

    sim = chi.Simulator(mol=(central_mon, neigh_mon), box= box)
    chi.gen_conf_fan_algorithm_c(sim, nonbondedEvaluation, iseed=iseed)

    seed = iseed + 1

    z_coord = 1
    itry = 0
    while True:
        # Get a copy of the last monomer
        m_new = copy(sim.get_ith_molecule(-1))
        # Add the monomer to the simulator
        sim.add_molecule(m_new)
        _, _, i_max, _ = chi.gen_conf_fan_algorithm_c(sim, nonbondedEvaluation, iseed=seed)
        if i_max[0] in [central_mon._dummy_head_atom, central_mon._dummy_tail_atom]:
            sim.remove_molecule(-1)
            seed += 1
            continue
        else:
            isoverlapping = check_overlap_z_number(sim, nonbondedEvaluation=nonbondedEvaluation)

        if isoverlapping:
            sim.remove_molecule(-1)
            seed += 1
            itry += 1
            if itry>ntry: break
            continue
        else:
            itry = 0
            z_coord +=1
            if z_coord > 100:
                print("========= WARNING ============")
                print("Coordination number (Z) >100 ")
                print("Likely something is going wrong")
                print("Please, check calculate_single_z_number function.")
                print("========= WARNING ============")
                break
        seed += 1

    return z_coord, sim


# #########################################################################
def calculate_average_z_number(c_mon, n_mon, samples_Z=20,
                               putTrialsMonomer=2000, iseed=None, debug=False,
                               idir=".", box=chi.CubicBox(30.0, 30.0, 30.0),
                               nonbondedEvaluation="truhlar"):

    # Copy of the segments to be sure that they are different instances
    # The function works with the copies.
    central_mon = copy(c_mon)
    neigh_mon = copy(n_mon)

    # Check if the monomers contain type atom information when
    # using the reevaluated distances proposed by Okuwaki et al.
    if nonbondedEvaluation.upper() == "OKUWAKI_CORRECTION":
        typeelements1 = central_mon._typeelements
        typeelements2 = neigh_mon._typeelements
        settype1 = None
        settype2 = None
        if typeelements1 is not None: settype1 = "Not None"
        if typeelements2 is not None: settype2 = "Not None"
        if typeelements1 is None or typeelements2 is None:
            print("========= ERROR ============")
            print("Overlap cannot be evaluated.")
            print("Method: {:} requires information of type atoms.".format(nonbondedEvaluation))
            print("Set type atoms through the set_typeatoms function of the Segment class")
            print("Central Segment: {:} --> {:}".format(central_mon._filecoord, settype1))
            print("Neigh   Segment: {:} --> {:}".format(neigh_mon._filecoord, settype2))
            print("========= ERROR ============")
            sys.exit()

    # Print each 10% the time
    # stride = int(samples_Z / 10)
    # if stride == 0: stride = 1
    # print("\t{:d} trials for coordination number {} and {}.\n\tPrinting info each {} trials".
    #       format(samples_Z, central_mon._filecoord, neigh_mon._filecoord, stride))

    z_list = []
    start_time = datetime.datetime.now()
    for itry in range(samples_Z):
        tmp, s = calculate_single_z_number(central_mon, neigh_mon,
                                           ntry=putTrialsMonomer, iseed = iseed, box=box,
                                           nonbondedEvaluation=nonbondedEvaluation)
        z_list.append(tmp)

        # if itry == 0 or (itry+1) % stride == 0.0 or itry == samples_Z-1:
        #     end_time = datetime.datetime.now()
        #     elapsed_time = end_time - start_time
        #     print("\t {} of {} trials in {:.1f} seconds. PID: {}".format(itry+1,samples_Z, elapsed_time.total_seconds(), os.getpid()))
        #     print("{} {}\n".format(central_mon._filecoord, neigh_mon._filecoord))

        if debug:
            path_Z_PDB = "{:s}/coordinationXX2b_{:04d}.pdb".format(idir, itry)
            s.write_PDB_simulator(filenamePDB=path_Z_PDB, box=s.box, with_vmd_tcl=True, dir=idir)

    Z_avg = np.mean(z_list)
    Z_std = np.std(z_list)
    a = min(z_list)
    b = max(z_list)
    #Z_hist  = np.histogram(z_list, bins=b-a+1, range=[a,b], density=False)
    return Z_avg, Z_std , z_list


# #########################################################################
def check_overlap_z_number(s, nonbondedEvaluation="truhlar"):

    overlap = False

    # The current molecule to be put in the coordination sphere
    imol_current = s.get_ith_molecule(-1)
    ref_coord = imol_current._coords
    ref_elements = imol_current._elements
    iatoms = imol_current._natoms
    if nonbondedEvaluation.upper() == "OKUWAKI_CORRECTION":
        ref_typeelements = imol_current._typeelements

    # Molecules just placed in the coordinate sphere
    nmols_to_check = s.get_nmolecules() - 1
    jatoms = 0
    conf_coord = np.ndarray((0,3))
    conf_elements = np.ndarray(0)
    conf_typeelements = np.ndarray(0)
    for j in range(1,nmols_to_check):
        jmol_placed =  s.get_ith_molecule(j)
        tmp_coord = jmol_placed._coords
        tmp_elements = jmol_placed._elements
        conf_coord = np.concatenate((conf_coord,tmp_coord),axis=0)
        conf_elements = np.concatenate((conf_elements,tmp_elements),axis=0)
        jatoms += jmol_placed._natoms
        tmp_typeelements = jmol_placed._typeelements
        if nonbondedEvaluation.upper() == "OKUWAKI_CORRECTION":
            conf_typeelements = np.concatenate((conf_typeelements,tmp_typeelements),axis=0)

    # All distances between the atoms in the molecule to be put and the rest of atoms
    # already placed.
    # TODO: This can be optimized due to is not neccesary to calculate all distances
    # in the case of overlapping structures.
    dij, _, _, _ = chi.distance_array(ref_coord, conf_coord)

    if nonbondedEvaluation.upper() == "TRUHLAR":
        for i in range(iatoms):
            Rvdw_i = chi.element_vdw_truhlar_radius[ref_elements[i]]
            for j in range(jatoms):
                Rvdw_j = chi.element_vdw_truhlar_radius[conf_elements[j]]
                if dij[i,j] < Rvdw_i+Rvdw_j:
                    overlap = True

    elif nonbondedEvaluation.upper() == "OKUWAKI_CORRECTION":

        for i in range(iatoms):
            itype = ref_typeelements[i]
            for j in range(jatoms):
                jtype = conf_typeelements[j]
                str1 = itype+"-"+jtype
                str2 = jtype+"-"+itype
                try:
                    d = chi.distances_revaluated_Okuwaki[str1]
                except KeyError:
                    try:
                        d = chi.distances_revaluated_Okuwaki[str2]
                    except KeyError:
                        print("========= ERROR ============")
                        print("Overlap cannot be evaluated.")
                        print("Method: {:} is not defined".format(nonbondedEvaluation))
                        print("Available methods are: truhlar and okuwaki_correction")
                        print("========= ERROR ============")
                        sys.exit()

                if dij[i,j] < d:
                    overlap = True

    else:
        print("========= ERROR ============")
        print("Overlap cannot be evaluated.")
        print("Method: {:} is not defined".format(nonbondedEvaluation))
        print("Available methods are: truhlar and okuwaki_correction")
        print("========= ERROR ============")
        sys.exit()

    return overlap


# #########################################################################
def calculate_average_anisotropy(c_mon, n_mon, orientations=100, nconf=10000,
                        box= chi.CubicBox(30.0, 30.0, 30.0), idir=".",
                        iseed=None, debug=False, nonbondedEvaluation="truhlar" ):


    # Seed for the random generator ========
    if iseed is None:
        iseed = int(982628732*random.random())

    # Copy of the segments to be sure that they are different instances
    # The function works with the copies.
    central_mon = copy(c_mon)
    neigh_mon = copy(n_mon)

    # Print each 25% the time
    stride = int(orientations / 4)
    if stride == 0: stride = 1
    print("\t{:d} orientations (conf: {:d}) for anisotropy factor of {} and {}.\n\tPrinting info each {} trials".
          format(orientations, nconf, central_mon._filecoord, neigh_mon._filecoord, stride))

    S_array = np.zeros(orientations)
    start_time = datetime.datetime.now()
    for i in range(orientations):
        central_mon.euler_orientation(iseed=iseed)
        S_avg, S_std = calc_single_anisotropy(central_mon, neigh_mon, nconf=nconf,
                                   iseed=iseed, debug=debug,
                                   nonbondedEvaluation=nonbondedEvaluation, itry=i)
        S_array[i] = S_avg

        if i == 0 or (i+1) % stride == 0.0 or i == orientations-1:
            end_time = datetime.datetime.now()
            elapsed_time = end_time - start_time
            print("\t {} of {} trials in {:.1f} seconds".format(i+1,orientations,
                                                                elapsed_time.total_seconds()))


    return S_array.mean(), S_array.std()

