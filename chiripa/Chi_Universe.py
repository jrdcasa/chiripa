import numpy as np
import os
import datetime
import glob
import sqlite3
import multiprocessing as mp
import chiripa as chi
import pandas as pd
from collections import defaultdict

class Chi_Universe(object):

    """Create an universe to perform the calculations to get the |chi| parameter

    .. |chi| replace:: :math:`{\chi}`

    Args:
        * ``input_dict`` (dict): A dictionary of the required keywords to perform the calculation.
        * ``header`` (boolean, default=True): If True then prints the program header on the log file.
        * ``output_screen`` (boolean, default=False): If True then prints the log file on screen.

    The chi_universe creates the environment to calculate the |chi| parameter
    from the structure of the segments. A **input_dict** must be given in the signature.
    This dictionary codifies all information needed to perform calculations
    The allowed **key:values** pairs in this dictionary
    are the following:

        * **'names'** (required): A list of strings for the molecule names.
        * **'filecoords'** (required): A list of strings for filenames to read coordinates.
        * **'filetop'** (required): A list of strings for filenames to read topology.
        * **'coordination_numbers_Z'** (boolean, optional): If not present the coordination numbers (Z\ :sub:`ij`) are not calculated.
        * **'Z-parameters'** (dictionary, optional) : A dictionary defining the features of the Z-calculation.
        * **'calculate_volume'** (boolean, optional) : If True it calculates the monomer volume.

    .. warning::
        The lists **'names'**, **'filecoords'**, **'filetop'** in the dictionary must have the same length. The order is also important,
        for example, the position 0 must correspond to the info of the Segment 1 and the position
        1 to the Segment 2.

    .. attention::
        Not all keys are defined.

    Properties:
        * self._names (list of string): Name of the segments
        * self._segments (list of segments): List of the two segments.
        * self._logger (Logger): Log to output the files.

    .. attention::
        Not all properties are defined.

    """
    # ========================================================================
    def __init__(self, input_dict, header=True, logname="output.log", output_screen=False, clean=None):

        fmt = "%a %d/%m/%Y, %H:%M:%S"

        # Clean files and exit ====================================================================
        #   full : Delete of files including the database
        #   light: Delete of files except the database
        self._clean = clean
        if not self._clean is None:
            self._clean_files()

        # List of names ============================================================================
        self._names = []

        # Create a list with segments (molecules) ==================================================
        self._segments = []

        # Create a dictionaty for results ==========================================================
        self._dictresults = dict()

        # Logger for CHIRIPA =======================================================================
        # If Z_results.log file does not exist and it is not running or the  database does not exists then
        # chiripa assumes that this is the first time that we launch
        # the script, thus the sampling_run is deactivated. If QM calculations are running or finished,
        # two paths can be taken: i) perform sampling of the conformacions or check the runnning calculations
        dblist = glob.glob("*.db")
        if (not os.path.exists("Z_results.log") and not os.path.exists("Z_running_tmp.txt"))\
            or len(dblist) == 0:
            self._logger = chi.init_logger("Output", fileoutput=logname, append=False, inscreen=output_screen)
            if header:
                self.__print_header()
            input_dict['sampling_run'] = False
            msg1 = "Sampling run has been deactivated!!!!."
        elif input_dict['sampling_run'] and os.path.exists("Done"):
            self._logger = chi.init_logger("Output", fileoutput = "output_results.log", append=False, inscreen=output_screen)
            if header:
                self.__print_header()
        else:
            self._logger = chi.init_logger("Output", fileoutput = "output_check.log", append=False, inscreen=output_screen)
            if header:
                self.__print_header()
        # Print ===================================================================================
        self._logger.info("\n\tJob Information\n\t---------------")
        self._logger.info("\tStarting: \t {}\n".format(datetime.datetime.now().strftime(fmt)))

        # Parse the input dictionary ==============================================================
        self._parse(input_dict)

        if input_dict['sampling_run'] and os.path.exists("Done"):
            #  Collect energies in the database and perform sampling
            #self._collect_results_sequential()
            s = chi.Sampling(self._database_name, self._Samplingmode, self._Samplingtemperatures,
                             self._SamplingMCcycles, self._SamplingMC_iseed, self._SamplingMC_nworkers,
                             self._boxl, self._SamplingMC_delta, Tref=self._SamplingMC_Tref, logger=self._logger)
        else:
            # Setup Segments ==========================================================================
            self._setupSegments(input_dict['names'],
                                input_dict['filecoords'],
                                input_dict['filetop'], input_dict['filetypes'], calculate_volume=input_dict['calculate_volume'])

            # TASKS in parallel =======================================================================
            # Tasks for z-calculation and qm claculations are independent, so they can be performed in parallel.
            # The task_zcoord launches the calculation of all pairs in the local machine.
            # The task_qm sends or checks the calculations on the remote or local server.
            data = [[0, input_dict['coordination_numbers_Z']], [1, input_dict['interaction_energy']]]

            task_zcoord = mp.Process(name='daemon', target=self.__zcalc, args=(data[0],))
            task_zcoord.start()

            task_qm = mp.Process(name='non-daemon', target=self.__qmserver, args=(data[1],))
            task_qm.start()

            task_zcoord.join()
            task_qm.join()

        self._logger.info("\n\tFinnish: \t {}".format(datetime.datetime.now().strftime(fmt)))

    # ========================================================================
    def __zcalc(self, list_of_parameters):

        task = list_of_parameters[0] #0
        parm = list_of_parameters[1] #input_dict['coordination_numbers_Z']

        message = ("\n\t# Coordination number Z calculation (localhost)"
                   "\n\t-----------------------------------------------")
        if task == 0:
            # TASK 0: Calculation of the coordination numbers Z ==============
            # parm is input_dict['coordination_numbers_Z']
            if parm and not os.path.exists("Z_results.log") and not os.path.exists("Z_running_tmp.txt"):
                print(message) if self._logger is None else self._logger.info(message)
                Z_results = chi.Z_calc_group(self._names, self._filecoords, self._filetop, self._filetypes,
                                             self._logger, self._Z_samples,
                                             self._Z_putTrialsMonomer,
                                             self._Z_nonbondedMethod,
                                             self._Z_debug,
                                             workers=4)
                # Write Z results to a file
                chi.write_Z_results(Z_results)
            else:
                print(message) if self._logger is None else self._logger.info(message)
                message = "\tCoordination numbers won't be calculated.\n"
                if os.path.exists("Z_results.log"):
                    message += "\tZ values are already available in Z_results.log.\n"
                print(message) if self._logger is None else self._logger.info(message)

    # ========================================================================
    def __qmserver(self, list_of_parameters):

        task = list_of_parameters[0] #1
        parm = list_of_parameters[1] #input_dict['interaction_energy']

        message = ("\n\t# QM calculations at {}"
                   "\n\t-----------------------------------------------".format(self._servername))
        print(message) if self._logger is None else self._logger.info(message)
        if task == 1:

            if os.path.isfile(self._database_name):
                # If file Done exists do not either download nor update the database
                if not os.path.isfile("Done"):
                    # Check the qm calculations
                    self._check_state_interaction_energy_evaluation(self._database_name)

            # Calculation of the interaction energies
            elif parm:

                self._interaction_energy_evaluation(self._names, self._filecoords, self._filetop, self._filetypes)

    # ========================================================================
    def _parse(self, dict_input):

        R""" Parser for the input dictionary.

        The format of the input has been defined in :py:class:`Chi_Universe`

        Args:
            dict_input (dict): Input dictionary.

        Returns:
            True when the input is well-formed. Otherwise, return False.


        All list must have the same length
        """

        #==================
        message  = "\t# Parse input dictionary\n"
        message += "\t------------------------\n"

        if not 'names' in dict_input.keys():
            m = "\n\tField '%s' in input dictionary must exist." % 'names'
            m += "\tThe current dictionary is: \n\t{}".format(dict_input)
            print(m) if self._logger is None else self._logger.error(m)
            exit()
        elif len(dict_input['names']) != 2:
            m = "\nField '%s' in input dictionary must have exactly two items." % 'names'
            m += "Current value is: \n\t{}".format(dict_input['names'])
            print(m) if self._logger is None else self._logger.error(m)
            exit()
        else:
            self._names =dict_input['names']

        # ==================
        if not 'filecoords' in dict_input.keys():
            m = "\nField '%s' in input dictionary must exist." % 'filecoords'
            m += "The current dictionary is: \n\t{}".format(dict_input)
            print(m) if self._logger is None else self._logger.error(m)
            exit()
        elif len(dict_input['filecoords']) != 2:
            m = "\nField '%s' in input dictionary must have exactly two items." % 'filecoords'
            m += "Current value is: \n\t{}".format(dict_input['filecoords'])
            print(m) if self._logger is None else self._logger.error(m)
            exit()
        else:
            self._filecoords = dict_input['filecoords']

        # ==================
        if not 'filetop' in dict_input.keys():
            m = "\nField '%s' in input dictionary must exist." % 'filetop'
            m += "The current dictionary is: \n\t{}".format(dict_input)
            print(m) if self._logger is None else self._logger.error(m)
            exit()
        elif len(dict_input['filetop']) != 2:
            m = "\nField '%s' in input dictionary must have exactly two items." % 'filetop'
            m += "Current value is: \n\t{}".format(dict_input['filetop'])
            print(m) if self._logger is None else self._logger.error(m)
            exit()
        else:
            self._filetop = dict_input['filetop']

        # ==================
        if not 'filetypes' in dict_input.keys():
            self._filetypes = [None, None]
        else:
            self._filetypes = dict_input['filetypes']

        # ==================
        if not 'boxl' in dict_input.keys():
            self._boxl = 30.0 #in angstroms
        else:
            self._boxl = float(dict_input['boxl'])  # in angstroms

        # ==================
        if not 'sampling_run' in dict_input.keys():
            self._sampling_run = False
            message += "\Sampling run selection is FALSE.\n"
        else:
            if dict_input['sampling_run']:
                self._sampling_run = True
                message += "\tSampling run selection is TRUE.\n"
                try:
                    if dict_input['sampling_parameters']:
                        try:
                            self._Samplingtemperatures = dict_input['sampling_parameters']['temperatures']  # in K
                        except KeyError:
                            self._Samplingtemperatures = [300, 300, 1]  # in K
                        try:
                            self._Samplingmode = dict_input['sampling_parameters']['mode']
                            allowed_values = ["metropolis_sequential", "boltzmann_reweight", "metropolis_random"]
                            if not self._Samplingmode.lower() in allowed_values:
                                self._logger.error("Sampling mode '{}' unknown".format(self._Samplingmode))
                                self._logger.error("Sampling mode must be {}".format(allowed_values))
                                exit()
                        except KeyError:
                            self._logger.error("Sampling mode '{}' unknown".format(self._Samplingmode))
                            self._logger.error("Sampling mode must be {}".format(allowed_values))
                            exit()
                        try:
                            self._Samplingpairs = dict_input['sampling_parameters']['pairs']  # can be a list or a string
                        except KeyError:
                            self._Samplingpairs = 'all'
                        try:
                            self._SamplingMCcycles = dict_input['sampling_parameters']['MC_cycles']  # can be a list or a string
                        except:
                            self._SamplingMCcycles = 1
                        try:
                            self._SamplingMC_iseed = dict_input['sampling_parameters']['iseed_MC']  # can be a list or a string
                        except:
                            self._SamplingMC_iseed = None
                        try:
                            self._SamplingMC_nworkers = dict_input['sampling_parameters']['n_workers']  # can be a list or a string
                        except:
                            self._SamplingMC_nworkers = 1
                        try:
                            self._SamplingMC_delta = dict_input['sampling_parameters']['delta']  # width of the histograms in kcal/mol
                        except:
                            self._SamplingMC_delta = 0.2
                        try:
                            self._SamplingMC_Tref = dict_input['sampling_parameters']['T_reference']  # Reference temepratre for reweighting
                        except:
                            self._SamplingMC_Tref = None
                except KeyError:
                    self._logger.error("You must specify sampling_parameters input.")
                    exit()

            else:
                self._sampling_run = False
                message += "\tSampling run selection is FALSE.\n"

        # ==================
        if not 'coordination_numbers_Z' in dict_input.keys() or self._sampling_run:
            dict_input['coordination_numbers_Z'] = False
            message += "\tCoordination numbers is FALSE.\n"
        else:
            # Default or user parameters for the Z-coordination calculation
            message += "\tCoordination numbers is TRUE.\n"
            try:
                if dict_input['coordination_numbers_Z']:
                    try:
                        self._Z_samples = dict_input['Z_parameters']['Z_samples']
                    except KeyError:
                        self._Z_samples = 20

                    try:
                        self._Z_putTrialsMonomer = dict_input['Z_parameters']['Z_puttrialmonomers']
                    except KeyError:
                        self._Z_putTrialsMonomer = 1000

                    try:
                        self._Z_debug = dict_input['Z_parameters']['Z_debug']
                    except KeyError:
                        self._Z_debug = False

                    try:
                        self._Z_nonbondedMethod = dict_input['Z_parameters']['Z_nonbonded']
                    except KeyError:
                        self._Z_nonbondedMethod = 'truhlar'
            except KeyError:
                self._Z_samples = 20
                self._Z_putTrialsMonomer = 1000
                self._Z_debug = False
                self._Z_nonbondedMethod = 'truhlar'

        # ==================
        if not 'calculate_volume' in dict_input.keys() or self._sampling_run:
            message += "\tVolume of the segments is FALSE.\n"
            dict_input['calculate_volume'] = False
        else:
            if not dict_input['calculate_volume']:
                message += "\tVolume of the segments is FALSE.\n"
            else:
                message += "\tVolume of the segments is TRUE.\n"

        # ==================
        if not 'interaction_energy' in dict_input.keys() or self._sampling_run:
            message += "\tDFT Interaction energies between segments is FALSE.\n"
            dict_input['interaction_energy'] = False
            self._database_name = glob.glob("*.db")[0]
            self._qm_engine = dict_input['energy_parameters']['qm_engine'].lower()
        else:
            try:
                if dict_input['energy_parameters']:
                    dict_input['interaction_energy'] = True
                    try:
                        self._qm_engine = dict_input['energy_parameters']['qm_engine'].lower()
                        self._database_name = self._qm_engine+'_sp.db'
                    except KeyError:
                        self._qm_engine = "gaussian"
                        self._database_name = self._qm_engine+'_sp.db'
                    try:
                        self._qm_conformations = dict_input['energy_parameters']['number_configurations']
                    except KeyError:
                        self._qm_conformations = 4
                    try:
                        self._qm_method = dict_input['energy_parameters']['qm_method']
                    except KeyError:
                        self._qm_method = 'hf'
                    try:
                        self._qm_basisset = dict_input['energy_parameters']['qm_basisset']
                    except KeyError:
                        self._qm_basisset = '3-21g'
                    try:
                        self._qm_charge = dict_input['energy_parameters']['qm_charge']
                    except KeyError:
                        self._logger.error("You must specify the charge of the system ")
                        exit()
                    try:
                        if dict_input['energy_parameters']['qm_scratch_dir'].lower() == 'none':
                            self._qm_scratch_dir = "./"
                        else:
                            self._qm_scratch_dir = dict_input['energy_parameters']['qm_scratch_dir']
                            if self._qm_scratch_dir[-1] != "/":
                                self._qm_scratch_dir+="/"

                    except KeyError:
                        self._qm_scratch_dir = None
                    try:
                        self._qm_multiplicity = dict_input['energy_parameters']['qm_multiplicity']
                    except KeyError:
                        m ="\tYou must specify the multiplicity of the system "
                        print(m) if self._logger is None else self._logger.error(m)

                        exit()
                    try:
                        self._qm_task = dict_input['energy_parameters']['qm_task']
                    except KeyError:
                        self._qm_task = "energy"
                    try:
                        self._qm_otherkeys = dict_input['energy_parameters']['qm_otherkeys']
                    except KeyError:
                        self._qm_otherkeys = ""
                    try:
                        self._qm_path_exe = dict_input['energy_parameters']['qm_path_exe']
                    except KeyError:
                        m = "You must specify the path to {0:s} in " \
                                           "[energy_parameters input][qm_path_exe].".format(self._qm_engine)
                        print(m) if self._logger is None else self._logger.error(m)
                        exit()
                    try:
                        self._qm_memory_mb = dict_input['energy_parameters']['qm_memory_mb']
                    except KeyError:
                        self._qm_memory_mb = None #Mb

                else:
                    message += "\tDFT Interaction energies between segments is FALSE.\n"
                    dict_input['interaction_energy'] = False
                    self._database_name = glob.glob("*.db")[0]
                    self._qm_engine = dict_input['energy_parameters']['qm_engine'].lower()
            except KeyError:
                self._logger.error("You must specify energy_parameters input.")
                exit()

        # ==================
        if not 'server' in dict_input.keys():
            message += "\tServer: localhost.\n"
            dict_input['server'] = {'name': 'localhost', 'queue': None, 'username': None,
                                    'key_file': None}
            self._servername = dict_input['server']['name']
            self._username = dict_input['server']['username']
            self._key_file = dict_input['server']['key_file']
            self._local_dir = "./calculations_local_dir/"
            self._remote_dir = "./calculations_remote_dir/"
            self._nodelist = None
            self._nodelistmaster = None
            self._partition = None
            self._partitionmaster = None
            self._memory_per_task = None
            self._cpus_per_task = mp.cpu_count()
            # Check if the local directory exists, if False make the dir.
            if not os.path.isdir(self._local_dir):
                try:
                    os.mkdir(self._local_dir)
                except OSError:
                    m = "Creation of the directory %s failed" % self._local_dir
                    print(m) if self._logger is None else self._logger.error(m)
                    exit()
        else:
            try:
                self._servername = dict_input['server']['name']
                message += "\tServer in {}.\n".format(self._servername)
            except KeyError:
                m = "You must specify a name for the server in 'server:name' input."
                print(m) if self._logger is None else self._logger.error(m)
                exit()
            try:
                if self._servername.lower() == "localhost":
                    self._username = None
                else:
                    self._username = dict_input['server']['username']
            except KeyError:
                m = "You must specify a username for the server in 'server:username' input."
                print(m) if self._logger is None else self._logger.error(m)
                exit()
            try:
                self._key_file = dict_input['server']['key_file']
            except KeyError:
                self._key_file = None
            try:
                self._local_dir = dict_input['server']['local_dir']
                # local_dir must be a directory
                if self._local_dir[-1] != "/":
                    self._local_dir += "/"
                # Check if the directory exists, if False make the dir.
                if not os.path.isdir(self._local_dir):
                    try:
                        os.mkdir(self._local_dir)
                    except OSError:
                        m = "Creation of the directory %s failed" % self._local_dir
                        print(m) if self._logger is None else self._logger.error(m)
                        exit()
            except KeyError:
                self._local_dir = os.getcwd()
            try:
                self._remote_dir = dict_input['server']['remote_dir']
                # local_dir must be a directory
                if self._remote_dir[-1] != "/":
                    self._remote_dir += "/"
                if self._servername.lower() == "localhost":
                    if not os.path.exists(self._remote_dir):
                        m = "\tThe remote_dir: {0:s} does not exist in localhost\n".format(self._remote_dir)
                        m += "\tCreating remote dir: ./calculations_remote_dir/\n"
                        print(m) if self._logger is None else self._logger.warning(m)
                        self._remote_dir = "./calculations_remote_dir/"
                        try:
                            os.mkdir(self._remote_dir)
                        except OSError:
                            m = "Creation of the directory %s failed" % self._remote_dir
                            print(m) if self._logger is None else self._logger.error(m)
                            exit()
            except KeyError:
                m = "\tYou must specify a remote dir for the server in {0:}:{1:} input."\
                        .format(self._servername, self._username)
                print(m) if self._logger is None else self._logger.error(m)
                exit()
            try:
                self._cpus_per_task = dict_input['server']['ncpus']
            except KeyError:
                self._cpus_per_task = mp.cpu_count()
            try:
                self._memory_per_task = dict_input['server']['memory']
            except KeyError:
                self._memory_per_task = None
            try:
                self._nodelist = dict_input['server']['nodelist']
            except (KeyError, AttributeError):
                self._nodelist = None

            try:
                self._nodelistmaster = dict_input['server']['nodelistmaster']
            except (KeyError, AttributeError):
                self._nodelistmaster = None

            try:
                self._partition = dict_input['server']['partition']
            except (KeyError, AttributeError):
                self._partition = "simacro"

            try:
                self._partitionmaster = dict_input['server']['partitionmaster']
            except (KeyError, AttributeError):
                self._partitionmaster = self._partition

            try:
                self._max_queue_jobs = dict_input['server']['max_queue_jobs']
            except (KeyError, AttributeError):
                self._max_queue_jobs = 40

        print(message) if self._logger is None else self._logger.info(message)

    # ========================================================================
    def _setupSegments(self, listnames, listfilecoords, listfiletop, listfiletypes, calculate_volume=False):

        R""" Segments to perform the calculations

        Args:
            listnames (list of strings): List of names for the segments
            listfilecoords (list of strings): List of the paths for the coordinates files
            listfiletop (list of strings): List of the paths for the topology files
            calculate_volume (boolean, default=False): Calculate the volume
                using :py:meth:`chiripa.Segment.Segment.calc_vdw_volume_VABC()`.


        This method setup the **self._segments** property

        """

        message  = "\t# Summary of the segments\n"
        message += "\t-------------------------\n"

        # Create segments
        nsegments = len(listfilecoords)
        for i in range(0,nsegments):
            s = chi.Segment(self._names[i], filecoord=listfilecoords[i], filetop=listfiletop[i], filetypeatoms=listfiletypes[i])
            self._segments.append(s)
            if self._names is None:
                self._names.append(listnames[i])
            message += "\t==========\n"
            message += "\tSegment {}\n".format(i+1)
            message += "\t==========\n"
            message += "\tName        : {}\n".format(listnames[i])
            message += "\tCoordinate  : {}\n".format(listfilecoords[i])
            message += "\tTopology    : {}\n\n".format(listfiletop[i])

        if calculate_volume:
            message += "\t# Volume of the  segments\n"
            message += "\t-------------------------\n"
            for i in range(0,nsegments):
                s = self._segments[i]
                v_vdw, v_tsar = s.calc_vdw_volume_VABC()
                message += "\t VdW Volume of segment {:20s} = {:.1f} angstrom^3/molecule\n".\
                                  format(self._names[i], v_vdw)

        del self._segments
        print(message) if self._logger is None else self._logger.info(message)

    # ========================================================================
    def _interaction_energy_evaluation(self, listnames, listfilecoords, listfiletop, listfiletypes):

        # Write parameters  ===================================
        message  = "\t QM engine       = {} \n".format(self._qm_engine)
        message += "\t Number of pairs = {} \n".format(self._qm_conformations)
        message += "\n\t# Server info"\
                   "\n\t-------------\n"
        message += "\t Server Name     = {} \n".format(self._servername)
        message += "\t User Name       = {} \n".format(self._username)
        message += "\t Key file        = {} \n".format(self._key_file)
        message += "\t Local dir       = {} \n".format(self._local_dir)
        message += "\t Remote dir      = {} \n".format(self._remote_dir)
        print(message) if self._logger is None else self._logger.info(message)

        # QM Keywords
        qm_keywords = {'qm_engine': self._qm_engine.lower(),
                       'qm_method': self._qm_method,
                       'charge': self._qm_charge,
                       'scratch_dir': self._qm_scratch_dir,
                       'basis_set': self._qm_basisset,
                       'multiplicity': self._qm_multiplicity,
                       'qm_path': self._qm_path_exe,
                       'qm_task': self._qm_task.lower(),
                       'other_keys': self._qm_otherkeys,
                       'mem_mb': self._qm_memory_mb,
                       }

        if self._servername.lower() != 'localhost':

            # Connect to the server *********************************
            server = chi.ServerSlurm(self._servername, self._database_name,
                                     self._username, self._key_file)
            server.connection(self._logger)

        else:

            server = chi.ServerLocal(self._servername, self._database_name)


        # Inputs for monomers =================================================================
        for i in range(2):
            s2 = chi.Segment(self._names[i], filecoord=listfilecoords[i], filetop=listfiletop[i], filetypeatoms=listfiletypes[i])
            s1 = None

            c = chi.CubicBox(self._boxl, self._boxl, self._boxl)
            sim = chi.Simulator(mol=(s1, s2), box= c)

            zipfile = "{}_inputs_{:1d}_{:1d}.zip".format(qm_keywords['qm_engine'], 0, i+1)

            server.qm_calc(0, i+1, 1,
                           qm_keywords, sim,
                           self._local_dir, self._remote_dir,
                           self._cpus_per_task,
                           self._Z_nonbondedMethod,
                           nodelist=self._nodelist, nodelistmaster= self._nodelistmaster,
                           partition=self._partition, partitionmaster=self._partitionmaster,
                           localzipfile=zipfile, memory=self._memory_per_task,
                           logger=self._logger)
            del s1
            del s2
            del sim

        # Inputs for segments =================================================================
        for i in range(2):
            for j in range(2):
                s1 = chi.Segment(self._names[i], filecoord=listfilecoords[i], filetop=listfiletop[i], filetypeatoms=listfiletypes[i])
                s2 = chi.Segment(self._names[j], filecoord=listfilecoords[j], filetop=listfiletop[j], filetypeatoms=listfiletypes[j])
                c = chi.CubicBox(self._boxl, self._boxl, self._boxl)
                sim = chi.Simulator(mol=(s1, s2), box= c)

                zipfile = "{}_inputs_{:1d}_{:1d}.zip".format(qm_keywords['qm_engine'], i+1, j+1)

                server.qm_calc(i+1, j+1,self._qm_conformations,
                               qm_keywords, sim,
                               self._local_dir, self._remote_dir,
                               self._cpus_per_task,
                               self._Z_nonbondedMethod,
                               nodelist=self._nodelist, nodelistmaster=self._nodelistmaster,
                               partition=self._partition, partitionmaster=self._partitionmaster,
                               localzipfile=zipfile, memory=self._memory_per_task,
                               logger=self._logger)

                del s1
                del s2
                del sim

        if self._servername.lower() != 'localhost':
            # Send qm jobs to slurm
            server.qm_allsend(self._local_dir, self._remote_dir,
                              partitionmaster=self._partitionmaster,
                              nodelistmaster=self._nodelistmaster,
                              jobname="sendAll",
                              logger=self._logger, maxjobsslurm=self._max_queue_jobs)
        else:
            # Send qm jobs to slurm
            server.qm_allsend(qm_keywords['qm_engine'], self._remote_dir,
                              logger=self._logger)

        return None

    # ========================================================================
    def _check_state_interaction_energy_evaluation(self, databasename):

        fmt = "%a %d/%m/%Y, %H:%M:%S"
        self._logger.info("\n\t#QM# Checking the energy calculations in the remote server ( {} )".format(datetime.datetime.now().strftime(fmt)))

        message  = "\n\t# Server info"\
                   "\n\t-------------\n"
        message += "\t Server Name     = {} \n".format(self._servername)
        message += "\t User Name       = {} \n".format(self._username)
        message += "\t Key file        = {} \n".format(self._key_file)
        message += "\t Local dir       = {} \n".format(self._local_dir)
        message += "\t Remote dir      = {} \n".format(self._remote_dir)
        print(message) if self._logger is None else self._logger.info(message)

        start_time = datetime.datetime.now()
        if self._servername.lower() != 'localhost':

            server = chi.ServerSlurm(self._servername, databasename, self._username, self._key_file)
            server.connection(self._logger)

        else:

            server = chi.ServerLocal(self._servername, databasename=databasename)

        # 1. Extract energy from the log files.
        fmt = "%H:%M:%S"
        self._logger.info("\t Remote Server: Extract energies ( {} )".format( datetime.datetime.now().strftime(fmt)))
        # Monomers
        for i in range(2):
            zipdir = self._qm_engine+"_inputs_{:1d}_{:1d}/".format(0, i+1)
            server.extract_energy_calculations(self._local_dir,
                                                      self._remote_dir+zipdir, qm_engine=self._qm_engine)
        fmt = "%H:%M:%S"
        self._logger.info("\t\t Monomer energies extracted ( {} )".format(datetime.datetime.now().strftime(fmt)))
        # Pairs
        for i in range(2):
            for j in range(2):
                zipdir = self._qm_engine+"_inputs_{:1d}_{:1d}/".format(i+1, j+1)
                server.extract_energy_calculations(self._local_dir,
                                                          self._remote_dir+zipdir, qm_engine=self._qm_engine)
                fmt = "%H:%M:%S"
                self._logger.info("\t\t Pair {}-{} energies extracted ( {} )".format(i+1, j+1, datetime.datetime.now().strftime(fmt)))

        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        self._logger.info("\t Remote Server: Extract energies for all pairs in {:.1f} seconds\n"
                          .format(elapsed_time.total_seconds()))

        # 2. Get the files containing the energies.
        start_time = datetime.datetime.now()
        self._logger.info("\t Remote Server: Get summary files ( {} )".format(start_time.strftime(fmt)))
        server.get_energy_from_calculations(self._qm_engine, self._local_dir, self._remote_dir)
        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        self._logger.info("\t Remote Server: Get summary files for all pairs in {:.1f} seconds\n"
                          .format(elapsed_time.total_seconds()))

        # 3. Get the gaussian logs for future references.
        start_time = datetime.datetime.now()
        self._logger.info("\t Remote Server: Get log files ( {} )".format(start_time.strftime(fmt)))
        server.get_logs(self._local_dir, self._remote_dir)
        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        self._logger.info("\t Remote Server: Get log files for all pairs in {:.1f} seconds"
                          .format(elapsed_time.total_seconds()))

        if self._servername.lower() != 'localhost':
            server.close_connection()

        # 4. Check the database if all calculations have been perfomed and write results in the log
        total, inserver, pending, complete, error, running = server._db_base.account_qm("qm_jobs")

        # Monomers
        msg = "\n\t NAME         TOTAL_JOBS INSERVER_JOBS  PENDING_JOBS  RUNNING_JOBS  COMPLETE_JOBS  ERROR_JOBS\n"
        msg += "\t================================================================================================"
        self._logger.info(msg)
        k = 0
        msg = ""
        for i in range(2):
            msg += "\tMonomer {0:1d}  :\t{1:4d}\t\t{2:4d}\t\t{3:4d}\t\t\t{4:4d}\t\t\t{5:4d}\t\t{6:4d}\n".\
                format(k+1, total[k], inserver[k], pending[k], running[k], complete[k], error[k])
            k += 1

        # Pairs
        for i in range(2):
            for j in range(2):
                msg += "\tSegment {0:1d}_{1:1d}: \t{2:4d}\t\t{3:4d}\t\t{4:4d}\t\t\t{5:4d}\t\t\t{6:4d}\t\t{7:4d}\n".\
                    format(i+1, j+1, total[k], inserver[k], pending[k], running[k], complete[k], error[k])
                k += 1
        msg += "\t================================================================================================"
        self._logger.info(msg)

        is_done = True
        for k in range(len(total)):
            if total[k] != complete[k]+error[k]:
                is_done = False
        if is_done:
            with open("Done", 'w') as f:
                fmt = "%a %d/%m/%Y, %H:%M:%S"
                l = "QM work is done {} \n".format(start_time.strftime(fmt))
                f.writelines(l)

    # ========================================================================
    def __print_header(self):

        msg ="""                       
        ***********************************************************************
                         CHI inteRactIon PArameter (CHIRIPA)
                         -----------------------------------
                         
                                    Version 1.0
                         
                                  Dr. Javier Ramos
                          Macromolecular Physics Department
                    Instituto de Estructura de la Materia (IEM-CSIC)
                                   Madrid (Spain)
                                   
                Chiripa is an open-source python library to calculate
                the Flory-Huggins interaction parameter (chi) between
                two molecular segments. Chiripa is a coloquial spanish
                word meaning "stroke of luck"
                
                This software is distributed under the terms of the
                GNU General Public License v3.0 (GNU GPLv3). A copy of 
                the license (LICENSE.txt) is included with this distribution 
                
        ***********************************************************************
                     
        """

        self._logger.info(msg)

    # ========================================================================
    def _clean_files(self):

        if self._clean is not None:
            if self._clean.lower() == "full" or self._clean.lower() == "light":
                # Delete input files
                filelist1 = glob.glob("*_inputs_?_?.zip")
                # Delete summary files
                filelist2 = glob.glob("summary_*.txt")
                # Delete logs file
                filelist3 = glob.glob("logs.zip")
                # Delete sh scripts
                filelist4 = glob.glob("*.sh")
                # Delete sh scripts
                filelist5 = glob.glob("*.log")
                filelist = filelist1 + filelist2 + filelist3 +filelist4 +filelist5
                for ifile in filelist:
                    try:
                        os.remove(ifile)
                    except FileNotFoundError:
                        pass
                if self._clean.lower() == "full":

                    # Delete Z_results
                    try:
                        os.remove("./Z_results.log")
                    except FileNotFoundError:
                        pass

                    # Delete the database
                    ifile = glob.glob("*_sp.db")
                    try:
                        os.remove(ifile[0])
                    except FileNotFoundError:
                        pass
                    except IndexError:
                        pass

                return




