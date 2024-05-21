import sqlite3
import itertools
import pandas as pd
import numpy as np
from collections import defaultdict
import chiripa as chi
import datetime
import random
import time
import os
import multiprocessing as mp


class Sampling(object):

    """Samppling object"""

    # ========================================================================
    def __init__(self, dbname, sampling_mode, sampling_temperatures, mc_cycles, iseed,
                 nworkers, boxl, delta, Tref=None, logger=None):

        # Logger ========================================================================
        self._logger = logger
        # Connect the database ==========================================================
        self._dbname = dbname
        self._con = None
        self._cursor = None
        # Setup the self._con and self._cursor of the database ==========================
        self.__get_db_cursor(dbname)
        self._sampling_mode = sampling_mode
        # Create a datframe (self._dfl) =================================================
        self._dfl = None
        self.__create_dataframe()
        # Histogram using all QM calculations
        self._qmhist = defaultdict()
        # Create histograms of the interaction energies calculated for each pair ========
        self._energy_qmsampling(delta=delta)

        # Structures to store the results ===============================================
        self._chi12_avg = defaultdict()
        self._chi12_noovlp_avg = defaultdict()
        self._wij12_avg = defaultdict()
        self._chi_avg = defaultdict()
        self._chi_noovlp_avg = defaultdict()
        self._wij_avg = defaultdict()
        # 2 for pair 1-1, 3 for pair 1-2, 4 for pair 2-1, 5 for pair 2-2
        self._Ttargethist = {2: [], 3: [], 4: [], 5: []}
        self._energy_Tref = {2: [], 3: [], 4: [], 5: []}

        # Sampling ======================================================================
        if self._sampling_mode.upper() == "METROPOLIS_SEQUENTIAL":
            self._metropolis_results_sequential(sampling_temperatures, boxl, workers=nworkers, mc_cycles=mc_cycles, iseed=iseed)
            self.__write_results(logger)
        elif self._sampling_mode.upper() == "METROPOLIS_RANDOM":
            self._metropolis_results_random(sampling_temperatures, boxl, workers=nworkers, mc_cycles=mc_cycles, iseed=iseed)
            self.__write_results(logger)
        elif self._sampling_mode.upper() == "BOLTZMANN_REWEIGHT":
            if Tref == None:
                target_temp = 0
            else:
                target_temp = Tref
            self._boltzmann_reweight(sampling_temperatures, boxl, delta=delta, mc_cycles=mc_cycles,
                                     target_temp=target_temp, iseed=iseed)
            self.__write_results(logger)

        # Commit and close the database =================================================
        self.__update_db()

    # # ========================================================================
    def _energy_qmsampling(self, delta=0.2):

        energy_1 = self._dfl[0]['energy'].values[0]  # energy_monomer_1 (hartrees)
        energy_2 = self._dfl[1]['energy'].values[0]  # energy_monomer_2 (hartrees)

        # Energy monomer1, momomer2, monomer1+monomer1(1_1), monomer1+monomer2(1_2), monomer2+monomer1(2_1), monomer2+monomer2(2_2)
        energy_ref = [energy_1, energy_2, energy_1 * 2, energy_1 + energy_2, energy_2 + energy_1, energy_2 * 2]

        eij_avg = defaultdict()
        eij_min = defaultdict()
        eij_max = defaultdict()

        # NON-MONTE CARLO STATISTICS (ALL ITEMS ARE ACCEPTED) ===============================================
        samp_full = defaultdict()
        for ipair in range(2, 6):
            eij_avg[ipair] = 0.0
            samp_full[ipair] = 0
            for index, row in self._dfl[ipair].iterrows():
                try:
                    row['delta_energy'] = (row['energy'] - energy_ref[ipair]) * 627.5095
                    eij_avg[ipair] += row['delta_energy']
                    samp_full[ipair] += 1
                except TypeError:
                    # If energy has not been calculated
                    row['delta_energy'] = None

        for ipair in range(2, 6):

                eij_avg[ipair] /= samp_full[ipair]
                eij_max[ipair] = self._dfl[ipair]['delta_energy'].max()
                eij_min[ipair] = self._dfl[ipair]['delta_energy'].min()

                hmin = np.floor(eij_min[ipair])
                hmax = np.ceil(eij_max[ipair])
                nbins = int((hmax - hmin)/delta)
                self._qmhist[ipair] = np.histogram(self._dfl[ipair]['delta_energy'], bins=nbins, range=(hmin, hmax))

        self._write_histogram(self._qmhist, "eijqm_histogram_blocks.dat")

    # ========================================================================
    def _write_histogram(self, h, filename):

        # Check type of object
        if isinstance(h, dict):
            pass
        elif isinstance(h, tuple) and len(h)==2:
            # Conver tuple to dict
            tmp= {0: h}
            h = tmp
        else:
            return False

        # Histograms file
        with open(filename, "w") as f:
            l = "# All values in blocks of pairs. 2, 3, 4 and 5 corresponds to 1-1, 1-2, 2-1 and 2-2 pairs\n"
            f.writelines(l)

            for key, value in h.items():
                l = "# Histogram {}. ibin freq density\n".format(key)
                f.writelines(l)
                s = sum(value[0])
                for i in range(len(value[0])):
                    half = (value[1][i]+value[1][i+1])/2.0
                    l = "{0:.2f} {1:d} {2:.3f}\n".format(half, value[0][i], value[0][i]/s)
                    f.writelines(l)
                f.writelines("\n\n")

        return True

    # ========================================================================
    @staticmethod
    def _MC_worker(list_of_parameters):

        # Get data ==================================================
        dfl_local = list_of_parameters[0]
        T = list_of_parameters[1]
        seed = list_of_parameters[2]
        mc_cycles = list_of_parameters[3]
        boxl = list_of_parameters[4]
        logger = list_of_parameters[5]
        shuffle_list = list_of_parameters[6]

        # Seed for the random generator =============================
        if seed is None:
            iseed = int(982628732 * random.random())
        else:
            iseed = seed

        fmt = "%H:%M:%S"
        m = "\t\tStart T={} K. Number of MC runs = {} ( {} ) PID: {}. seed: {}".format(
                                      T, mc_cycles, datetime.datetime.now().strftime(fmt), os.getpid(), iseed)
        print(m) if logger is None else logger.info(m)


        energy_1 = dfl_local[0]['energy'].values[0]  # energy_monomer_1 (hartrees)
        energy_2 = dfl_local[1]['energy'].values[0]  # energy_monomer_2 (hartrees)

        # Energy monomer1, momomer2, monomer1+monomer1(1_1), monomer1+monomer2(1_2), monomer2+monomer1(2_1), monomer2+monomer2(2_2)
        energy_ref = [energy_1, energy_2, energy_1 * 2, energy_1 + energy_2, energy_2 + energy_1, energy_2 * 2]

        # MONTECARLO STATISTICS =========================================
        eij_mc = {2: [], 3: [], 4: [], 5: []}
        samp_mc = {2: [], 3: [], 4: [], 5: []}
        sij_mc = {2: [], 3: [], 4: [], 5: []}
        wij = []
        wij_12 = []
        max_sampling_energyvalues = {2: [], 3: [], 4: [], 5: []}
        samp_max = {2: 0, 3: 0, 4: 0, 5: 0}


        for icycle in range(0, mc_cycles):
            E_dictaccepted = {2: [], 3: [], 4: [], 5: []}

            for ipair in range(2, 6):
                E_current = 0.0; E_sum = 0.0; accepted = 0; E_listgenerated = [];  cog_list = []
                l_deltaenergy = dfl_local[ipair]['delta_energy'].tolist()
                l_cog = dfl_local[ipair]['cog'].tolist()
                if shuffle_list:
                    l = list(zip(l_deltaenergy,l_cog))
                    if seed is None:
                        iseed2 = int(982628732 * random.random())
                    else:
                        iseed2 = iseed
                    random.Random(iseed2).shuffle(l)
                    l_deltaenergy, l_cog = zip(*l)

                for i in range(len(l_deltaenergy)):
                    E_listgenerated.append(l_deltaenergy[i])
                    isaccepted, RT = chi.metropolis_MC(E_current, l_deltaenergy[i], T, iseed=iseed)
                    iseed += 1
                    if isaccepted:
                        E_current = l_deltaenergy[i]
                        E_dictaccepted[ipair].append(E_current)
                        cog_list.append([float(s) for s in  l_cog[i].split()])

                eij = np.mean(E_dictaccepted[ipair])
                eij_mc[ipair].append(eij)
                samp_mc[ipair].append(len(E_dictaccepted[ipair]))
                sij, sij_std = chi.testPointInsidePyramid(boxl, cog_list, iseed=iseed)
                sij_mc[ipair].append(sij)

                c_samp = len(E_dictaccepted[ipair])
                if c_samp > samp_max[ipair]:
                    max_sampling_energyvalues[ipair] = [E_dictaccepted[ipair], eij, sij]
                    samp_max[ipair] = c_samp

        # print(max_sampling_energyvalues)
        # for ipair in range(2, 6):
        #     print(icycle, len(max_sampling_energyvalues[ipair][0]))
        #     print(icycle, max_sampling_energyvalues[ipair][1], max_sampling_energyvalues[ipair][2])
        # print(samp_mc)
        fmt = "%H:%M:%S"
        m = "\t\t\t\tEnd T={} K. ( {} ) PID: {}".format(
                                      T, datetime.datetime.now().strftime(fmt), os.getpid())
        print(m) if logger is None else logger.info(m)

        return T, eij_mc, samp_mc, sij_mc, max_sampling_energyvalues

    # ========================================================================
    def _metropolis_results_sequential(self, sampling_temp, boxl, mc_cycles, workers=6, iseed=9837724721):

        # Processors info =====================================
        PROCESSES_AVAILABLE = int(mp.cpu_count())
        NUMBER_OF_TASKS = workers

        # Start time ==========================================
        start_time = datetime.datetime.now()
        fmt = "%a %d/%m/%Y, %H:%M:%S"

        # Temperatures ========================================
        Ti = sampling_temp[0]
        Tf = sampling_temp[1]+sampling_temp[2]
        deltaT = sampling_temp[2]
        nTemps = int((Tf-Ti)/deltaT)

        # Write parameters  ===================================
        message  = "\t# Sampling Algorithm parameters\n"
        message += "\t-------------------------------\n"
        message += "\t Starting time                                         = {} \n".format(start_time.strftime(fmt))
        message += "\t Number of MC Cyles                                    = {} \n".format(mc_cycles)
        message += "\t Number of temperatures from {}K to {}K in {}K steps   = {} \n".format(Ti, Tf, deltaT, nTemps)
        message += "\t Number of processes available                         = {} \n".format(PROCESSES_AVAILABLE)
        message += "\t Number of tasks in parallel                           = {} \n".format(NUMBER_OF_TASKS)
        message += "\t Metropolis mode                                       = sequential \n"

        print(message) if self._logger is None else self._logger.info(message)

        # Starting mp.pool ====================================
        pool_tasks = mp.Pool(processes=NUMBER_OF_TASKS)

        data = []
        for itemp in range(Ti, Tf, deltaT):
             data.append([self._dfl, itemp, iseed, mc_cycles, boxl, self._logger, False])

        results = pool_tasks.map(self._MC_worker, data)
        pool_tasks.close()

        self._collect_results_parallel(results, mc_cycles)

    # ========================================================================
    def _metropolis_results_random(self,sampling_temp, boxl, mc_cycles, workers=6, iseed=None):

        # Processors info =====================================
        PROCESSES_AVAILABLE = int(mp.cpu_count())
        NUMBER_OF_TASKS = workers

        # Start time ==========================================
        start_time = datetime.datetime.now()
        fmt = "%a %d/%m/%Y, %H:%M:%S"

        # Temperatures ========================================
        Ti = sampling_temp[0]
        Tf = sampling_temp[1]+sampling_temp[2]
        deltaT = sampling_temp[2]
        nTemps = int((Tf-Ti)/deltaT)

        # Write parameters  ===================================
        message  = "\t# Sampling Algorithm parameters\n"
        message += "\t-------------------------------\n"
        message += "\t Starting time                                         = {} \n".format(start_time.strftime(fmt))
        message += "\t Number of MC Cyles                                    = {} \n".format(mc_cycles)
        message += "\t Number of temperatures from {}K to {}K in {}K steps   = {} \n".format(Ti, Tf, deltaT, nTemps)
        message += "\t Number of processes available                         = {} \n".format(PROCESSES_AVAILABLE)
        message += "\t Number of tasks in parallel                           = {} \n".format(NUMBER_OF_TASKS)
        message += "\t Metropolis mode                                       = random\n"
        print(message) if self._logger is None else self._logger.info(message)

        # Starting mp.pool ====================================
        pool_tasks = mp.Pool(processes=NUMBER_OF_TASKS)

        data = []
        for itemp in range(Ti, Tf, deltaT):
             data.append([self._dfl, itemp, iseed, mc_cycles, boxl, self._logger, True])

        results = pool_tasks.map(self._MC_worker, data)
        pool_tasks.close()

        self._collect_results_parallel(results, mc_cycles)

    # ========================================================================
    def _collect_results_parallel(self, results, mc_cycles):

        ntemps = len(results)
        R = 1.987 * 1e-03  # kcal/molK

        # Get Z values from file Z_results.log
        Z_values = [0, 0]
        Z_std = [0, 0]
        with open('Z_results.log', 'r') as f:
            lines = f.readlines()
            for item in lines:
                Z_values.append(float(item.split("=")[1].split()[0]))
                Z_std.append(float(item.split("=")[1].split()[2]))

        # Create the file to save the results
        with open("eij_allvalues_blocks.dat", "w") as f:
            l = "# All values in blocks of temperatures\n"
            f.writelines(l)

        # df_results_dict
        df_results_dict = {"Pair":[], "T(K)":[],
                           "E_ij(kcal/mol)":[], "E_ij_std":[],
                           "S_ij":[], "S_ij_std": [],
                           "Z_ij":[], "Z_ij_std": [] }

        # Loop over the results
        for item in sorted(results, key=lambda x:x[0]):

            itemp = item[0]
            eij_mc = item[1]
            samp_mc = item[2]
            sij_mc = item[3]

            # Create a dataframe to store the results
            eij_avg = {2: [], 3: [], 4: [], 5: []}
            eij_std = {2: [], 3: [], 4: [], 5: []}
            sij_avg = {2: [], 3: [], 4: [], 5: []}
            sij_std = {2: [], 3: [], 4: [], 5: []}

            RT = itemp * R  # kcal/mol

            pair_list = []; itemp_list = []
            eij_avg_list = []; eij_std_list = []
            sij_avg_list = []; sij_std_list = []
            zij_list = []; zij_std_list = []

            for ipair in range(2, 6):
                eij_avg[ipair] = np.mean(eij_mc[ipair])
                eij_std[ipair] = np.std(eij_mc[ipair])
                sij_avg[ipair] = np.mean(sij_mc[ipair])
                sij_std[ipair] = np.std(sij_mc[ipair])
                i, j = self._dfl[ipair].iloc[0, 1].split("_")[-2:]

                df_results_dict["Pair"] += ["{}-{}".format(i, j)]
                df_results_dict["T(K)"] += [itemp]
                df_results_dict["E_ij(kcal/mol)"] += [eij_avg[ipair]]
                df_results_dict["E_ij_std"] += [eij_std[ipair]]
                df_results_dict["S_ij"] += [sij_avg[ipair]]
                df_results_dict["S_ij_std"] += [sij_std[ipair]]
                df_results_dict["Z_ij"] += [Z_values[ipair]]
                df_results_dict["Z_ij_std"] += [Z_std[ipair]]

            self._chi12_avg[itemp] = 2.0 * eij_avg[3] * Z_values[3] * sij_avg[3] - \
                                     (eij_avg[2] * Z_values[2] * sij_avg[2] + eij_avg[5] * Z_values[5] * sij_avg[5])
            if itemp != 0:
                self._chi12_avg[itemp] = 0.5 * self._chi12_avg[itemp] / RT
            else:
                self._chi12_avg[itemp] = 0.0
            self._wij12_avg[itemp] = 2.0 * eij_avg[3] - (eij_avg[2] + eij_avg[5])
            self._chi12_noovlp_avg[itemp] = 2.0 * eij_avg[3] * Z_values[3] - \
                                            (eij_avg[2] * Z_values[2] + eij_avg[5] * Z_values[5])
            if itemp != 0:
                self._chi12_noovlp_avg[itemp] = 0.5 * self._chi12_noovlp_avg[itemp] / RT
            else:
                self._chi12_noovlp_avg[itemp] = 0.0

            self._chi_avg[itemp] = (eij_avg[3] * Z_values[3] * sij_avg[3] + eij_avg[4] * Z_values[4] * sij_avg[4]) - \
                                   (eij_avg[2] * Z_values[2] * sij_avg[2] + eij_avg[5] * Z_values[5] * sij_avg[5])
            if itemp != 0:
                self._chi_avg[itemp] = 0.5 * self._chi_avg[itemp] / RT
            else:
                self._chi_avg[itemp] = 0.0
            self._wij_avg[itemp] = (eij_avg[3] + eij_avg[4]) - (eij_avg[2] + eij_avg[5])
            self._chi_noovlp_avg[itemp] = (eij_avg[3] * Z_values[3] + eij_avg[4] * Z_values[4]) - \
                                          (eij_avg[2] * Z_values[2] + eij_avg[5] * Z_values[5])
            if itemp != 0:
                self._chi_noovlp_avg[itemp] = 0.5 * self._chi_noovlp_avg[itemp] / RT
            else:
                self._chi_noovlp_avg[itemp] = 0.0

            # Write files
            if self._sampling_mode.upper() == "METROPOLIS_SEQUENTIAL" or \
               self._sampling_mode.upper() == "METROPOLIS_RANDOM":
                self.__write_energy_blocks(itemp, mc_cycles, eij_mc, samp_mc, sij_mc, Z_values)

        self._df_results = pd.DataFrame.from_dict(df_results_dict)

    # ========================================================================
    def __write_energy_blocks(self, itemp, nseeds_run, eij_mc, samp_mc, sij_mc, Z_values):

        R = 1.987 * 1e-03
        RT = itemp * R

        # Write all values
        with open("eij_allvalues_blocks.dat", "a") as f:
            l = "#iseed\te_11\tsamp_11\te_12\tsamp_12\te_21\tsamp_21\te_22\tsamp_22\ts11\ts12\ts21\ts22" \
                "\twij\twij_12\tFlory12_noovlp\tFlory12_ovlp\tFlory_noovlp\tFlory_ovlp (at {}K)\n".format(
                itemp)
            f.writelines(l)
            for iseed in range(0, nseeds_run):

                eij_Zij_3 = eij_mc[3][iseed] * Z_values[3]
                eij_Zij_4 = eij_mc[4][iseed] * Z_values[4]
                eij_Zij_2 = eij_mc[2][iseed] * Z_values[2]
                eij_Zij_5 = eij_mc[5][iseed] * Z_values[5]
                chi_12_novlp = 2.0 * eij_Zij_3 -  ( eij_Zij_2 + eij_Zij_5 )
                chi_12_ovl= 2.0 * eij_Zij_3 * sij_mc[3][iseed] -  ( eij_Zij_2 * sij_mc[2][iseed] + eij_Zij_5 * sij_mc[5][iseed] )
                chi_all_novlp = (eij_Zij_3 + eij_Zij_4) - (eij_Zij_2 + eij_Zij_5)
                chi_all_ovl = (eij_Zij_3 * sij_mc[3][iseed] + eij_Zij_4 * sij_mc[4][iseed]) - (eij_Zij_2 * sij_mc[2][iseed] + eij_Zij_5 * sij_mc[5][iseed])

                l = "{0:5d}\t{1:.3f}\t{2:5d}\t{3:.3f}\t{4:5d}\t{5:.3f}\t{6:5d}\t{7:.3f}\t{8:5d}\t{9:.3f}\t{10:.3f}" \
                    "\t{11:.3f}\t{12:.3f}\t{13:.3f}\t{14:.3f}\t{15:.3f}\t{16:.3f}\t{17:.3f}\t{18:.3f}\n". \
                    format(iseed, eij_mc[2][iseed], samp_mc[2][iseed],
                           eij_mc[3][iseed], samp_mc[3][iseed],
                           eij_mc[4][iseed], samp_mc[4][iseed],
                           eij_mc[5][iseed], samp_mc[5][iseed],
                           sij_mc[2][iseed], sij_mc[3][iseed],
                           sij_mc[4][iseed], sij_mc[5][iseed],
                           eij_mc[3][iseed] + eij_mc[4][iseed] - eij_mc[2][iseed] - eij_mc[5][iseed],
                           eij_mc[3][iseed] - 0.5 * (eij_mc[2][iseed] + eij_mc[5][iseed] ),
                           0.5*chi_12_novlp/RT,
                           0.5*chi_12_ovl/RT,
                           0.5*chi_all_novlp/RT,
                           0.5*chi_all_ovl/RT)

                f.writelines(l)
            f.writelines("\n\n")

    # ========================================================================
    def __write_results(self, logger):

        m = "\n\t#Sampling results\n"
        m += "\t#===============\n"
        m += '\t{}\n'.format(self._df_results.to_string(index = False, formatters={"#Pair": "\t{:}".format,
                                                                              "T(K)": "\t{:6,.1f}".format,
                                                                              "E_ij(kcal/mol)": "\t{:6,.4f}".format,
                                                                              "E_ij_std": "\t{:6,.4f}".format,
                                                                              "S_ij": "\t{:6,.4f}".format,
                                                                              "S_ij_std": "\t{:6,.4f}".format,
                                                                              "Z_ij": "\t{:6,.2f}".format,
                                                                              "Z_ij_std": "\t{:6,.2f}".format,}))
        m += "\t#===================\n"
        logger.info(m)

        m = "\t#Flory-Hugging parameters all pairs\n"
        m += "\t#==================================\n"
        m += "\t#Temp(K)\t\tchi value_noovlp\t\tchi value_ovlp\t\tw_ij\n"

        for itemp, ichi in self._chi_avg.items():
            m+="\t{0:>5.1f}\t\t{1:>8.3f}\t\t{2:>8.3f}\t\t{3:>8.3f}\n".\
                format(float(itemp), self._chi_noovlp_avg[itemp], self._chi_avg[itemp], self._wij_avg[itemp])
        m += "\t#========================\n"
        logger.info(m)

        m = "\t#Flory-Hugging parameters 11, 12, 22 pairs\n"
        m += "\t#==================================\n"
        m += "\t#Temp(K)\t\tchi value_noovlp_12\t\tchi value_ovlp_12\t\tw_ij_12\n"
        for itemp, ichi in self._chi12_avg.items():
            m+="\t{0:>5.1f}\t\t{1:>8.3f}\t\t{2:>8.3f}\t\t{3:>8.3f}\n".\
                format(float(itemp), self._chi12_noovlp_avg[itemp], self._chi12_avg[itemp], self._wij12_avg[itemp])
        m += "\t#========================\n"
        logger.info(m)

    # ========================================================================
    def _boltzmann_reweight(self, sampling_temp, boxl, delta=0.2, target_temp=0, mc_cycles=0, iseed=None):

        # Processors info =====================================
        PROCESSES_AVAILABLE = int(mp.cpu_count())
        NUMBER_OF_TASKS = 1

        # Start time ==========================================
        start_time = datetime.datetime.now()
        fmt = "%a %d/%m/%Y, %H:%M:%S"

        # Temperatures ========================================
        Ti = sampling_temp[0]
        Tf = sampling_temp[1]+sampling_temp[2]
        deltaT = sampling_temp[2]
        nTemps = int((Tf-Ti)/deltaT)

        # Write parameters  ===================================
        message  = "\t# Sampling Algorithm parameters\n"
        message += "\t-------------------------------\n"
        message += "\t Starting time                                          = {} \n".format(start_time.strftime(fmt))
        message += "\t Number of MC Cyles                                     = {} \n".format(mc_cycles)
        message += "\t Number of temperatures from {0:3d}K to {1:3d}K in {2:3d}K steps = {3:<3d} \n".format(Ti, Tf, deltaT, nTemps)
        message += "\t Number of processes available                          = {} \n".format(PROCESSES_AVAILABLE)
        message += "\t Number of tasks in parallel                            = {} \n".format(NUMBER_OF_TASKS)
        if target_temp != 0:
            message += "\t Boltzmann Reweight Target Temperature                  = {0:3d}K\n".format(target_temp)
        else:
            message += "\t Boltzmann Reweight (from QM histograms)                = {0:3d}K\n".format(target_temp)

        print(message) if self._logger is None else self._logger.info(message)

        if target_temp != 0:

            itemp = target_temp
            # Starting mp.pool ====================================
            pool_tasks = mp.Pool(processes=NUMBER_OF_TASKS)
            data = []
            mc_cycles = 5
            data.append([self._dfl, itemp, iseed, mc_cycles, boxl, self._logger, True])

            results = pool_tasks.map(self._MC_worker, data)
            pool_tasks.close()

            # For each pair
            sij_list = []
            for ipair in range(2, 6):

                self._energy_Tref[ipair] = results[0][4][ipair][0]
                max_v = np.max(self._energy_Tref[ipair])
                min_v = np.min(self._energy_Tref[ipair])
                hmin = np.floor(min_v)
                hmax = np.ceil(max_v)
                nbins = int((hmax - hmin) / delta)
                self._Ttargethist[ipair] = np.histogram(self._energy_Tref[ipair], bins=nbins, range=(hmin, hmax))
                sij_list.append(np.mean(results[0][3][ipair]))

        else:
            sij_list = []
            for ipair in range(2,6):
                cog_list = []
                l_cog = self._dfl[ipair]['cog'].tolist()
                for i in range(len(l_cog)):
                    cog_list.append([float(s) for s in l_cog[i].split()])

                sij, sij_std = chi.testPointInsidePyramid(boxl, cog_list, iseed=iseed)
                sij_list.append(sij)
            self._Ttargethist = self._qmhist

        # Write histograms and calculate averages for other temperatures
        eij_avg_temp = self.__write_reweight_values(target_temp, sampling_temp, delta)

        # # Create a dataframe to store the results
        if target_temp == 0:
            results = []
            # Temperature 0K
            a = (target_temp, {2: [eij_avg_temp[target_temp][0]], 3: [eij_avg_temp[target_temp][1]],
                         4: [eij_avg_temp[target_temp][2]], 5: [eij_avg_temp[target_temp][3]]},
                 {2: [], 3: [], 4: [], 5: []},
                 {2: [sij_list[0]], 3: [sij_list[1]], 4: [sij_list[2]], 5: [sij_list[3]]},
                 {2: [], 3: [], 4: [], 5: []})
            results.append(a)

        for itemp in range(Ti, Tf, deltaT):

            a = (itemp, {2:[eij_avg_temp[itemp][0]], 3:[eij_avg_temp[itemp][1]],
                         4:[eij_avg_temp[itemp][2]], 5:[eij_avg_temp[itemp][3]]},
                 {2: [], 3: [], 4: [], 5: []},
                 {2: [sij_list[0]], 3: [sij_list[1]], 4: [sij_list[2]], 5: [sij_list[3]]},
                 {2: [], 3: [], 4: [], 5: []})
            results.append(a)


        self._collect_results_parallel(results, mc_cycles)

    # ========================================================================
    def __write_reweight_values(self, Tref, sampling_temp, delta):

        # Temperatures ========================================
        Ti = sampling_temp[0]
        Tf = sampling_temp[1]+sampling_temp[2]
        deltaT = sampling_temp[2]
        nTemps = int((Tf-Ti)/deltaT)

        # R value ========================================
        R = 1.987 * 1e-03  # kcal/molK
        if Tref != 0:
            beta_ref = 1.0/(R*Tref) #kcal/mol
        else:
            beta_ref = 0

        # Accumulate averages for target temperature
        #
        #   <f(E)> = Sum(i=1..N) f(Ei)*Pbeta(Ei) / Sum(i=1..N) Pbeta(Ei)
        #  where Pbeta(E) is the frequency value for energy Ei
        eij_avg = defaultdict(list)
        for ipair in range(2, 6):
            #  Pbeta(Ei) * Ei
            n = np.sum(self._Ttargethist[ipair][0][0:] * self._Ttargethist[ipair][1][0:-1])
            # Pbeta(Ei)
            m = np.sum(self._Ttargethist[ipair][0][0:])
            eij_avg[Tref].append(n / m)

        # Write reference histogram
        with open("rw_histograms_Tref_{0:03d}K.dat".format(Tref), 'w') as f:
            f.writelines("# Histograms for reweighting boltzmann\n")
            f.writelines("# e11\tfreqe11\tHe11\te12\tfreqe12\tHe12\te21\tfreqe21\tHe21\te22\tfreqe22\tHe22\t "
                         "(at reference T: {}K, delta: {}kcal/mol)\n".format(Tref, delta))
            # Get maximum number of items in the histograms
            max_items = 0
            for ipair in range(2, 6):
                l = len(self._Ttargethist[ipair][0])
                if l>max_items:
                    max_items = l
            # Write histograms
            lines = []
            for item in range(0, max_items):
                l = ""
                for ipair in range(2, 6):
                    try:
                        n = np.sum(self._Ttargethist[ipair][0])
                        l += "{0:>6.1f}  {1:>6d}  {2:>5.3f}  ".format(self._Ttargethist[ipair][1][item],
                                                                     self._Ttargethist[ipair][0][item],
                                                                     self._Ttargethist[ipair][0][item]/n)
                    except IndexError:
                        l += "{0:^6s}  {1:>6s}  {2:^4s}   ".format("nan", "nan", "nan")
                l += "\n"
                f.writelines(l)

        # Calculate and write reweigthed energy values for histograms
        p_beta = {2: [], 3: [], 4: [], 5: []}
        for itemp in range(Ti, Tf, deltaT):
            if itemp == Tref: continue
            beta = 1.0/(itemp * R)  # kcal/mol
            delta_beta = beta - beta_ref #(kcal/mol^-1)
            for ipair in range(2,6):
                p_beta[ipair] = self._Ttargethist[ipair][0][0:]*np.exp(-delta_beta*self._Ttargethist[ipair][1][0:-1])
                sum_p_beta = np.sum(p_beta[ipair])
                sum_e_p_beta = np.sum(self._Ttargethist[ipair][1][0:-1]*p_beta[ipair])
                eij_avg[itemp].append(sum_e_p_beta/sum_p_beta)


            # Write reference histogram
            with open("rw_histograms_Tref_{0:03d}K.dat".format(Tref), 'a') as f:
                f.writelines("\n\n")
                f.writelines("# e11\tfreqe11\tHe11\te12\tfreqe12\tHe12\te21\tfreqe21\tHe21\te22\tfreqe22\tHe22\t "
                             "(at reference T: {}K, delta: {}kcal/mol)\n".format(itemp, delta))
                # Get maximum number of items in the histograms
                max_items = 0
                for ipair in range(2, 6):
                    l = len(p_beta[ipair])
                    if l>max_items:
                        max_items = l
                # Write histograms
                lines = []
                for item in range(0, max_items):
                    l = ""
                    for ipair in range(2, 6):
                        try:
                            n = np.sum(p_beta[ipair])
                            l += "{0:>6.1f}  {1:>6.1f}  {2:>5.3f}  ".format(self._Ttargethist[ipair][1][item],
                                                                         p_beta[ipair][item],
                                                                         p_beta[ipair][item]/n)

                        except IndexError:
                            l += "{0:^6s}  {1:>6s}  {2:^4s}   ".format("nan", "nan", "nan")
                    l += "\n"
                    f.writelines(l)



        print("cJJ:", eij_avg)

        return eij_avg

    # ========================================================================
    def __get_db_cursor(self, dbname):

        self._con = sqlite3.connect(dbname)
        self._cursor = self._con.cursor()

        # Add column delta energy if does not exist ===================================
        try:
            self._cursor.execute("ALTER TABLE qm_jobs ADD delta_energy FLOAT")
        except sqlite3.OperationalError:
            pass
        try:
            self._cursor.execute("ALTER TABLE qm_jobs ADD Sampling_accepted VARCHAR")
        except sqlite3.OperationalError:
            pass

        # Update delta_energy values in the database
        # Create a dataframe ==========================================================
        sql = "SELECT *  FROM qm_jobs"
        df = pd.read_sql_query(sql, self._con)

        # Divide the dataframe in six pieces for monomer1, monomer2, 1_1, 1_2, 2_1 and 2_2 pairs =========
        dfl = list()
        dfl.append(df[df['name_job'].str.contains("_0_1")])
        dfl.append(df[df['name_job'].str.contains("_0_2")])
        dfl.append(df[df['name_job'].str.contains("_1_1")])
        dfl.append(df[df['name_job'].str.contains("_1_2")])
        dfl.append(df[df['name_job'].str.contains("_2_1")])
        dfl.append(df[df['name_job'].str.contains("_2_2")])

        energy_1 = dfl[0]['energy'].values[0]  # energy_monomer_1 (hartrees)
        energy_2 = dfl[1]['energy'].values[0]  # energy_monomer_2 (hartrees)

        # Energy monomer1, momomer2, monomer1+monomer1(1_1), monomer1+monomer2(1_2), monomer2+monomer1(2_1), monomer2+monomer2(2_2)
        energy_ref = [energy_1, energy_2, energy_1 * 2, energy_1 + energy_2, energy_2 + energy_1, energy_2 * 2]

        # Delta-energy row in database
        for ipair in range(2, 6):
            for index, row in dfl[ipair].iterrows():
                try:
                    row['delta_energy'] = (row['energy'] - energy_ref[ipair]) * 627.5095
                except TypeError:
                    # If energy has not been calculated
                    row['delta_energy'] = 1000.0
                namejob = row['name_job']
                sql = """UPDATE qm_jobs SET delta_energy = '{}' WHERE name_job = '{}'""".format(row['delta_energy'], namejob)
                self._cursor.execute(sql)

        self._con.commit()
        del(dfl)
        del(df)

    # # ========================================================================
    def __update_db(self):

        self._con.commit()
        self._con.close()

    # ========================================================================
    def __create_dataframe(self):

        # Create a dataframe ==========================================================
        sql = "SELECT *  FROM qm_jobs"
        df = pd.read_sql_query(sql, self._con)

        # Divide the dataframe in six pieces for monomer1, monomer2, 1_1, 1_2, 2_1 and 2_2 pairs =========
        self._dfl = list()
        self._dfl.append(df[df['name_job'].str.contains("_0_1")])
        self._dfl.append(df[df['name_job'].str.contains("_0_2")])
        self._dfl.append(df[df['name_job'].str.contains("_1_1")])
        self._dfl.append(df[df['name_job'].str.contains("_1_2")])
        self._dfl.append(df[df['name_job'].str.contains("_2_1")])
        self._dfl.append(df[df['name_job'].str.contains("_2_2")])

