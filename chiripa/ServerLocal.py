import paramiko
import os
import datetime
import subprocess
import shutil
from socket import gaierror
import chiripa as chi

"""
Some functions are inspired or directly copy from slurmqueen project

(https://github.com/Kasekopf/SlurmQueen)

"""

class ServerLocal(chi.Server):

    # ===========================================================================================
    def __init__(self, nameserver, databasename, logger=None):

        """
        Initialize a local server without queue system

        Args:
            nameserver (str): The name of the server (i.e: trueno.csic.es or localhost).
            databasename (str) :  File name of the database.

        """

        super().__init__(nameserver, username=None, keyfile=None,
                         dbname=databasename, logger=None)

        if not self._db_base is None:
            self._db_base.create_table_qmjobs()

    # ===========================================================================================
    def execute_cmd(self, command, other_input=None, timeout=10):

        """
        Execute a command (from the WORKING directory) on this server.
        Throws an exception if the connection fails.

        Args:
            command (str): The name to execute.
            other_input (str): Text to be written to stdin.
            timeout (int): A timeout to use when executing the command, in seconds.

        Returns:
            tuple: All text written to stdout by the command
            as a tuple of strings (stdout, stderr).

        Examples:
            A typical use of the method.

            >>> remote_dir = "./"
            >>> zipfilename = "file.zip"
            >>> cmd2 = 'cd %s; rm -f %s'%(remote_dir,zipfilename)
            >>> stdout, stderr = self.execute_cmd(cmd2)

        """

        p = subprocess.Popen(command, shell=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        exit_status = 0
        if len(err) != 0:
                exit_status = 2

        if exit_status != 0:
            print("\tWarning!!! Command {} cannot be run in server localhost".format(command))

        return out.decode("utf8"), err.decode("utf8")

    # ===========================================================================================
    def qm_calc(self, iseg, jseg, ncalc,
                qm_keywords, simulator,
                localdir, remotedir,
                cpuspertask,
                nonbondedEvaluation,
                nodelist=None, nodelistmaster=None,
                partition=None, partitionmaster=None,
                localzipfile=None, memory=None,
                logger=None):
        """
        Prepare the data to send Gaussian16 jobs on a localserver

        #. Create the database if it does not exist
        #. Define the number of nodes used in the calculations
        #. For each pair segement replica:

            a. Write gaussian input
            b. Write script to send the calculation
            c. Insert data in the gaussian database

        #.  Send gaussian jobs to localhost

        Args:
            iseg (int): Index for the Base Segment
            jseg (int): Index for the Screening Segment
            ncalc (int): Number of replicas for the pair iseg-jseg.
            gaussian_keywords (str): Route section (# lines) for `gaussian 16 <https://gaussian.com/input/>`_
                (Specify desired calculation type, model chemistry, and other options)
            simulator (Simulator): A :doc:`Simulator <./simulator>` instance.
            localdir (str): Path to store files in the local server.
            remotedir (str): Path to store files in the remote server.
            g16path (str): Path the to the gaussian executable.
            scratch_dir (str): Path to the scratch dir used for gaussian
            cpuspertask (int): Number of cpus to use in the gaussian calculation
            nodelist (list): List of nodes. For localhost should be ``None`` **(default=None)**.
            localzipfile (str): Path for a zip file of the inputs **(default="gaussian_inputs.zip")**
            memory (str): Memory (example: 500Mb)
            logger (Logger): A :doc:`Logger <./setup_logger>` instance

        Examples:
            Use example

            >>> gaussian_keywords = '#p 3-21g hf int=grid=ultrafine scf=(maxcycle=5000)'
            >>> self._local_dir = '/tutorials/tutorial01 gauss local/calculations_local_dir/'
            >>> self._remote_dir = './calculations_remote_dir/'
            >>> self._qm_path_exe = '/opt/g16/g16'
            >>> memory = "500Mb"
            >>> logger = chi.init_logger("Output", fileoutput = "output_check.log", append=True, inscreen=output_screen)
            >>> server.gaussian_calc(0, 1, self._qm_conformations,
                                     gaussian_keywords,
                                     sim, self._local_dir, self._remote_dir,
                                     self._qm_path_exe,
                                     self._qm_scratch_dir,
                                     self._cpus_per_task,
                                     nodelist=None,
                                     localzipfile=zipfile,
                                     memory=str(self._memory_per_task) + "Mb", logger=None)


        """

        file_base = localzipfile.split("/")[-1].split(".")[0]
        cog_list = list()

        # Create the database
        if self._db_base is None:
            self._db_base = DBjobs("{}_sp.db".format(qm_keywords['qm_engine']), logger=logger)
            self._db_base.create_table_qmjobs()

        #print("1:", datetime.datetime.now().time())

        # Create gaussian and script files in the local area
        inode = 0
        try:
            nnodes = len(nodelist)
        except (AttributeError, TypeError):
            nnodes = 1

        for iconf in range(ncalc):

            if ncalc > 1:
                chi.gen_conf_fan_algorithm_c(simulator, nonbondedEvaluation)

            cog_list = simulator.get_ith_molecule(-1).center_of_geom()
            inputname = "{0:s}_{1:05d}_{2:1d}_{3:1d}".format(qm_keywords['qm_engine'], iconf, iseg, jseg)
            idir = remotedir+"/"+file_base+"/{}".format(inputname)
            simulator.write_qminput_simulator(idir, qm_keywords,
                                              inputname,
                                              delete_ok=True,
                                              mem=memory,
                                              nproc=cpuspertask,
                                              is_send_script="localhost")

            # Insert data in the gaussian database
            sql_insert = "INSERT INTO qm_jobs VALUES ({0:d}, '{1:s}', 'CREATED', 'NULL', 'NULL', '{2:s}', 'NULL')".\
                format(self._db_index, inputname, str(cog_list)[1:-1])
            self._db_base.insert_data(sql_insert)
            self._db_index += 1
            inode = (inode+1) % nnodes

        #print("6:",datetime.datetime.now().time())

        d = {'status_job' : 'INSERVER'}
        self._db_base.update_allcolumns_tosamevalue("qm_jobs", d)
        self._db_base.commit_db()

        #print("8:",datetime.datetime.now().time())
        # Send gaussian jobs to localhost
        idir_remote = remotedir+file_base+"/"

    # ===========================================================================================
    def qm_allsend(self, qm_engine, remote_dir, logger=None):

        """
        Send all gaussian jobs. This only has been tested in Linux

        Args:
            remote_dir (str): Path to store files in the remote server.

        """

        self.generate_bashscript_send_local(qm_engine,  remote_dir, logger=logger)

        cmd = "cd {}; bash full_send.sh".format(remote_dir)
        os.system(cmd)

    # ===========================================================================================
    def get_logs(self, local_dir, remote_dir):

        """
        This method is used to get all gaussian logs in a
        zip file called "logs.zip"

        Args:
            local_dir (str): Path to store files in the local server.
            remote_dir (str): Path to store files in the remote server.

        """

        zipfilename = "logs.zip"
        #os.mkdir(local_dir+"gaussian_logs")

        cmd1 = 'cd %s; zip -R %s \'*.log\''%(remote_dir,zipfilename)
        stdout, stderr = self.execute_cmd(cmd1)
        #
        cmd3 = 'cd %s; mv %s %s'%(remote_dir, zipfilename, local_dir)
        stdout, stderr = self.execute_cmd(cmd3)

    # ===========================================================================================
    def get_energy_from_calculations(self, qm_engine, local_dir, remote_dir):

        """
        Get the summary files and update the database with the energies from log files.
        This is for localhosts

        Args:
            qm_engine (str): Name of the QM package (example: gaussian)
            local_dir (str): Path to store files in the local server.
            remote_dir (str): Path to store files in the remote server.

        """

        energyfilename = "summary_energy.txt"
        summaryfilename = "summary.txt"

        zipfilename = "summary.zip"
        cmd1 = 'cd %s; ls -R *_inputs*/summary*'%remote_dir

        stdout, stderr = self.execute_cmd(cmd1)

        # Update the database
        # Monomers
        for i in range(2):
            subdir = qm_engine+"_inputs_{:1d}_{:1d}/".format(0, i+1)
            name2 = "summary_energy.txt".format(0, i+1)
            name3 = "summary.txt".format(0, i+1)
            fulldir = remote_dir+subdir
            if os.path.exists(fulldir+name2):
                with open(fulldir+name2, 'r') as f:
                    lines = f.readlines()
                    cmd = 'cp %s %s_0_%d.txt'%(fulldir+name2, local_dir+name2.split(".")[0],i+1)
                    stdout, sterr = self.execute_cmd(cmd)
                    for item in lines:
                        try:
                            name_item,energy_h, energy_rel, time_s = item.split()
                        except ValueError:
                            continue
                        self._db_base.update_data_row("qm_jobs", "energy", energy_h, "name_job", name_item)
                        self._db_base.update_data_row("qm_jobs", "cpu_time", time_s, "name_job", name_item)
            if os.path.exists(fulldir+name3):
                with open(fulldir+name3, 'r') as f:
                    lines = f.readlines()
                    cmd = 'cp %s %s_0_%d.txt'%(fulldir+name3, local_dir+name3.split(".")[0],i+1)
                    stdout, sterr = self.execute_cmd(cmd)
                    for item in lines:
                        name_item, state, name_output = item.split()
                        self._db_base.update_data_row("qm_jobs", "status_job", state, "name_job", name_item)
        # Pairs
        for i in range(2):
            for j in range(2):
                subdir = qm_engine + "_inputs_{:1d}_{:1d}/".format(i+1, j+1)
                name2 = "summary_energy.txt".format(i+1, j+1)
                name3 = "summary.txt".format(i+1, j+1)
                fulldir = remote_dir+subdir
                if os.path.exists(fulldir+name2):
                    with open(fulldir+name2, 'r') as f:
                        lines = f.readlines()
                        cmd = 'cp %s %s_%d_%d.txt' % (fulldir + name2, local_dir + name2.split(".")[0], i + 1, j + 1)
                        stdout, sterr = self.execute_cmd(cmd)
                        for item in lines:
                            try:
                                name_item,energy_h, energy_rel, time_s = item.split()
                            except ValueError:
                                continue
                            self._db_base.update_data_row("qm_jobs", "energy", energy_h, "name_job", name_item)
                            self._db_base.update_data_row("qm_jobs", "cpu_time", time_s, "name_job", name_item)
                if os.path.exists(fulldir+name3):
                    with open(fulldir+name3, 'r') as f:
                        lines = f.readlines()
                        cmd = 'cp %s %s_%d_%d.txt' % (fulldir + name3, local_dir + name3.split(".")[0], i + 1, j + 1)
                        stdout, sterr = self.execute_cmd(cmd)
                        for item in lines:
                            name_item, state, name_output = item.split()
                            self._db_base.update_data_row("qm_jobs", "status_job", state, "name_job", name_item)
        self._db_base.commit_db()

    # ===========================================================================================
    def extract_energy_calculations(self, local_dir, remote_dir, qm_engine="gaussian"):

        """
        Create a bash file to get all energy values from logs in the remote_dir.
        Files summary_energy.txt and summary.txt are created

        Args:
            local_dir (str): Path to store files in the local server.
            remote_dir (str): Path to store files in the remote server.
            qm_engine (str): Name of the QM package (example: gaussian)

        """

        scriptname = "check_local_dir.sh"

        # Generate the bash-script to extract the energy from gaussian calculations
        self.generate_bashscript_check_jobs(qm_engine, remote_dir, inputname=scriptname)

        # Run the script --> Extract energy to a file named: summary_energy.txt
        try:
            cmd1 = 'cd %s; bash ./%s'%(remote_dir, scriptname)
            stdout, stderr = self.execute_cmd(cmd1)
        except:
            pass

    # ===========================================================================================
    @staticmethod
    def generate_bashscript_check_jobs(qm_engine, localdir, inputname="check_remote_dir.sh"):

        """Check_remote_dir.sh

        Args:
            localdir:
            inputname:

        Returns:

        """
        if qm_engine.lower() == "gaussian":
            extract_energy = "E\("
            extract_time = "`egrep 'Elapsed' $idir/$output  | awk '{print $3*86400+$5*3600+$7*60+$9}'`\n"
            is_complete_calc = "Normal term"
        elif qm_engine.lower() == "nwchem":
            extract_energy = "Total DFT energy"
            extract_time="`egrep wall $idir/$output |tail -1|awk '{print $6}'|sed 's/.$//'`\n"
            is_complete_calc = "Total times"
        elif qm_engine.lower() == "gamess":
            is_complete_calc = "EXECUTION OF GAMESS TERMINATED"
            extract_energy = "TOTAL ENERGY ="
            extract_time = "`egrep 'TOTAL WALL CLOCK TIME' $idir/$output | tail -1 | awk '{print $5'}`"

        if localdir[-1] == "/":
            localfile = localdir + inputname
        else:
            localfile = localdir + "/" + inputname

        with open(localfile, 'w') as f:
            l='#!/bin/bash\n'
            l+='rm -f "summary.txt"\n'
            l+='rm -f "summary_energy.txt"\n'
            l+='input="jobs.txt"\n'
            l+='index=0\n'
            l+='if [[ -f "$input" ]]; then\n'
            l+='  while IFS= read -r line\n'
            l+='  do\n'
            l+='    input=`echo "$line" | awk \'{print $1}\'`\n'
            l+='    idir=`echo "$line" | awk \'{print $2}\'`\n'
            l+='    output=`echo "$line" | awk \'{print $3}\'`\n'
            l+='    if [[ -f "$idir/$output" && ! -z `egrep "{}" $idir/$output` ]]; then\n'.format(is_complete_calc)
            l+='\n'
            l+='       echo "$idir COMPLETED $output" >>"summary.txt"\n'
            if qm_engine.lower() == "gamess":
                l+='       en=`egrep "{}"  $idir/$output | awk \'{{print $4}}\'`\n'.format(extract_energy)
            else:
                l+='            en=`egrep "{}"  $idir/$output | awk \'{{print $5}}\'`\n'.format(extract_energy)
                l+='            enmp2=`egrep "UMP2" $idir/$output | awk \'{{print $6}}\'|tr \'D\' \'E\'`\n'
                l+='            [[ -z "$enmp2" ]] && en=$en || tmp=$enmp2\n'
                l+='            [[ ! -z "$tmp" ]] && en=$(printf "%.14f" $tmp)\n'
            l+='       time_s={}\n'.format(extract_time)
            l+='       if [ $index -eq 0 ]; then\n'
            l+='          e0=$en\n'
            l+='       fi\n'
            l+='       index=`echo $index+1| bc -l`\n'
            l+='       erel=`printf "%.3f" $(echo "scale=3; (($en)-($e0))*627.5095"|bc) `\n'
            l+='       echo "$idir $en $erel $time_s" >>summary_energy.txt\n'
            l+='     else\n'
            l+='       echo "$idir RUNNING $output" >>"summary.txt"\n'
            l+='\n'
            l+='     fi\n'
            l+='   done < "$input"\n'
            l+='\n'
            l+='   [[ ! -z `egrep "PENDING|RUNNING" summary.txt` ]]  || echo -n >done\n'
            l+='fi\n'
            l+=' #echo -n >done\n'

            txt =l
            f.writelines(txt)

    # ***********************************************************************************
    @staticmethod
    def generate_bashscript_send_local(qm_engine, remotedir, logger=None):

        """full_send.sh

        Args:
            localdir:

        Returns:

        """

        if qm_engine.lower() == "gaussian":
            extension = "com"
        elif qm_engine.lower() == "nwchem":
            extension = "nw"
        elif qm_engine.lower() == "gamess":
            extension = "inp"
        else:
            m = "ERROR: QM Engine {} is unknown.".format(qm_engine.lower())
            print(m) if logger is None else logger.error(m)
            exit()


        localfile = remotedir+"full_send.sh"

        l = '#!/bin/bash\n'
        l += '\n'
        l += 'TOTALJOBS=`ls -ld *_[0-9]*/*_[0-9]* |wc -l`\n'
        l += 'WK=`pwd`'
        l += '\n'
        l += 'for idir in `ls -d *_[0-9]*/*_[0-9]*`; do\n'
        l += '\n'
        l += '    cd $idir\n'
        l += '\n'
        l += '    if [[ ! -e ../jobs.txt ]]; then\n'
        l += '        echo -n >../jobs.txt\n'
        l += '    fi\n'
        l += '    inp=`ls *.{}`\n'.format(extension)
        l += '    out=$inp.log\n'
        l += '\n'
        l += '    b=`basename $idir`\n'
        l += '    echo "$inp $b $out" >>../jobs.txt\n'
        l += '    bash send.sh\n\n'
        l += '    cd $WK\n'
        l += 'done\n'
        txt=l

        with open(localfile, 'w') as f:
            f.writelines(txt)