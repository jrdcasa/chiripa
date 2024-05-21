import paramiko
import os
import datetime
from socket import gaierror
import chiripa as chi

"""
Some functions are inspired or directly copy from slurmqueen project

(https://github.com/Kasekopf/SlurmQueen)

"""

class ServerSlurm(chi.Server):

    __version__ = "1.0.0"

    # ===========================================================================================
    def __init__(self, nameserver,  databasename, username, key_file, logger=None):

        """
        Initialize a server based on SLURM queue system

        Args:
            nameserver (str): The name of the server (i.e: trueno.csic.es or localhost).
            databasename (str):  File name of the database.
            username (str): User name in the server
            key_file (str): Name of the key file generated with "ssh-keygen -t rsa". \
            The passphrase must be empty and the public part must be installed in the server \
            (see `web <https://serverpilot.io/docs/how-to-use-ssh-public-key-authentication/>`_)


        """

        super().__init__(nameserver, username=username,
                         keyfile=key_file, dbname=databasename)

        if not self._db_base is None:
            self._db_base.create_table_qmjobs()

    # ===========================================================================================
    def connection(self, logger=None):
        """
        Open a connection to this server. Throws an exception if the connection fails.

        Returns:
            A `SSHClient`_ connection to this server

        .. _SSHClient: http://docs.paramiko.org/en/stable/api/client.html
        """

        if not self.is_connected():
            self._client = paramiko.SSHClient()
            self._client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            key = paramiko.RSAKey.from_private_key_file(self._key_file)

            try:
                self._client.connect(self._nameserver, username=self._username, pkey=key)

            except gaierror:
                self._client = None
                raise Exception("Unable to connect to server " + self._nameserver)

            fmt = "%a %d/%m/%Y, %H:%M:%S"
            m = "\t#QM# Connected to " + self._nameserver + " at {}\n".format(datetime.datetime.now().strftime(fmt))
            print(m) if logger is None else logger.info(m)

        return self._client

    # ===========================================================================================
    def is_connected(self):
        """
        Check if we are actively connected to the server.

        Returns:
            bool: True if we have a connection open to the server, and False otherwise.

        """
        transport = self._client.get_transport() if self._client else None
        return transport and transport.is_active()

    # ===========================================================================================
    def close_connection(self):
        """
        Check if we are actively connected to the server.

        Returns:
            bool: True if we have a connection open to the server, and False otherwise.

        """
        self._client.close()

    # ===========================================================================================
    def ftp_connect(self):

        """
        Open an sftp connection to this server. Throws an exception if the connection fails.

        Returns:
             An `sftp`_ connection to this server.

        .. _sftp: http://docs.paramiko.org/en/stable/api/client.html
        """
        return self.connection().open_sftp()

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
        #stdin, stdout, stderr = self.connection().exec_command(command, timeout=timeout) # Non-blocking call
        stdin, stdout, stderr = self._client.exec_command(command, timeout=timeout) # Non-blocking call
        exit_status = stdout.channel.recv_exit_status()   # Blocking call

        if exit_status != 0:
            print("Command {} cannot be run in server {}".format(command, self._nameserver ))

        if other_input is not None:
            stdin.write(other_input)
            stdin.flush()

        return stdout.read().decode("utf8"), stderr.read().decode("utf8")

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
        Prepare the data to send Gaussian16 jobs

        #. Create the database if it does not exist
        #. Define the number of nodes used in the calculations
        #. For each pair segement replica:

            a. Write qm inputs
            b. Write script to send the calculation
            c. Insert data in the database

        #. Zip all qm jobs.
        #. FTP qm inputs.
        #. Unzip the remote file.

        Args:
            iseg (int): Index for the Base Segment
            jseg (int): Index for the Screening Segment
            ncalc (int): Number of replicas for the pair iseg-jseg.
            qm_keywords (str): TODO
            simulator (Simulator): A :doc:`Simulator <./simulator>` instance.
            localdir (str): Path to store files in the local server.
            remotedir (str): Path to store files in the remote server.
            qmpath (str): Path the to the qm executable.
            scratch_dir (str): Path to the scratch dir used for gaussian
            cpuspertask (int): Number of cpus to use in the gaussian calculation
            nodelist (list): List of nodes in the remote server. **(default=None)**.
            nodelistmaster(list): Node to run the general script
            partition (str): Name of the partition to run the calculations
            partitionmaster (str): Name of the partition run the general script
            localzipfile (str): Path for a zip file of the inputs **(default="gaussian_inputs.zip")**
            memory (str): Memory (example: 500Mb)
            logger (Logger): A :doc:`Logger <./setup_logger>` instance

        Examples:
            Use example

        """

        file_base = localzipfile.split("/")[-1].split(".")[0]
        cog_list = list()

        # Create the database
        if self._db_base is None:
            self._db_base = DBjobs("{}_sp.db".format(qm_keywords['qm_engine']), logger=logger)
            self._db_base.create_table_qmjobs()

        # Directory to store the inputs
        prev_dir = os.getcwd()
        os.chdir(localdir)

        #print("1:", datetime.datetime.now().time())

        # Create gaussian and script files in the local area ===================
        inode = 0
        try:
            nnodes = len(nodelist)
        except (AttributeError, TypeError):
            nnodes = 1

        fmt = "%H:%M:%S"
        m = "\t#QM# Prepare, Zip, Send and Unzip QM inputs for {}-{} segments ( {} )".format(iseg, jseg, datetime.datetime.now().strftime(fmt))
        print(m) if logger is None else logger.info(m)

        for iconf in range(ncalc):

            if ncalc > 1:
                chi.gen_conf_fan_algorithm_c(simulator, nonbondedEvaluation)

            cog_list = simulator.get_ith_molecule(-1).center_of_geom()
            inputname = "{0:s}_{1:05d}_{2:1d}_{3:1d}".format(qm_keywords['qm_engine'], iconf, iseg, jseg)
            idir = localdir+file_base+"/{}".format(inputname)

            try:
                n = nodelist
            except IndexError:
                n = None
            simulator.write_qminput_simulator(idir, qm_keywords,
                                              inputname,
                                              delete_ok=True,
                                              mem=memory,
                                              nproc=cpuspertask,
                                              is_send_script="slurm",
                                              partition=partition,
                                              nodelist=n)

            # Insert data in the gaussian database
            sql_insert = "INSERT INTO qm_jobs VALUES ({0:d}, '{1:s}', 'CREATED', 'NULL', 'NULL', '{2:s}', 'NULL')".\
                                                      format(self._db_index, inputname, str(cog_list)[1:-1])
            try:
                self._db_base.insert_data(sql_insert)
            except:
                pass
            self._db_index += 1
            try:
                inode = (inode+1) % nnodes
            except:
                pass
        self._db_base.commit_db()

        # ================ Compressing input files ===================
        # Zip all gaussian jobs
        chi.compress_files(dirbase=localdir+file_base, zipname=localzipfile)

        # FTP gaussian inputs *********************************
        server_sftp = self.ftp_connect()

        try:
            server_sftp.put(localdir+localzipfile, remotedir+localzipfile)
        except FileNotFoundError:
             m = "\nEither file {} does not exist in localhost\n".format(localdir + localzipfile)
             m += "Or directory {} does not exist in the remote server\n".format(remotedir)
             print(m) if logger is None else logger.error(m)
             os.remove("{}_sp.db".format(qm_keywords['qm_engine']))
             exit()
        server_sftp.close()

        d = {'status_job' : 'INSERVER'}
        self._db_base.update_allcolumns_tosamevalue("qm_jobs", d)
        self._db_base.commit_db()

        # Unzip the remote file *********************************
        # "unzip -n" --> never overwrite existing files
        stdout, stderr = self.execute_cmd("cd {}; unzip -n {}; pwd".
                                          format(remotedir,
                                          localzipfile.split(".")[0],
                                          remotedir+localzipfile))

        m = "\t#QM# Finnished at {}\n".format(datetime.datetime.now().strftime(fmt))
        print(m) if logger is None else logger.info(m)

        os.chdir(prev_dir)

    # ===========================================================================================
    def qm_allsend(self, localdir, remotedir,
                   partitionmaster=None,
                   nodelistmaster=None,
                   jobname="all_sh", logger=None, maxjobsslurm=60):

        """
        Send all qm jobs to a SLURM server.

        Args:
            localdir (str): Path to store files in the local server.
            remotedir (str): Path to store files in the remote server.
            partitionmaster (str): Name of the node to send the script to send the gaussian jobs to the queue server
            nodelistmaster (List): List containing the name of the nodes where gaussian calculation will run
            jobname (str): Name of the job in SLURM syst
            logger (Logger): A :doc:`Logger <./setup_logger>` instance


        """

        fmt = "%H:%M:%S"
        m = "\t#QM# Send Jobs to remote server ( {} )".format( datetime.datetime.now().strftime(fmt))
        print(m) if logger is None else logger.info(m)


        self.generate_bashscript_send_slurm(localdir, remotedir, maxjobsslurm=maxjobsslurm)

        # Put the script in the server to send the jobs
        server_sftp = self.ftp_connect()
        try:
            server_sftp.put(localdir+"full_send.sh", remotedir+"full_send.sh")
        except:
            m = "\nEither file {} does not exist in localhost\n".format(localdir +"full_send.sh")
            m += "Or directory {} does not exist in the remote server\n".format(remotedir)
            print(m) if logger is None else logger.error(m)
            os.remove("{}_sp.db".format(qm_keywords['qm_engine']))
            exit()

        server_sftp.close()

        if nodelistmaster is None:
            cmd1="cd %s; sbatch --partition=%s --job-name=%s  full_send.sh;"\
                 %(remotedir, partitionmaster, jobname)
        else:
            cmd1="cd %s; sbatch --partition=%s --nodelist=%s --job-name=%s  full_send.sh;"\
                 %(remotedir, partitionmaster, nodelistmaster, jobname)

        self.execute_cmd(cmd1)

        m = "\t#QM# Finnished at     {} ".format( datetime.datetime.now().strftime(fmt))
        print(m) if logger is None else logger.info(m)

        # Get the number of jobs by FTP *********************************
        # trueno_sftp = self.ftp_connect()
        # trueno_sftp.get(remotedir+"jobs.txt", localdir+"jobs.txt")
        # trueno_sftp.close()
        # #
        # # list_pid_jobs = [int(item.split()[-1]) for item in stdout.split("\n") if len(item)>10]
        # dict_pid_idjob = {}
        # with open(localdir+"jobs.txt", 'r') as f:
        #     lines = f.readlines()
        #     for l in lines:
        #         label1 = int(l.split()[0])
        #         label2 = int(l.split()[1])
        #         label3 = l.split()[2]
        #         dict_pid_idjob[label1] = [label2, label3]
        #
        # return dict_pid_idjob

    # ===========================================================================================
    def get_logs(self, local_dir, remote_dir):

        zipfilename = "logs.zip"
        #os.mkdir(local_dir+"gaussian_logs")

        cmd1 = 'cd %s; zip -R %s \'*.log\''%(remote_dir,zipfilename)
        stdout, stderr = self.execute_cmd(cmd1)

        try:
            server_sftp = self.ftp_connect()
            server_sftp.get(remote_dir+zipfilename, local_dir+zipfilename)
            server_sftp.close()
        except:
            pass

        cmd2 = 'cd %s; rm -f %s'%(remote_dir,zipfilename)
        stdout, stderr = self.execute_cmd(cmd2)

    # ===========================================================================================
    def get_energy_from_calculations(self, qm_engine, local_dir, remote_dir):

        """
        Get the summary files and update the database with the energies from log files.
        This is for SLURM servers

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

        # Get files from FTP
        server_sftp = self.ftp_connect()
        # Monomers
        for i in range(2):
                name1 = "{:s}_inputs_{:1d}_{:1d}/".format(qm_engine, 0, i+1)
                name2 = "summary_energy_{:1d}_{:1d}.txt".format(0, i+1)
                name3 = "summary_{:1d}_{:1d}.txt".format(0, i+1)
                try:
                    server_sftp.get(remote_dir+name1+energyfilename, local_dir+name2)
                except FileNotFoundError:
                    pass
                try:
                    server_sftp.get(remote_dir+name1+summaryfilename, local_dir+name3)
                except FileNotFoundError:
                    pass
        # Pairs
        for i in range(2):
            for j in range(2):
                name1 = "{:s}_inputs_{:1d}_{:1d}/".format(qm_engine, i+1, j+1)
                name2 = "summary_energy_{:1d}_{:1d}.txt".format(i+1, j+1)
                name3 = "summary_{:1d}_{:1d}.txt".format(i+1, j+1)
                try:
                    server_sftp.get(remote_dir+name1+energyfilename, local_dir+name2)
                except FileNotFoundError:
                    pass
                try:
                    server_sftp.get(remote_dir+name1+summaryfilename, local_dir+name3)
                except FileNotFoundError:
                    pass
        server_sftp.close()

        # Update the database
        # Monomers
        for i in range(2):
                name2 = "summary_energy_{:1d}_{:1d}.txt".format(0, i+1)
                name3 = "summary_{:1d}_{:1d}.txt".format(0, i+1)
                if os.path.exists(local_dir+name2):
                    with open(local_dir+name2, 'r') as f:
                        lines = f.readlines()
                        for item in lines:
                            try:
                                name_item, pid_slurm, energy_h, energy_rel, time_s = item.split()
                            except ValueError:
                                continue
                            self._db_base.update_data_row("qm_jobs", "pid_slurm", pid_slurm, "name_job", name_item)
                            self._db_base.update_data_row("qm_jobs", "energy", energy_h, "name_job", name_item)
                            self._db_base.update_data_row("qm_jobs", "cpu_time", time_s, "name_job", name_item)
                if os.path.exists(local_dir+name3):
                    with open(local_dir+name3, 'r') as f:
                        lines = f.readlines()
                        for item in lines:
                            pid_slurm, name_item, state = item.split()
                            self._db_base.update_data_row("qm_jobs", "status_job", state, "name_job", name_item)

        # Pairs
        for i in range(2):
            for j in range(2):
                name2 = "summary_energy_{:1d}_{:1d}.txt".format(i+1, j+1)
                name3 = "summary_{:1d}_{:1d}.txt".format(i+1, j+1)
                if os.path.exists(local_dir+name2):
                    with open(local_dir+name2, 'r') as f:
                        lines = f.readlines()
                        for item in lines:
                            try:
                                name_item, pid_slurm, energy_h, energy_rel, time_s = item.split()
                            except ValueError:
                                continue
                            self._db_base.update_data_row("qm_jobs", "pid_slurm", pid_slurm, "name_job", name_item)
                            self._db_base.update_data_row("qm_jobs", "energy", energy_h, "name_job", name_item)
                            self._db_base.update_data_row("qm_jobs", "cpu_time", time_s, "name_job", name_item)

                if os.path.exists(local_dir+name3):
                    with open(local_dir+name3, 'r') as f:
                        lines = f.readlines()
                        for item in lines:
                            try:
                                pid_slurm, name_item, state = item.split()
                                self._db_base.update_data_row("qm_jobs", "status_job", state, "name_job", name_item)
                            except ValueError:
                                pass

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

        scriptname = "check_remote_dir.sh"

        # Generate the bash-script to extract the energy from qm calculations in slurm
        self.generate_bashscript_check_jobs(qm_engine, local_dir, inputname=scriptname)

        # Transfer the script to the remote machine
        try:
            server_sftp = self.ftp_connect()
            server_sftp.put(local_dir+scriptname, remote_dir+scriptname)
            server_sftp.close()
        except:
            pass

        # Run the script --> Extract energy to a file named: summary_energy.txt
        try:
            cmd1 = 'cd %s; bash %s'%(remote_dir, scriptname)
            stdout, stderr = self.execute_cmd(cmd1)
        except:
            pass

    # ===========================================================================================
    @staticmethod
    def generate_bashscript_check_jobs(qm_engine, localdir, inputname="check_remote_dir.sh"):

        """Generate a bash script to check the jobs in the server. The name of the bash script is given by
        the parameter inputname

        Args:
            qm_engine (str):Name of the QM package (gaussian or nwchem)
            localdir (str): Path to store files in the local server.
            inputname (str):Name of the script

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
            l = '#!/bin/bash\n'
            l+='rm -f "summary.txt"\n'
            l+='rm -f "summary_energy.txt"\n'
            l+='input="jobs.txt"\n'
            l+='index=0\n'
            l+='if [[ -f "$input" ]]; then\n'
            l+='    while IFS= read -r line\n'
            l+='    do\n'
            l+='        pid=`echo "$line" | awk \'{print $1}\'`\n'
            l+='        idir=`echo "$line" | awk \'{print $2}\'`\n'
            l+='        output=`echo "$line" | awk \'{print $3}\'`\n'
            l+='        p=`sacct -j $pid --format=JobID,JobName%20,State| egrep -v "batch|----|JobID" | head -1`\n'
            l+='        echo $p >>"summary.txt"\n'
            l+='        if [[ ! -z `echo $p | egrep "COMPLETED"` ]]; then\n'
            if qm_engine.lower() == "gamess":
                l+='        en=`egrep "{}"  $idir/$output | awk \'{{print $4}}\'`\n'.format(extract_energy)
            else:
                l+='            en=`egrep "{}"  $idir/$output | awk \'{{print $5}}\'`\n'.format(extract_energy)
                l+='            enmp2=`egrep "UMP2" $idir/$output | awk \'{{print $6}}\'|tr \'D\' \'E\'`\n'
                l+='            [[ -z "$enmp2" ]] && en=$en || tmp=$enmp2\n'
                l+='            [[ ! -z "$tmp" ]] && en=$(printf "%.14f" $tmp)\n'
            l+='            time_s={}\n'.format(extract_time)
            l+='            if [ $index -eq 0 ]; then\n'
            l+='                 e0=$en\n'
            l+='            fi\n'
            l+='            index=`echo $index+1| bc -l`\n'
            l+='            erel=`printf "%.3f" $(echo "scale=3; (($en)-($e0))*627.5095"|bc) `\n'
            l+='            echo "$idir $pid $en $erel $time_s" >>summary_energy.txt\n'
            l+='         #echo "$idir $pid $en $erel"\n'
            l+='        fi\n'
            l+='    done < "$input"\n'
            l+='\n'
            l+='    [[ ! -z `egrep "PENDING|RUNNING" summary.txt` ]]  || echo -n >done\n'
            l+='fi\n'
            l+='#echo -n >done\n'
            txt =l
            f.writelines(txt)

    #===========================================================================================
    @staticmethod
    def generate_bashscript_send_slurm(localdir, remotedir, maxjobsslurm=50):

        """Generate the script **full_send.sh** in order to send jobs to a SLURM server

        Args:
            localdir (str): Path to store files in the local server.
            remotedir (str): Path to store files in the remote server.
            maxjobsslurm (str): Number of maximum jobs send to slurm

        """

        localfile = localdir+"full_send.sh"
        with open(localfile, 'w') as f:
            l='#!/bin/bash\n'
            l+='\n'
            l+='MAXJOBSINSLURM={}\n'.format(maxjobsslurm)
            l+='\n'
            l+='NJOBS=`squeue -h |wc -l`\n'
            l+='\n'
            l+='while true; do\n'
            l+='    for idir in `ls -d *_[0-9]*/*_[0-9]*`; do\n'
            l+='        echo $idir, ${NJOBS}, ${JOBSEND}, ${TOTALJOBS}\n'
            l+='        cd $idir\n'
            l+='        base=`basename $idir`\n'
            l+='\n'
            l+='        if [[ ! -e ../jobs.txt ]]; then\n'
            l+='            echo -n >../jobs.txt\n'
            l+='        fi\n'
            l+='\n'
            l+='        if [[ $NJOBS -lt $MAXJOBSINSLURM  ]]; then\n'
            l+='            if [[  -z `egrep $base ../jobs.txt` ]]; then\n'
            l+='                sbatch send.sh 1>tmp.txt;\n'
            l+='                jobid=`awk \'{print $NF}\' tmp.txt`\n'
            l+='                echo "${jobid} ${base} ${base}.log" >>../jobs.txt\n'
            l+='                rm tmp.txt\n'
            l+='            fi\n'
            l+='        fi\n'
            l+='\n'
            l+='        NJOBS=`squeue -h |wc -l`\n'
            l += '        #echo $NJOBS\n'
            l += '        JOBSEND=`cat ../jobs.txt | wc -l`\n'
            l += '        TOTALJOBS=`ls -ld ../*_[0-9]* |wc -l`\n'
            l += '        cd ../..\n'
            l += '        echo "JOBSEND: ${JOBSEND}, TOTALJOBS: ${TOTALJOBS}"\n'
            l+='\n'
            l += '    done\n'
            l += '    if [[ ${JOBSEND} -ge ${TOTALJOBS} ]]; then\n'
            l += '        break\n'
            l += '    fi\n'
            l += '    # Each 60 seconds checks the jobs\n'
            l += '    sleep 60\n'
            l+='done\n'

            txt = l

            f.writelines(txt)