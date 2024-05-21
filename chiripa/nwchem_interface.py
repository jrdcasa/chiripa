import datetime

# ***********************************************************************************
def nwchem_write_optm(ffullname, dict_elements, coords,
                      nwchem_keywords, chkfullpath=None,
                      nproc=1, mem=None,
                      title="Input file generated by CHIRIPA",
                      charge=0, multiplicity=1):

    """

    Args:
        ffullname:
        dict_elements:
        coords:
        nwchem_keywords:
        chkfullpath:
        nproc:
        mem:
        title:
        charge:
        multiplicity:

    Returns:

    """

    # Check input parameters ******
    if coords.shape[0] != len(dict_elements):
        print("Warning!!!. The {} cannot be written. "
              "Number of atoms in element list and number of coordinates must be the same".format(ffullname))
        return False

    prefix = ffullname.split("/")[-1].split(".")[0]

    with open(ffullname,'w') as fin:

        # Write echo
        fin.writelines("echo\n\n")

        # Write start and title
        fin.writelines("start {}\n\n".format(prefix))
        fin.writelines("title \"{}\"\n".format(title+": "+prefix))
        fin.writelines("charge {}\n".format(nwchem_keywords['charge']))

        # Write coordinates
        fin.writelines("\ngeometry units angstroms print xyz autosym noautoz\n")
        # Coordinates ******
        index = 1
        for icoord in coords:
            iel = dict_elements[index-1]
            fin.writelines("{0:s} {1:>12.6f} {2:>12.6f} {3:>12.6f} \n".
                            format(iel, icoord[0], icoord[1], icoord[2]))
            index += 1
        fin.writelines("end\n")

        # scratch directory
        if not nwchem_keywords['scratch_dir'] is None:
            fin.writelines("\nscratch_dir {}\n".format(nwchem_keywords['scratch_dir']))
        fin.writelines("\n")

        # Basis set
        fin.writelines("basis\n")
        fin.writelines("\t* library {}\n".format(nwchem_keywords['basis_set']))
        fin.writelines("end\n")
        fin.writelines("\n")

        # DFT method
        fin.writelines("dft\n")
        fin.writelines("\txc {}\n".format(nwchem_keywords["qm_method"]))
        fin.writelines("\tNOIO\n")
        fin.writelines("\tDIRECT\n")
        fin.writelines("\tmult {}\n".format(nwchem_keywords["multiplicity"]))
        fin.writelines("\tgrid fine\n".format(nwchem_keywords["multiplicity"]))
        fin.writelines("end\n")
        fin.writelines("\n")

        # TASK method
        fin.writelines("driver\n")
        fin.writelines("\tmaxiter 100\n")
        fin.writelines("end\n")
        fin.writelines("\n")

        fin.writelines("task dft {}\n".format(nwchem_keywords["qm_task"]))

    return True

# ***********************************************************************************
def nwchem_basic_slurm_script(maindir,
                              inputname,
                              nwchempath,
                              partition,
                              nodelist=None,
                              numbernodes=1,
                              cpuspertask=1, memory=None, jobname=None,
                              fnamescript="send.sh"):

    """
    Write a script to run Nwchem in server with slurm queue system

    ``Parameters``::
        * **fnamescript** (type str) : Name of the script to run with sbatch
        * **g16path** (type str) : Path to the gaussian executable
        * **inputname** (type str) : name of the input file to run Gaussian
        * **partition** (type str) : name of the partition in the Slurm system
        * **nodelist** (type list) : list of nodes to run the job
        * **numbernodes** (type integer) : Number of nodes
        * **cpuspertask** (type integer) : Number of cores within the nodes
        * **memory** (type integer): Memory in Gigabytes
        * **jobname** (type str) : Name in the job list

    ``Returns``::
        * **None**

    """

    if jobname is None:
        jobname = "nwchem"

    with open(maindir+"/"+fnamescript, 'w') as f:

        f.writelines("#!/bin/bash\n")
        f.writelines("#SBATCH --partition={}\n".format(partition))
        if not nodelist is None:
            l = ""
            for item in nodelist:
                l+=item+", "
            l = l[:-1]
            f.writelines("#SBATCH --exclude=\"{}\"\n".format(l))
        f.writelines("#SBATCH -N 1\n")
        f.writelines("#SBATCH -n {}\n".format(cpuspertask))
        if not memory is None:
            f.writelines("#SBATCH --mem={}G\n".format(memory))
        f.writelines("#SBATCH --job-name={}\n".format(jobname))
        f.writelines("\n")
        f.writelines("WD=`pwd`\n")
        f.writelines("\n")
        f.writelines("cd $SLURM_SUBMIT_DIR\n")
        f.writelines("\n")
        f.writelines("export OMP_NUM_THREADS=1\n")
        f.writelines("\n")
        f.writelines("echo \"Job Started!!! `date`\" >time.dat\n")
        f.writelines("\n")
        f.writelines("NWCHEM={}\n".format(nwchempath))
        f.writelines("NPROC={}\n".format(cpuspertask))
        f.writelines("INP={}\n".format(inputname+".nw"))
        f.writelines("OUT={}\n".format(inputname+".log"))
        f.writelines("\n")
        f.writelines("mpirun -np $NPROC $NWCHEM $INP &>$OUT\n")
        f.writelines("\n")
        f.writelines("echo \"Job Done!!! `date`\" >>time.dat\n")

# ***********************************************************************************
def nwchem_basic_local_script(maindir, inputname, nwchempath,
                              cpuspertask, fnamescript="send.sh"):

    with open(maindir+"/"+fnamescript, 'w') as f:
        f.writelines("#!/bin/bash\n")
        f.writelines("export OMP_NUM_THREADS=1\n")
        f.writelines("echo \"Job Started!!! `date`\" >time.dat \n")
        f.writelines("\n")
        f.writelines("NWCHEM={}\n".format(nwchempath))
        f.writelines("NPROC={}\n".format(cpuspertask))
        f.writelines("INP={}\n".format(inputname+".nw"))
        f.writelines("OUT={}\n".format(inputname+".nw.log"))
        f.writelines("\n")
        f.writelines("mpirun -np $NPROC $NWCHEM $INP &>$OUT\n")
        f.writelines("\n")
        f.writelines("rm -f *.db *.movecs\n")
        f.writelines("echo \"Job Done!!! `date`\" >>time.dat\n")