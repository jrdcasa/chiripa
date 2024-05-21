import os
import shutil
import numpy as np
import chiripa as chi

class Simulator(object):

    """
    This class collects all information neccesary to make the simulation

    """

    # ***********************************************************************************
    def __init__(self, mol=(), box=None):

        """
        Constructor. Build a simulator object.

        .. warning::
            If the constructor is called with only one molecule, the correct way to do it is:\n
            **Simulator(mol=(s1,))** instead of **Simulator(mol=(s1))**.\n
            In the last case, the parameter is not interpreted as a tuple

        Args:
            mol (tuple): Molecules in the simulator
            box (Box): A box object.
        """

        self.molecules = []
        self.natoms = 0
        # First add the box and after the molecules
        self.box = None
        if not box is None:
            self.box = box

        if mol:
            for imol in mol:
                if imol is None: continue
                self.add_molecule(imol)
            self.get_natoms()

    # ***********************************************************************************
    def __str__(self):

        """ Returns the state of the attributtes of an instance"""

        objstr = str(self.__repr__())+"\n"
        for key in self.__dict__:
            try:
                value = getattr(self,key)
                objstr +=  str(key) +": "+str(value) + "\n"
            except AttributeError:
                objstr += str(key) +": NOT SET" + "\n"
        return objstr

    # ***********************************************************************************
    def add_molecule(self, m):

        """

        Add the ``Segment`` m to the ``Simulator``. If the parameter m is not an instance of ``Segment``
        nothing is added to the ``Simulator`` and it returns ``False``

        Args:
            m (Segment): A ``Segment`` to be added to the ``Simulator``

        Returns:
            ``True`` if the segment is succesfully added, otherwise ``False``

        """

        if not isinstance(m, chi.Segment):
            self.__warningerror_class(1)
            return False

        self.molecules.append(m)

        return True

    # ***********************************************************************************
    def get_coordinates(self):

        """Return the coordinates of the ``Simulator``

        Returns:
            A ndarray[natoms,3]: Coordinates of all atoms in the system.

        """

        index_mol = 0
        for imol in self.molecules:
            if index_mol == 0:
                coords_full = imol._coords
            else:
                coords_full = np.concatenate((coords_full, imol._coords),axis=0)
            # print(type(coords_full))
            index_mol +=1

        return coords_full

    # ***********************************************************************************
    def get_elements(self):
        """Get the elements of the ``Simulator``

        Returns:
            dict: Elements in the ``Simulator``. d[index]=element

        .. code-block::

            {0: 'C', 1: 'H', 2: 'H', 3: 'C', 4: 'H', 5: 'H', 6: 'H', 7: 'H', 8: 'C', 9: 'H', 10: 'H', 11: 'H', 12: 'C', 13: 'H', 14: 'H', 15: 'C', 16: 'H', 17: 'H', 18: 'H'}


        """

        elements = {}
        index = 0
        for imol in self.molecules:
            for i in imol._elements:
                elements[index] = i
                index += 1

        return elements

    # ***********************************************************************************
    def get_ith_molecule(self, index=0):
        """
        Get the ith molecule of the simulator

        Args:
            index (int) : Index of the molecule

        Returns:
            A segment instance or None if index is not correct

        """

        if not isinstance(index,int):
            return None

        if index >= len(self.molecules):
            return None
        else:
            return self.molecules[index]

    # ***********************************************************************************
    def get_natoms(self):
        """Get the number of atoms

        Returns:
            int: Number of atoms

        """
        self.natoms = 0
        for imol in self.molecules:
            self.natoms += imol._natoms

        return self.natoms

    # ***********************************************************************************
    def get_nmolecules(self):

        """Get the number of molecules of the simulator

        Returns:
            int: Number of segments in the ``Simulator``

        """

        return len(self.molecules)

    # ***********************************************************************************
    def get_unique_elements(self):
        """Get a list with the chemical elements present in the ``Simulator``

        Returns:
            A set with the name of the elements

        """

        elements = ()
        for imol in self.molecules:
            t = tuple(imol._elements)
            elements += t

        elements = set(elements)

        return elements

    # ***********************************************************************************
    def remove_molecule(self, m):

        """Remove the segment with index m from the Simulator

        Args:
            m (int): Index of the segment in the molecule list

        Return:
            True if the segment is removed from the simulator, otherwise False

        """

        if m < len(self.molecules):

            del self.molecules[m]
            return True
        else:
            return False

    # ***********************************************************************************
    def write_PDB_simulator(self,filenamePDB='simulator.pdb',
                            title='no title', box=None, writeconnect=True,
                            with_vmd_tcl=False, appendpdb=False):

        """Write the molecules in the Simulator object to a pdb file.

        If the PDB file exists, the file will be overwritten.

        Args:
            filenamePDB (str): Name of the pdb (or path)
            title (str): Title in the PDB file section
            box (Box): The Box information to be written in the PDB file CRYST section
            writeconnect(bool): If TRUE the "CONECT" section will be written
            with_vmd_tcl (str): Create a default (template) TCL script to load the PDB in VMD
            appendpdb (bool): Append the structure to an existing PDB file

        """

        if not box is None:
            self.box = box

        header="HEADER    Generated by Chipar\n"
        title="TITLE    {0:s}\n".format(title)
        if self.box is None:
            cryst="CRYST1{0:9.3f}{1:9.3f}{2:9.3f}{3:7.2f}{4:7.2f}{5:7.2f} {6:s}\n".\
                format(100.0, 100.0, 100.0, 90.0, 90.0, 90.0, 'P 1')
        else:
            cryst="CRYST1{0:9.3f}{1:9.3f}{2:9.3f}{3:7.2f}{4:7.2f}{5:7.2f} {6:s}\n".\
                format(self.box.Lx, self.box.Ly, self.box.Lz, self.box.alpha, self.box.beta, self.box.gamma, 'P 1')

        if appendpdb:
            label = 'a'
        else:
            label = 'w'

        with open(filenamePDB, label) as f:

            f.writelines(header)
            f.writelines(title)
            f.writelines(cryst)

            # ATOMS
            chain_identification=['A','B','C','D','E']
            imol = 0
            iatom_global = 0
            for iseg in self.molecules:
                imol += 1
                c = iseg.get_coords()
                j = imol%(len(chain_identification)) - 1
                for iatom in range(iseg._natoms):
                    iatom_global += 1
                    line ="ATOM  {0:5d} {1:4s}{2:1s}{3:3s} {4:1s}{5:4d}{6:1s}   {7:8.3f}{8:8.3f}{9:8.3f}{10:6.2f}{11:6.2f}\n".\
                        format(iatom_global, iseg._elements[iatom],' ','MOL',chain_identification[j],imol,' ',c[iatom,0],c[iatom,1],c[iatom,2],1.0,1.0)
                    f.writelines(line)

            # CONNECTIVITY
            offset = 0
            line = ""
            if writeconnect:
                for iseg in self.molecules:
                    for key, item in iseg._topology._graphdict.items():
                        if len(item) > 0:
                            line="CONECT{0:5d}".format(key+1+offset)
                            for ielem in item:
                                line += "{0:5d}".format(ielem+1+offset)
                        line += "\n"
                        f.writelines(line)
                    offset += iseg._natoms
            f.writelines("END\n")


        # Write TCL for vmd --> vmd -e <name of the tcl>
        onlyfilenamePDB = "./"+filenamePDB.split("/")[-1]
        filenameTCL = os.path.splitext(filenamePDB)[0]+".tcl"
        if with_vmd_tcl:
            with open(filenameTCL,'w') as f:

                loadpdb = "mol new {0:s} autobonds no waitfor all\n".format(onlyfilenamePDB)
                representation ="pbc box -center origin\n" \
                                 "mol delrep 0 top\n" \
                                 "mol selection all\n" \
                                 "mol representation CPK\n" \
                                 "mol addrep top\n" \
                                 "mol showrep top 0 on\n" \
                                 "mol modcolor 0 top \"Chain\"\n"

                f.writelines(loadpdb)
                f.writelines("\n")
                f.writelines(representation)

    # ***********************************************************************************
    def write_qminput_simulator(self, maindir, qm_keywords,
                                inputname, delete_ok=False,
                                nproc="1", mem=None, is_send_script=None,
                                partition=None, nodelist=None):


        """Write an input for a QM package.
         Each package needs an interface. QM packages implemented:

            * Gaussian16 (https://gaussian.com/)
            * NwChem 6.8 and 6.6 (https://www.nwchem-sw.org/)
            * Gamess version (https://www.msg.chem.iastate.edu/gamess/)

        Args:
            maindir (str): Directory to store the input file
            qm_keywords (dict): A dictionary with QM keywords to run the calculation
            delete_ok: (bool): If True the maindir is overwritten
            inputname (str): Name of the gaussian file
            mem (str): Memory used by gaussian
            nproc (str): Number of processors

        Returns:
            True if the input file is succesfully created

        """

        isok = True
        # Create the directory to run the calculation
        fullpath = os.path.abspath(maindir) # Full path to the directory to run Gaussian16
        try:
            os.makedirs(fullpath)
        except FileExistsError:
            if delete_ok:
                shutil.rmtree(fullpath)
                os.makedirs(fullpath)
            else:
                isok = False

        # Create the input file for Gaussian package
        elements_full = self.get_elements()
        coords = self.get_coordinates()

        # =========================== WRITE SCRIPTS ==============================
        # Write inputs
        if qm_keywords['qm_engine'] == "gaussian":
            finputname = fullpath + "/" + inputname + ".com"
            chi.gaussian_write_optm(finputname, elements_full,
                                    coords, qm_keywords,
                                    chkfullpath=inputname.split(".")[0]+".chk",
                                    nproc = nproc,
                                    mem=mem)
        elif qm_keywords['qm_engine'] == "nwchem":
            finputname = fullpath + "/" + inputname + ".nw"
            chi.nwchem_write_optm(finputname, elements_full,
                                  coords, qm_keywords,
                                  nproc=nproc,
                                  mem=mem)
        elif qm_keywords['qm_engine'] == "gamess":
            finputname = fullpath + "/" + inputname + ".inp"
            chi.gamess_write_optm(finputname, elements_full,
                                  coords, qm_keywords,
                                  nproc=nproc,
                                  mem=mem)
        else:
            m = "ERROR SIMULATOR WRITE INNPUT QM FOR THE COMBINATION: {} {}".format(is_send_script. qm_keywords['qm_engine'])
            print(m)
            exit()

        # =========================== SLURM SCRIPTS ==============================
        # Write send bash scripts --> To run the job in the server
        if is_send_script == "slurm" and qm_keywords['qm_engine'] == "gaussian":
            chi.gaussian_basic_slurm_script(maindir,
                                            inputname,
                                            qm_keywords['qm_path'],
                                            partition,
                                            nodelist=nodelist,
                                            jobname=inputname,
                                            cpuspertask=nproc, memory=mem)
        elif is_send_script == "slurm" and qm_keywords['qm_engine'] ==  "nwchem":
            chi.nwchem_basic_slurm_script(maindir,
                                          inputname,
                                          qm_keywords['qm_path'],
                                          partition,
                                          nodelist=nodelist,
                                          jobname=inputname,
                                          cpuspertask=nproc, memory=mem)
        elif is_send_script == "slurm" and qm_keywords['qm_engine'] ==  "gamess":
            chi.gamess_basic_slurm_script(maindir,
                                          inputname,
                                          qm_keywords['qm_path'],
                                          qm_keywords['scratch_dir'],
                                          partition,
                                          nodelist=nodelist,
                                          jobname=inputname,
                                          cpuspertask=nproc, memory=mem)
        elif is_send_script == "localhost" and qm_keywords['qm_engine'] == "gaussian":
            chi.gaussian_basic_local_script(maindir, inputname,
                                            qm_keywords['qm_path'],
                                            qm_keywords['scratch_dir'],)

        elif is_send_script == "localhost" and qm_keywords['qm_engine'] == "nwchem":
            chi.nwchem_basic_local_script(maindir, inputname,
                                            qm_keywords['qm_path'],
                                            cpuspertask=nproc,)

        elif is_send_script == "localhost" and qm_keywords['qm_engine'] == "gamess":
            chi.gamess_basic_local_script(maindir, inputname,
                                          qm_keywords['qm_path'],
                                          qm_keywords['scratch_dir'],
                                          cpuspertask=nproc,)
        else:
            m = "ERROR SIMULATOR WRITE SEND QM FOR tHE COMBINATION: {} {}".format(is_send_script, qm_keywords['qm_engine'])
            print(m)
            exit()

        return isok

    # ***********************************************************************************
    def __warningerror_class(self, msg, txt=""):

        if msg == 1:
            line = "\n\t======== WARNING ==========\n" \
                   "\tThe paramater passed to Simulator.add_molecule must be an Segment instance.\n" \
                   "\tMolecule is not added to the Simulator object!!!.\n " \
                   "\t======== WARNING ==========\n"
            print(line)
            print(txt)
            print(self)
