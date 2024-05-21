import os
from copy import copy
from openbabel import openbabel
import numpy as np
import logging
import chiripa as chi

class Segment(object):

    """Class to represent molecules


    There are two ways to build a Segment object:

        1. Segment(filecoord='ethylene.pdb')
        2. Segment(natoms=2, xlist= [], ylist=[], zlist = [], elements = [])

    If **filecoord** is present, the paremeters **natoms**, **xlist**, **ylist**, **zlist** and **elementlist** are ignored.
    Otherwise, these parameters must be consistent.

    if **filetop** is not given, then the topology is guessed.

    if **filetypeatoms** is present, the types assignation is done according to the file.
    The order of the atoms must be the same that the order in the coordinates and topology files.
    The format of this file must be the following:

    .. code-block::

                <Number> <Atom> <type_of_atom>
                1 C c3
                2 C c3
                3 C c3
                ...
                14 H hc

    .. warning::
                Be careful when use np.transpose function. This function seems to change the order of the array to
                "Fortran-type" instead to C-Order. If use np.transpose you will use np.ascontiguousarray

                Example:

                a = np.tranpose(b) --> a in Fortran order irrespective of the order of C

                a = np.ascontiguousarray(np.transpose(b)) --> a in C order

                This issue is important when use mode="c" in pyx files for Cython

    .. warning::
                For big molecules (>1000 atoms) deactivate the guessing of topology (guesstopol=False)

    """

    # ***********************************************************************************
    def __init__(self, name, filecoord=None, filetop=None, filetypeatoms=None, natoms=0,
                 xlist=None, ylist=None, zlist=None, elementlist=None, guesstopol=True,
                 dummy_head = -1, dummy_tail = -1):

        """Constructor of a Segment object

        Args:
            name (str): Name of the segment
            filecoord (str): Path of the coordinates filecoord (Format: PDB, GRO, XYZ)
            filetop (str): Path of the topology (Format: PDB)
            filetypeatoms (str): Name of the file containing the matching between atoms and atomtypes.
                This is used mainly to assign the reevaluated distances
                by `Okuwaki_ et_ al_. (Table 2) <https://pubs.acs.org/doi/abs/10.1021/acs.jpcb.7b08461>`_
            natoms (int): Number of atoms
            xlist (list, float): Float numbers x component of the coordinates (in angstroms)
            ylist (list, float): y component of the coordinates (in angstroms)
            zlist (list, float): z component of the coordinates (in angstroms)
            elementlist (list, str): Element Symbol
            guesstopol (bool): If True activate the guessing of topology based in a distance algorithm.
            dummy_head (int): Index of the atom acting as dummy head atom to mimic polymer chain
            dummy_tail (int): Index of the atom acting as dummy tail atom to mimic polymer chain


        """
        self._name = name
        self._filecoord = filecoord
        self._filetop = filetop
        self._logger = logging.getLogger("Output") #""Segment", append=True, )
        self._typeelements = None

        if filecoord is not None:
            self.load_from_disk(filecoord)
        else:
            if xlist is None: xlist = []
            if ylist is None: ylist = []
            if zlist is None: zlist = []
            self._natoms = natoms
            self._coords = np.column_stack((np.asarray(xlist),
                                            np.asarray(ylist),
                                            np.asarray(zlist)))
            self._elements = np.asarray(elementlist, dtype=np.str)

        self.check_parameter_consistence()

        self._dummy_head_atom = dummy_head
        self._dummy_tail_atom = dummy_tail
        self._charge = 0
        self._isBOassigned = False

        if guesstopol:
            if filecoord is None and filetop is None:
                #self._topology = None
                self._topology = chi.Topology(nvert=self._natoms)
                if self._natoms == 0:
                    self._topology = None
                else:
                    self._topology.guess_bonds_topology(self._coords, self._elements)
            elif filetop is None:
                self._topology = chi.Topology(nvert=self._natoms)
                self._topology.guess_bonds_topology(self._coords, self._elements)
            else:
                self.set_topology_from_disk(filetop)
        else:
            self._topology = None

        self._filetypeatoms = filetypeatoms
        if filetypeatoms is not None:
            self.set_typeatoms(filetypeatoms)

    # ***********************************************************************************
    def __str__(self):

        """Returns the state of the attributtes of an instance"""

        objstr = str(self.__repr__())+"\n"
        for key in self.__dict__:
            try:
                value = getattr(self,key)
                objstr +=  str(key) +": "+str(value) + "\n"
            except AttributeError:
                objstr += str(key) +": NOT SET" + "\n"
        return objstr

    # ***********************************************************************************
    def __copy__(self):

        s = Segment(self._name, guesstopol=False)

        s._natoms = self._natoms
        s._topology = copy(self._topology)
        s._coords = self._coords.copy()
        s._elements = self._elements.copy()
        s._filecoord = self._filecoord
        s._filetop = self._filetop
        s._filetypeatoms  = self._filetypeatoms
        if self._typeelements is not None:
            s._typeelements = self._typeelements.copy()

        return s

    # ***********************************************************************************
    def __eq__(self, other):

        if other is None:
            return None

        res = True
        for key in self.__dict__:
            if isinstance(self.__dict__[key], np.ndarray):
                par = np.array_equal(self.__dict__[key], other.__dict__[key])
                res = res and par
            elif isinstance(self.__dict__[key],chi.Topology):
                par = self.__dict__[key] == other.__dict__[key]
                res = res and par
            elif isinstance(self.__dict__[key],chi.Segment):
                par = self.__dict__[key] == other.__dict__[key]
                res = res and par
            else:
                par = self.__dict__[key] == other.__dict__[key]
                res = res and par

        return res

    # ***********************************************************************************
    def assign_bond_orders(self):

        """This function assigns bond orders to the bonds according to the algorithm reported in:

        "Automated simultaneous assignment of bond orders and formal charges"
        Ivan D. Welsh and Jane R. Allison
        J. Cheminform (2019) 11:18

        https://doi.org/10.1186/s13321-019-0340-0

        The function uses the external software ``indigo-bondorders`` (located in thirdparty/indigo-bondorder).
        This code is compiled and installed in thirdparty/indigox


        .. warning::
            The structure to assign bonds needs to have all hydrogen bonds correctly placed.
            United atom models do not work with this function.

        """
        #sys.path.insert(0, "/home/jramos/PycharmProjects/CHIPAR_2/")
        import indigox as ix

        # Periodic Table data from indigox
        PT = ix.PeriodicTable()
        # Build a molecule in the indigox framework
        mol = ix.Molecule()
        mol.SetTotalCharge(self._charge)

        # Add all atoms in a dictionary and get the bonds in the
        # framework of indigox program
        all_atoms = dict()
        bonds_topo = self._topology.get_edges()
        for i, j in bonds_topo:
            if i not in all_atoms:
                # Element of i
                e = self._elements[i]
                all_atoms[i] = mol.NewAtom(PT[e])
                index = all_atoms[i].SetIndex(i)
                name = e+str(index)
                all_atoms[i].SetName(name)
            if j not in all_atoms:
                # Element of j
                e = self._elements[j]
                all_atoms[j] = mol.NewAtom(PT[e])
                index = all_atoms[j].SetIndex(j)
                name = e+str(index)
                all_atoms[j].SetName(name)

            mol.NewBond(all_atoms[i], all_atoms[j])

        # Setup to use the FPT algorithm with single electrons without preplacing
        # to calculate bond orders and formal charges
        opts = ix.Options.AssignElectrons
        opts.ALGORITHM = opts.Algorithm.FPT
        opts.FPT.ADD_EDGES_TO_TD = False
        opts.FPT.MINIMUM_PROPAGATION_DEPTH = 1
        opts.USE_ELECTRON_PAIRS = False

        # Calculate bond orders and formal charges.
        # Count have the total number of resonance structures
        nresonances = mol.AssignElectrons()
        #print("{} resonace structure(s) calculated with a score of {}.".format(nresonances, mol.GetMinimumElectronAssignmentScore()))

        # Sum all order bonds for the resonace structures.
        for iresonance in range(nresonances):

            mol.ApplyElectronAssignment(iresonance)

            for ibond in mol.GetBonds():
                i = ibond.GetSourceAtom().GetIndex()
                j = ibond.GetTargetAtom().GetIndex()
                bo =  ibond.GetOrder()
                if bo == bo.SINGLE_BOND:
                    self._topology._orderbonds[i,j] += 1.0
                    self._topology._orderbonds[j,i] += 1.0
                elif bo == bo.DOUBLE_BOND:
                    self._topology._orderbonds[i,j] += 2.0
                    self._topology._orderbonds[j,i] += 2.0
                elif bo == bo.TRIPLE_BOND:
                    self._topology._orderbonds[i,j] += 3.0
                    self._topology._orderbonds[j,i] += 3.0
                else:
                    print("=================================================")
                    print("Warning!!!! -> Bond order cannot be assigned "
                          "between {} and {} atoms".format(i,j))
                    print(bo)
                    print("=================================================")
                    self._isBOassigned = False

        # Correct for aromaticity
        nrows = self._topology._orderbonds.shape[0]
        ncols = self._topology._orderbonds.shape[1]
        for i in range(nrows):
            for j in range(i,ncols):
                m = self._topology._orderbonds[i,j] % nresonances
                if m != 0:
                    self._topology._orderbonds[i,j] = 1.5
                    self._topology._orderbonds[j,i] = 1.5
                else:
                    self._topology._orderbonds[i,j] /= nresonances
                    self._topology._orderbonds[j,i] /= nresonances

    # ***********************************************************************************
    def calc_vdw_volume_VABC(self):

        """Calculation of the van der Waals volume using the method reported by Zhao et al.

        "Fast Calculation of van der Waals Volume as a Sum of Atomic and
        Bond Contributions and Its Application to Drug Compounds",  J. Org. Chem. 2003, 68, 7368-7373
        https://pubs.acs.org/doi/10.1021/jo034808o.

        The VdW radii and volume are taken from
        element_vdw_vmd_radius_bondi and element_vdw_vmd_volume_bondi, respectively.

        The formula (4) of the article will be used in this function:

        .. image:: ../imgs/volume_vdw.png
            :width: 400pt


        Returns:
            (tuple): tuple containing:

                ``volume_vdw`` (float) Van der waals volume using equation 4 in (angstroms^3/molecule).
                ``volume_tsar`` (float) Van der waals volume using equation 6 in (angstroms^3/molecule).

        """

        s1 = 0.0
        for iatom in self._elements:
            s1 += chi.element_vdw_vmd_volume_bondi[iatom]

        mol = openbabel.OBMol()

        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "fix")


        obConversion.ReadFile(mol, self._filecoord)

        Rg = len(mol.GetSSSR())
        RA = sum([ 1 for i in mol.GetSSSR() if i.IsAromatic() is True])
        RNR = Rg - RA
        NB = self._natoms - 1 + Rg

        volume_vdw = s1 - 5.92*NB -14.7*RA -3.8*RNR

        volume_tsar = 0.801*volume_vdw + 0.18

        return volume_vdw, volume_tsar

    # ***********************************************************************************
    def center_of_geom(self, atomlist = None):

        """
        Calculate the center of geometry of the current coordinates

        Args:
            atomlist (list): Index of the atoms to calculate the center of geometry. If ``None`` all atoms are used

        Returns:
            ``cog``: ndarray[3] of float64. Coordinates of the ``Segment`` geometry center

        """

        if atomlist is None:
            natoms = self._natoms
        else:
            natoms = len(atomlist)

        tmp = np.zeros(3)
        c = self.get_coords(atomlist=atomlist)
        for iatom in range(natoms):
            tmp += c[iatom,:]

        cog = tmp/natoms
        return cog

    # ***********************************************************************************
    def center_of_mass(self, atomlist = None):

        """
        Calculate the center of mass of the current coordinates

        Args:
            atomlist (list): Index of the atoms to calculate the center of mass. If ``None`` all atoms are used

        Returns:
            ``com``: ndarray[3] of float64. Coordinates of the ``Segment`` center of mass

        """

        if atomlist is None:
            natoms = self._natoms
        else:
            natoms = len(atomlist)

        mtotal = 0.0
        tmp = np.zeros(3)
        c = self.get_coords(atomlist=atomlist)
        for iatom in range(natoms):
            if atomlist is None:
                m = chi.atomic_mass[self._elements[iatom]]
            else:
                m = chi.atomic_mass[self._elements[atomlist[iatom]]]
            mtotal += m
            tmp += c[iatom,:]*m

        com = tmp/mtotal
        return com

    # ***********************************************************************************
    def check_parameter_consistence(self):

        """It checks the length of the parameters passed through the constructor.

        The length of the x, y, z and element arrays must be equal to the number ot atoms.
        If there is not consistency raises a ValueError otherwise return True.

        Returns:
            True: if all parameters are consistent

        """

        condition = (self._coords.shape[0] == self._natoms)

        if not condition:
            try:
                line = "\n\t======== ERROR ==========\n" \
                   "\tCoordinates arrays must have equal length and equal to number of atoms\n" \
                   "\tLength coords: %d\n" \
                   "\tLength Elements: %d\n"\
                   "\tNumber of atoms: %d\n" \
                    "\tCoordfile: %s\n" \
                    "\tTopofile: %s\n" \
                   "\t======== ERROR ==========\n"%(self._coords.shape[0] , len(self._elements), self._natoms, self._filecoord, self._filetop)
            except:
                line = "\n\t======== ERROR ==========\n" \
                       "\tCoordinates arrays must have equal length and equal to number of atoms\n" \
                       "\t======== ERROR ==========\n"
            self._logger.error(line)
            raise ValueError ("Coordinates arrays must have equal length and equal to number of atoms")

        return True

    # ***********************************************************************************
    def euler_orientation(self, iseed=None):
        """

        New coordinates of the atoms accordingly to random Euler angles.
        There are many definitions of the Euler angles
        (see: https://en.wikipedia.org/wiki/Euler_angles)

        The definition here used is that given in:

        .. code-block::

            MATHEMATICAL METHODS FOR PHYSICISTS
            SEVENTH EDITION
            George B. Arfken, Hans J. Weber, Frank E. Harris
            pag: 140-142

        .. image:: ../imgs/euler_book.png

        Parameters:
            iseed (int): Seed for the pseudo-random number generator. If ``None`` the seed is created from the current system time
                otherwise a deterministic random data using the ``iseed`` value is expected

        Returns:
            ``euler``: , list of floats, :math:`{\\alpha}`, :math:`{\\beta}` and :math:`{\\gamma}` values for the Euler angles in radians

        """

        # Generate euler angles ========
        if iseed is None:
            euler = chi.generate_random_euler_angles()
        else:
            euler = chi.generate_random_euler_angles(seed=iseed)

        # Create rotation matrix
        S = chi.euler_rotation_matrix(euler)
        # Take the transpose of the coordinates
        C = np.ascontiguousarray(np.transpose(self.get_coords()))
        # Change the coordinates (R) dot-product (column vector of the coordinates)
        #print(ref.flags)
        self._coords = np.ascontiguousarray(np.transpose(np.dot(S,C)))

        return euler

    # ***********************************************************************************
    def get_coords(self, atomlist = None):

        """
        Get the coordinates of the segment

        Parameters:
            atomlist (list): Index of the atoms to return. If ``None`` all atoms are returned

        Returns:
            ``coords`` : ndarray[natoms, 3], natoms = self._natoms or len(atomlist)). Coordinates array

        """

        if atomlist is None:
            return self._coords
        else:
            tmp_coords = np.zeros((len(atomlist),3))
            i = 0
            for item in atomlist:
                tmp_coords[i] = self._coords[item]
                i += 1
            return tmp_coords

    # ***********************************************************************************
    def load_from_disk(self, path):

        """
        Wrapper to load the structure from disk. Allowed files are pdb, xyz and gro files.
        The format of the file is taken from the extension of the file. If this extension is unkown then
        the method raises an exception.

        Args:
            path (str): Path to the file in the disk

        """

        ext = os.path.splitext (path)[1]

        if ext == ".pdb":
            self.read_pdb_from_scratch(path)
        elif ext == '.xyz':
            self.read_xyz_from_scratch (path)
        elif ext == '.gro':
            self.read_gro_from_scratch(path)
        elif ext == ".sdf":
            self.read_sdf_coordtopo_from_scratch(path)
        else:
            self.printerror ("Unkown molecular format for file: %s" % path)
            raise Exception("Unkown molecular format for file: %s"%path)

    # ***********************************************************************************
    def printerror(self, msg1):

        """
        Print the error message msg1 to the logger

        Args:
            msg1 (str): Message error

        """

        self._logger.error(msg1)

    # ***********************************************************************************
    def read_gro_from_scratch(self, gro_path):

        """
        Read the coordinates from a gro file. It checks if the file exist

        Args:
            gro_path (str): Path to the file in the disk

        """

        if os.path.isfile(gro_path):
            f = open(gro_path)
        else:
            self.printerror(msg1="GRO file must be provided")
            raise Exception('chipar.segment.read_gro_from_scratch requires an existing file as argument')

        f.readline().strip()
        self._natoms = int(f.readline().strip())

        iline = 1
        #chainid = 1
        xlist = []
        ylist = []
        zlist = []
        elist = []
        while iline <= self._natoms:
            line = f.readline()
            #resid = int(line[0:5].strip())
            #resname = line[5:10].strip()
            atomname = line[10:15] #force field name
            #tag = int(line[15:20].strip())
            xlist.append(float(line[20:28].strip())*10)
            ylist.append(float(line[28:36].strip())*10)
            zlist.append(float(line[36:44].strip())*10)
            elist.append(atomname.strip())
            iline += 1

        self._coords = np.column_stack((np.asarray(xlist),
                                        np.asarray(ylist),
                                        np.asarray(zlist)))
        self._elements = np.asarray(elist, dtype=np.str)

        f.close()

    # ***********************************************************************************
    def read_pdb_from_scratch(self, pdb_path):

        """
        Read thecoordinates from a pdb file. It checks if the file exist

        Args:
            pdb_path (str): Path to the file in the disk

        """

        if os.path.isfile(pdb_path):
            f = open(pdb_path)
        else:
            f = None
            self.printerror(msg1="PDB file must exist\n {} does not exist".format(pdb_path))
            exit()

        self._natoms = 0
        xlist = []
        ylist = []
        zlist = []
        elist = []
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                #tag = int(line[6:11].strip())
                name = line[12:16].strip()
                #resname = line[17:20].strip()
                #chainid = line[21]
                #resid = line[22:26].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                elem = line[76:78].strip()
                if elem == "":
                    elist.append(name.capitalize())
                else:
                    elist.append(elem.capitalize())
                xlist.append(x)
                ylist.append(y)
                zlist.append(z)
                self._natoms += 1

        self._coords = np.column_stack((np.asarray(xlist),
                                        np.asarray(ylist),
                                        np.asarray(zlist)))

        self._elements = np.asarray(elist, dtype=np.str)

        f.close()

    # ***********************************************************************************
    def read_sdf_coordtopo_from_scratch(self, sdf_path):

        """
        Read coordinates and topology from a SDF file

        Args:
            sdf_path (str): Path to the sdf file

        """

        if os.path.isfile(sdf_path):
            f = open(sdf_path)
        else:
            self.printerror(msg1="SDF file must be provided")
            raise Exception('chipar.segment.read_sdf_coordtopo_from_scratch requires an existing file as argument')

        try:
            for n in range(3):
                next(f)
        except StopIteration:
            raise PysimmError('pysimm.system.read_mol requires either '
                              'file or string as argument')
        version = None
        line = next(f)
        self._natoms = int(line.split()[0])
        nbonds = int(line.split()[1])
        if len(line.split()) >= 3:
            version = line.split()[-1]

        self._topology = chi.Topology(nvert=self._natoms)

        xlist = []
        ylist = []
        zlist = []
        elist = []
        if version == 'V2000':
            for iatom in range(self._natoms):
                line = next(f)
                x, y, z, elem, tmp, charge = line.split()[:6]
                xlist.append(float(x))
                ylist.append(float(y))
                zlist.append(float(z))
                elist.append(elem)

            for n in range(nbonds):
                line = next(f)
                iatom, jatom, order = list(map(int, line.split()[:3]))
                self._topology.add_edge([iatom-1, jatom-1])
                self._topology._orderbonds[iatom-1, jatom-1] = order
                self._topology._orderbonds[jatom-1, iatom-1] = order

        elif version == 'V3000':
            next(f)
            line = next(f)
            self._natoms = int(line.split()[0])
            nbonds = int(line.split()[1])
            next(f)
            for iatom in range(self._natoms):
                line = next(f)
                idf, elem, x, y, z, charge = line.split()[2:8]
                xlist.append(x)
                ylist.append(y)
                zlist.append(z)
                elist.append(elem)

            next(f)
            next(f)

            for n in range(nbonds):
                line = next(f)
                idf, order, iatom, jatom = list(map(int, line.split()[2:6]))
                self.add_edge([iatom-1, jatom-1])
                self._orderbonds[iatom-1, jatom-1] = order
                self._orderbonds[jatom-1, iatom-1] = order

        self._coords = np.column_stack((np.asarray(xlist),
                                        np.asarray(ylist),
                                        np.asarray(zlist)))
        self._elements = np.asarray(elist, dtype=np.str)

        f.close()

    # ***********************************************************************************
    def read_topology_from_pdb(self, path):

        """
        Try to set up the topology reading a PDB file. The "CONECT" section is
        used to yield the connectivity of the molecule, if present. Otherwise, the
        bonds are guessed.

        Args:
            path (str): Path to the PDB file

        """

        if os.path.isfile(path):
            f = open(path)
        else:
            self.printerror(msg1="PDB file for topology must be provided")
            raise Exception('chipar.segment.read_topology_from_pdb requires an existing file as argument')

        isthereconnect = False
        for line in f:
            if line.startswith('CONECT'):
                l = line.split()
                i = int(l[1])
                self._topology.add_vertex(i-1)
                for jj in l[2:]:
                    j = int(jj)
                    if j < i: continue
                    self._topology.add_vertex(j-1)
                    self._topology.add_edge([i-1,j-1])
                    isthereconnect = True

        if not isthereconnect:
            self._topology.guess_bonds_topology(self._coords, self._elements)

        f.close()

        # DEBUG
        # print (self._topology)
        # self._topology.draw_graph(title="kk")

    # ***********************************************************************************
    def read_topology_from_xyz(self, path):

        """
        Try to set up the topology reading a XYZ file, the
        bonds are guessed.

        Args:
            path (str): Path to the XYZ file

        """

        if os.path.isfile(path):
            f = open(path)
        else:
            self.printerror(msg1="XYZ file for topology must be provided")
            raise Exception('chipar.segment.read_topology_from_pdb requires an existing file as argument')

        isthereconnect = False

        if not isthereconnect:
               self._topology.guess_bonds(self._coords, self._elements)

        f.close()

        # DEBUG
        # print (self._topology)
        # self._topology.draw_graph(title="kk")

    # ***********************************************************************************
    def read_xyz_from_scratch(self, xyz_path):

        """
        Read the coordinates from a xyz file. It checks if the file exist

        Args:
            xyz_path (str): Path to the file in the disk

        """

        if os.path.isfile(xyz_path):
            f = open(xyz_path)
        else:
            self.printerror(msg1="XYZ file must be provided")
            raise Exception('chipar.segment.read_xyz_from_scratch requires an existing file as argument')

        nparticles = int(f.readline().strip())

        self._natoms = 0
        xlist = []
        ylist = []
        zlist = []
        elist = []
        f.readline().strip()
        for _ in range(nparticles):
            elem, x, y, z = f.readline().split()
            xlist.append(float(x))
            ylist.append(float(y))
            zlist.append(float(z))
            elist.append(elem.capitalize())
            self._natoms += 1

        self._coords = np.column_stack((np.asarray(xlist),
                                        np.asarray(ylist),
                                        np.asarray(zlist)))
        self._elements = np.asarray(elist, dtype=np.str)

        assert int(nparticles == self._natoms), \
            "Number of particles in the header is different to the read atoms in the xyz file"
        f.close()

    # ***********************************************************************************
    def set_charge(self, charge):

        """
        Setup the total charge of the segment

        Args:
            charge (int): Total charge of the segment

        """

        self._charge = charge
        return charge

    # ***********************************************************************************
    def set_topology_from_disk(self, path):

        """
        Wrapper to load the topology from disk. Allowed files are pdb, xyz and sdf files.
        The format of the file is taken from the extension of the file. If this extension is unkown then
        the method raises an exception.

        Args:
            path (str): Path to the file in the disk

        """

        ext = os.path.splitext (path)[1]

        if ext == ".pdb":
            self._topology = chi.Topology(nvert=self._natoms)
            self.read_topology_from_pdb(path)
        elif ext == ".xyz":
            self._topology = chi.Topology(nvert=self._natoms)
            self.read_topology_from_xyz(path)
        elif ext == ".gro":
            self._topology = chi.Topology(nvert=self._natoms)
            self.read_topology_from_gro(path)
        elif ext == ".sdf":
            "The topology is already set up in the " \
            "read_sdf_coordtopo_from_scratch method"
            pass
        else:
            self.printerror ("Unkown topology format for file: %s" % path)
            raise Exception("Unkown topology format for file: %s"%path)

    # ***********************************************************************************
    def set_typeatoms(self, filetypeatoms):

        """
        Set the type of atoms in the segment (self._typeelements)

        Args:
            filetypeatoms (str): File containing the type of atoms


        .. code-block::

                <Number> <Atom> <type_of_atom>
                1 C c3
                2 C c3
                3 C c3
                ...
                14 H hc

        """

        typelist = []
        with open(filetypeatoms, 'r') as f:
            for line in f:
                i, el, types = line.split()
                typelist.append(types)

        self._typeelements = np.asarray(typelist, dtype=np.str)
        return None

    # ***********************************************************************************
    def translate_vector(self, v):

        """
        Translate the segment along the vector v.
        This function changes the coordinates of the segment

        .. image:: ../imgs/translation.png
            :width: 200pt

        The vector p represents the coordinates of each atom.

        Args:
            v: (type: a list or (1,3)-ndarray): The atoms are translated along this vector

        """

        c = self.get_coords()
        for iatom in range(self._natoms):
            c[iatom,:] += v

        return None

    # ***********************************************************************************
    def write_xyz(self, filenamexyz="segment.xyz", atomlist = None):

        """
        Write the segment to a XYZ file

        Args:
            filenamexyz (str): Name of the file
            atomlist (list): Index of the atoms to be written in the cxyz file. If ``None`` all atoms are written

        """

        c_tmp = self.get_coords(atomlist = atomlist)
        natoms = c_tmp.shape[0]

        with open(filenamexyz, 'w') as w:

            w.writelines("{}\n".format(natoms))
            w.writelines("Molecule {}\n".format(filenamexyz))

            for i in range(natoms):
                if atomlist is None:
                    elem = self._elements[i]
                else:
                    elem = self._elements[atomlist[i]]

                line = "{0:s} {1:.4f} {2:.4f} {3:.4f}\n".format(elem, c_tmp[i,0], c_tmp[i,1], c_tmp[i,2])
                w.writelines(line)

    # ***********************************************************************************
    def read_topology_from_gro(self, path):

        """
        Try to set up the topology reading a XYZ file, the
        bonds are guessed.

        Args:
            path (str): Path to the XYZ file

        """

        if os.path.isfile(path):
            f = open(path)
        else:
            self.printerror(msg1="GRO file for topology must be provided")
            raise Exception('chipar.segment.read_topology_from_gro requires an existing file as argument')

        isthereconnect = False

        if not isthereconnect:
               self._topology.guess_bonds(self._coords, self._elements)

        f.close()

        # DEBUG
        # print (self._topology)
        # self._topology.draw_graph(title="kk")

