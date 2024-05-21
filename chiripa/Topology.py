from chiripa.MolecularGraph import MolecularGraph
from chiripa.internal_coordinates import distance_array
from chiripa.atomic_data import element_cov_radius, maximal_valences
import numpy as np

"""

Reference 1:    "A rule-based algorithm for automatic bond type perception"
                Qian Zhang, Wei Zhang, Youyong Li, Junmei Wang, Liling Zhang and Tingjun Hou
                Journal of Cheminformatics 2012, 4:26
                https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-4-26
                
Reference 2:    "Automatic Perception of Organic Molecules Based on Essential Structural Information"
                 Yuan Zhao, Tiejun Cheng, and Renxiao Wang*
                 J. Chem. Inf. Model. 2007, 47, 1379-1385

Reference 3:    "A New Algorithm for Exhaustive Ring Perception in a Molecular Graph"
                Th. Hanser, Ph. Jauffret, and G. Kaufmann
                J. Chem. Inf. Comput. Sci. 1996, 36, 1146-1152
"""

class Topology(MolecularGraph):

    """
    This is a derived class of MolecularGraph. It specializes the molecular graph to a topology
    """

    __slots__ = ['_orderbonds', '_nringsCauchy']

    # #########################################################################
    def __init__(self, nvert=-1, listbonds = None, undirected=True):

        """
            Constructor. It calls the super class to create a molecular graph. The

            ``Parameters``:
                * **nvert** (type: int, default = -1) -->
                * **listbonds** (type: list, default = None) -->
                * **undirected** (type: boolean, default = True) -->

            ``Return``:
                * **None**

        """

        super().__init__(nvert=nvert, listbonds=listbonds, undirected=undirected)

        self._orderbonds = np.zeros([self._natoms, self._natoms], dtype=float)
        self._nringsCauchy = 0

    # #########################################################################
    def __copy__(self):

        t = Topology()
        t._natoms = self._natoms
        t._nringsCauchy = self._nringsCauchy
        t._undirected = self._undirected
        t._bonds = self._bonds[:]
        t._cycles = self._cycles[:]
        t._nmols = self._nmols[:]
        t._graphdict = self._graphdict.copy()
        t._orderbonds = self._orderbonds.copy()
        t._iatch = self._iatch.copy()

        return t

    # #########################################################################
    def __eq__(self, other):

        """
        Overrides equal method

        ``Parameters``:
            * **other** (type: Topology) -->
        """

        if other is None:
            return None

        res = True
        # Get both attributtes from the super and sub clasess because the use of __slots__
        keys = super().__slots__+self.__slots__
        for key in keys:
            #print(key, "self."+key)
            #print(getattr(self,key))

            if isinstance(getattr(self,key), np.ndarray):
                par = np.array_equal(getattr(self,key), getattr(other,key))
                res = res and par
            elif isinstance(getattr(self,key), Topology):
                par = self.__dict__[key] == other.__dict__[key]
                res = res and par
            else:
                par = getattr(self,key) == getattr(other,key)
                res = res and par

        return res

    # #########################################################################
    def guess_bonds_topology(self, coords, elements):

        """
        Given a set of atoms, it guess if a bond exists between two atoms

        """

        natoms = self._natoms

        if np.shape(coords)[0] != natoms:
            raise ValueError('Coord must have same natoms rows. Natoms: {0:d}, Coords: {1:d}'
                             .format(natoms, np.shape(coords)[0]))
        if np.shape(elements)[0] != natoms:
            raise ValueError('Element must have same natoms rows. Natoms: {0:d}, Elements: {1:d}'
                             .format(natoms, np.shape(elements)[0]))

        # Calculate the atom distance matrix
        dist, tmp1, tmp2, tmp3 = distance_array(coords, coords)
        # Set up the connectivity of the molecule
        self.detect_connectivity(dist, elements)

        # Cauchy formula to detect the number of rings in the molecule
        nsegments = len(self.get_forest())
        nbonds = len(self.get_allbonds())
        self._nringsCauchy = nbonds - self._natoms + nsegments

    # #########################################################################
    def guess_nringsCauchy(self):
        # Cauchy formula to detect the number of rings in the molecule
        nsegments = len(self.get_forest())
        nbonds = len(self.get_allbonds())
        self._nringsCauchy = nbonds - self._natoms + nsegments
        return self._nringsCauchy

    # #########################################################################
    def get_bonds_topologyCONNECTPDB(self, filePDB):

        filePDB.seek(0)

        for line in filePDB:
            if not line.startswith('CONECT'):
                continue
            # The lines containing only the label does not take into account
            #CONNECT (without numbers)
            if line.split()[1]:
                iatom = int(line.split()[1])
                for jatom in line.split()[2:]:
                    self.add_edge([iatom-1, int(jatom)-1])
                    self._orderbonds[iatom-1, int(jatom)-1] = 0
                    self._orderbonds[int(jatom)-1, iatom-1] = 0

    # #########################################################################
    def detect_connectivity(self, distances, elements, test_max_valence=True):

        # Identification of bonded atoms using the method proposed in Reference 1
        # based on the distances of atoms
        isbonded = lambda dij,ri, rj : 0.8 < dij < ri+rj+0.4

        for iatom in range(self._natoms):
            for jatom in range(iatom+1, self._natoms):
                d  = distances[iatom, jatom]
                r1 = element_cov_radius[elements[iatom]]
                r2 = element_cov_radius[elements[jatom]]
                if isbonded(d, r1, r2):
                    self.add_edge([iatom, jatom])
                    self._orderbonds[iatom, jatom] = 0
                    self._orderbonds[jatom, iatom] = 0

        if test_max_valence:
            self.check_atom_max_valence(distances, elements)

    # #########################################################################
    def check_atom_max_valence(self, distances, elements):

        # Check number of covalently connected neighbors.
        # If the number of neighbors is greater that the value
        # given in maximal valence dictionary, remove the longest distance.
        for iatom in range(0, self._natoms):
            e = elements[iatom]
            if e in maximal_valences.keys():
                neigh =  self.get_neighbours(iatom)
                n_neigh = len(neigh)

                # For each neighbour
                while n_neigh > maximal_valences[e]:
                    neigh =  self.get_neighbours(iatom)
                    max_dist = 0.0
                    iatom_max = -1
                    jatom_max = -1
                    for jatom in neigh:
                        dij = distances[iatom, jatom]
                        if dij > max_dist:
                            max_dist = dij
                            iatom_max = iatom
                            jatom_max = jatom
                    self.remove_edge([iatom_max, jatom_max])
                    n_neigh -= 1

















