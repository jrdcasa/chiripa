from chiripa.Box import Box
import numpy as np
import sys

class CubicBox(Box):

    """
    Cubic Box class

    """

    # ***********************************************************************************
    def __init__(self, xhigh, yhigh, zhigh, xlow=0.0, ylow=0.0, zlow=0.0, units="real"):

        """Builds a cubic box

        The simulation box has its origin at (xlo, ylo, zlo) point. The lattice vectors are defined as:

            * Orthogonal cubic box:
                * **a** =(xhi-xlo,0,0)
                * **b** =(0, yhi-ylo,0)
                * **c** =(0, 0,zhi-zlo)
                * Lx=Ly=Lz
                * alpha=beta=gamma=90

        Args:
            xhigh: (type float) --> X Coordinate of the maximum vertex in distance units.
            yhigh: (type float) --> Y Coordinate of the maximum vertex in distance units.
            zhigh: (type float) --> Z Coordinate of the maximum vertex in distance units.
            xlow: (type float, default=0.0) --> X Coordinate of the box origin in distance units
            ylow: (type float, default=0.0) --> Y Coordinate of the box origin in distance units
            zlow: (type float, default=0.0) --> Z Coordinate of the box origin in distance units
            units: (type string, default="real") --> **real** or **reduced** units

        .. image:: ../imgs/simulation_box.png
            :height: 300px

        """
        Box.__init__(self, xhigh, yhigh, zhigh, xlow=xlow, ylow=ylow, zlow=zlow,
                     alpha=90.0, beta=90.0, gamma=90.0, units=units)

        if xlow >= xhigh or\
            ylow >= yhigh or\
            zlow >= zhigh:

            self.__error(1)

        if not (self.Lx == self.Ly == self.Lz):

            self.__error(2)

        self.a = np.array([self.Lx, 0.0, 0.0])
        self.b = np.array([0.0, self.Ly, 0.0])
        self.c = np.array([0.0, 0.0, self.Lz])

        self.isorthonormal = True
        self.typebox = "cubic"

    # ***********************************************************************************
    def center_of_box(self):

        """Calculate the geometrical center ot the simulation box

        Returns:
            An (1,3)-array with the coordinates of the geometric center

        """

        xc = self.xlo+(self.xhi-self.xlo)/2.0
        yc = self.ylo+(self.yhi-self.ylo)/2.0
        zc = self.zlo+(self.zhi-self.zlo)/2.0

        return np.array([xc,yc,zc])

    # ***********************************************************************************
    def linked_cell(self, max_rc, skin, logger=None):

        """This method initializes the linked cell grid for the cubic box.

        The simulation domain is divided into subcells with an edge length (rc+skin)
        greater than or equal to the maximum cut-off radius (max_rc) of the interaction
        to be computed.

        It changes the value of self.cellGrid attribute.

        Parameters:

            rc: (float) --> maximum cut-off radius.
            skin: (float) --> extra distance to create the subcells
            logger: (logger) --> A logger object to write

        .. image:: ../imgs/cutoff_linkedcell.png
            :height: 200px

        """

        if max_rc == 0.0:
            m1 = self.Lx
            m2 = self.Ly
            m3 = self.Lz
            max_rc = np.max([m1, m2, m3])

        nx = int(self.Lx / (max_rc + skin))
        ny = int(self.Ly / (max_rc + skin))
        nz = int(self.Lz / (max_rc + skin))
        if nx < 1 or ny < 1 or nz < 1:
            self.__error(3, {'nx':nx, 'ny':ny, 'nz':nz, 'rc':max_rc, 'skin':skin})

        rix = self.Lx/float(nx)
        riy = self.Ly/float(ny)
        riz = self.Lz/float(nz)

        self.cellGrid = np.array([nx, ny, nz])

        m =  "\t\t--------------------LINKED CELL INFO--------------------\n"
        m += "\t\tBox Length             : [{0:.3f}, {1:.3f}, {2:.3f}] {3:s}\n".format(self.Lx, self.Ly, self.Lz, self.nameunits)
        m += "\t\tMaximum cutoff distance: {0:.3f} {1:s}\n".format(max_rc, self.nameunits)
        m += "\t\tSkin distance          : {0:.3f} {1:s}\n".format(skin, self.nameunits)
        m += "\t\tNumber of subcells     : [{0:d},{1:d},{2:d}] cells\n".format(nx, ny, nz)
        m += "\t\tSubcell Lenghts        : [{0:.3f},{1:.3f},{2:.3f}] {3:s}\n".format(rix, riy, riz, self.nameunits)
        m += "\t\t--------------------------------------------------------\n"
        print(m) if logger is None else logger.info(m)

        self._define_regions()

        return self.cellGrid

    # ***********************************************************************************
    def _define_regions(self):

        """Setups the neighbours array (self.nearby) with periodic boundary conditions (PBC)

        The grid of cells starts in 0 and finnish in
        (self.cellGrid[0]*self.cellGrid[1]*self.cellGrid[2])-1

        Example:
            self.nearby[:,10] = All neighbors of the cell number 10

        .. image:: ../imgs/cell_grid_6x6x6.png
            :height: 200px

        """

        maxCell3  = self.cellGrid[0]*self.cellGrid[1]*self.cellGrid[2]
        self.nearby = np.zeros([27, maxCell3], dtype=np.int)

        # Find the set of cells that contain sites( = mers) that
        # interact with sites in cell icell.
        for iz in range(0, self.cellGrid[2]):
            for iy in range(0, self.cellGrid[1]):
                for ix in range(0, self.cellGrid[0]):
                    # Define the x,y,z address of central cell
                    icellx = ix
                    icelly = iy * self.cellGrid[0]
                    icellz = iz * self.cellGrid[0]*self.cellGrid[1]

                    # Combine to get the index of the central cell
                    icell = icellx + icelly + icellz

                    # Find the 27 nearest neighbor cells of icell
                    ith = 0
                    for jz in range(-1, 2):
                        for jy in range(-1, 2):
                            for jx in range(-1, 2):
                                # Define the x, y, z address of the jcell found by shifting
                                # the icell's x address by +/-1 or 0, its y address by
                                # + / -ncells or 0, and its z address by + / -ncells * ncells or 0
                                jcellx = icellx + jx
                                jcelly = icelly + jy * self.cellGrid[0]
                                jcellz = icellz + jz * self.cellGrid[0]*self.cellGrid[1]

                                #  Get all the minimum image cells
                                if (jcellx == self.cellGrid[0]):
                                    jcellx = 1
                                elif jcellx == -1:
                                    jcellx = self.cellGrid[0]-1

                                if jcelly == self.cellGrid[0]*self.cellGrid[1]:
                                    jcelly = 0
                                elif jcelly == -self.cellGrid[0]:
                                    jcelly = self.cellGrid[0]*(self.cellGrid[1]-1)

                                if (jcellz == self.cellGrid[0]*self.cellGrid[1]*self.cellGrid[2]):
                                    jcellz = 0
                                elif jcellz == -self.cellGrid[0]*self.cellGrid[1]:
                                    jcellz = self.cellGrid[0]*self.cellGrid[1]*(self.cellGrid[2]-1)

                                #  Combine to get the index of the interacting cell
                                # Number of elements in each of the 27 neighbour cells
                                #   --- nearby(0: 26, 0:maxcell3-1) ---
                                self.nearby[ith, icell] = jcellx + jcelly + jcellz
                                ith = ith + 1

    # ***********************************************************************************
    def __error(self, id_error, parameters=None):

        """Raises exceptions for the Cubic Box

        Args:
            id_error: (int) --> Id of the error
            parameters: (dict, optional) --> Parameters to print on the screen

        """

        if id_error == 1:

            print ("====== ERROR =====")
            print ("Bad setup for CUBIC BOX")
            print ("Lo coordinates cannot be greater than Hi coordinates.")
            print ("The simulation box cannot be created. The program STOPS.")
            print ("self.xlo: {0:.3f}, self.xhi: {1:.3f}".format(self.xlo, self.xhi))
            print ("self.ylo: {0:.3f}, self.yhi: {1:.3f}".format(self.ylo, self.yhi))
            print ("self.zlo: {0:.3f}, self.zhi: {1:.3f}".format(self.zlo, self.zhi))
            raise Exception("Lo coordinates cannot be greater than Hi coordinates.")

        elif id_error == 2:

            print ("====== ERROR =====")
            print ("Bad setup for CUBIC BOX")
            print ("The simulation box cannot be created. The program STOPS.")
            print ("self.Lx: {0:.2f}".format(self.Lx))
            print ("self.Ly: {0:.2f}".format(self.Ly))
            print ("self.Lz: {0:.2f}".format(self.Lz))
            print ("self.alpha: {0:.2f}".format(self.alpha))
            print ("self.beta: {0:.2f}".format(self.beta))
            print ("self.gamma: {0:.2f}".format(self.gamma))
            raise Exception("All lengths, Lx, Ly and Lz must be equals. Lx={0:.2f}, Ly={1:.2f}, Lz={2:.2f}".format(self.Lx,self.Ly,self.Lz))

        elif id_error == 3:

            print ("====== ERROR =====")
            print ("Bad setup for CUBIC BOX")
            print ("The cut-off + skin cannot be greater than the length of the cell.")
            for key, value in parameters.items():
                print(key, value)
            print ("Cell length: ",self.Lx, self.Ly, self.Lz)
            raise Exception("The cut-off + skin cannot be greater than the length of the cell.")


        else:

            print ("====== ERROR =====")
            print ("Bad setup for CUBIC BOX")
            print ("Unknown error.")
            print ("The simulation box cannot be created. The program STOPS.")
            raise Exception("The simulation box cannot be created. The program STOPS.")


