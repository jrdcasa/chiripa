from abc import ABC, abstractmethod

class Box(ABC):

    """

    .. warning:: This class cannot be directly instantiated  because is an abstract class.

    """

    __slots__ = ['units', 'Lx', 'Ly', 'Lz', 'xlo', 'ylo', 'zlo', 'xhi', 'yhi', 'zhi',
                 'alpha', 'beta', 'gamma', 'isorthonormal', 'a', 'b', 'c', 'typebox',
                 'cellGrid', 'nearby']

    def __init__(self, xhigh, yhigh, zhigh, xlow=0.0, ylow=0.0, zlow=0.0, alpha=90.0, beta=90.0, gamma=90.0, units="real"):

        """

        The simulation box has its origin at (xlo, ylo, zlo) point. The lattice vectors are defined as:

            * Orthogonal cubic box **a** =(xhi-xlo,0,0), **b** =(0, yhi-ylo,0), **c** =(0, 0,zhi-zlo), a=b=c, alpha=beta=gamma=90

        Args:

            xhigh: (type float) --> X Coordinate of the maximum vertex in distance units
            yhigh: (type float) --> Y Coordinate of the maximum vertex in distance units
            zhigh: (type float) --> Z Coordinate of the maximum vertex in distance units
            xlow: (type float, default=0.0) --> X coordinates of the box origin in distance units
            ylow: (type float, default=0.0) --> Y coordinates of the box origin in distance units
            zlow: (type float, default=0.0) --> Z coordinates of the box origin in distance units
            alpha: (type float, default=90.0) --> Angles between b and c in degrees
            beta: (type float, default=90.0) --> Angles between a and c in degrees
            gamma: (type float, default=90.0) --> Angles between a and b in degrees
            units: (type string, default="real") --> **real** or **reduced** units

        """

        super().__init__()

        if units.upper() == "REAL":
            self.units = "real"
            self.nameunits = "angstroms"
        elif units.upper() == "REDUCED":
            self.units = "reduced"
            self.nameunits = "sigma"
        else:
            print ("\nWARNING!!!!!: ")
            print ("{0:s} units are not reconigzed".format(units))
            print ("Allowed values are real or reduced")
            print ("Changing units parameter to real\n")
            self.units = "real"

        self.xlo = xlow
        self.ylo = ylow
        self.zlo = zlow

        self.xhi = xhigh
        self.yhi = yhigh
        self.zhi = zhigh

        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

        # Length of the edges in the simulation box
        self.Lx = self.xhi - self.xlo
        self.Ly = self.yhi - self.ylo
        self.Lz = self.zhi - self.zlo

        self.a = None
        self.b = None
        self.c = None

        self.typebox = None
        self.isorthonormal = None
        self.cellGrid = None
        self.nearby = None

    def __str__(self):

        """

            Returns the state of the attributtes of an instance. Override __str__ method

        """

        objstr = str(self.__repr__())+"\n"
        for key in self.__slots__:
            try:
                value = getattr(self,key)
                objstr +=  str(key) +": "+str(value).upper() + "\n"
            except AttributeError:
                objstr += str(key) +": NOT SET" + "\n"

        return objstr

    @abstractmethod
    def center_of_box(self):
        """
            Virtual method to be implemented in all derived classes
        """
        pass

    @abstractmethod
    def linked_cell(self, max_rc, skin):
        """
            Virtual method to be implemented in all derived classes
        """

        pass