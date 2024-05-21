Box
---
**Overview**
============

Define a simulation BOX. This is an Abstract Class.

**Example**
===========

This class cannot be directly instantiated because is an abstract class

**Attributes**
==============
    * ``units`` (string) : Allowed values are "real" or "reduced". Indicate the class of units to be used in the box definition.
    * ``nameunits`` (string) : For "real" --> "angstroms". For "reduced" --> "sigma".
    * ``a``  (ndarray(3)): Unit cell vector a
    * ``b``  (ndarray(3)): Unit cell vector b
    * ``c``  (ndarray(3)): Unit cell vector c
    * ``alpha`` (float): Axial angle :math:`\alpha`
    * ``beta``  (float): Axial angle :math:`\beta`
    * ``gamma`` (float): Axial angle :math:`\gamma`
    * ``xlo`` (float): X coordinates of the box origin in distance units
    * ``ylo`` (float): Y coordinates of the box origin in distance units
    * ``zlo`` (float): Z coordinates of the box origin in distance units
    * ``xhi`` (float): X coordinates of the box origin in distance units
    * ``yhi`` (float): Y coordinates of the box origin in distance units
    * ``zhi`` (float): Z coordinates of the box origin in distance units
    * ``Lx`` (float):  Length of the edge in the X-direction :math:`x_{hi}-x_{lo}`
    * ``Ly`` (float):  Length of the edge in the Y-direction :math:`y_{hi}-y_{lo}`
    * ``Lz`` (float):  Length of the edge in the Z-direction :math:`z_{hi}-z_{lo}`
    * ``typebox`` (string): "cubic".Only "cubic" is implemented in this version
    * ``isorthonormal`` (boolean): Simulation box orthogonal or not
    * ``cellGrid`` (ndarray[3]): Number of subcells for each dimension, [nx, ny, nz]
    * ``nearby`` (ndarray[27,nx*ny*nz]): Neigbours array


**API**
=======

.. autoclass:: chiripa.Box.Box
    :members:
    :show-inheritance:
    :special-members: __init__, __str__
    

