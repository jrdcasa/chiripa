Cubic Box
---------
**Overview**
============

A generalization of the class ``Box`` for cubic boxes. 
The dimensions of a cubic box must be equal ( :math:`L_x=L_y=L_z`), thus the user has to be
consistent when the object be created.

**Example**
===========

.. code-block::

    import chiripa as chi

    # Creating a simple CubicBox of lenth 10 angstroms
    box=chi.CubicBox(10, units="real")
    box=chi.CubicBox(10)

    # Creating a simple CubicBox of lenth 10 sigma (in reduced units) 
    box=chi.CubicBox(10, units="reduced")

    # Raising an exception
    box2 = chi.CubicBox(20., 20., 20., xlow=5., ylow=5., zlow=10)
    # Lx = 15, Ly = 15, Lz = 10

    # Creating a box and the subcells for linked link method
    chi.init_logger("Output", fileoutput="out/test_box.log", append=False, inscreen=False)
    box_simple = chi.CubicBox(10, 10, 10)
    box_simple.linked_cell(1.0, 0.5, logger=log)



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
    * ``isorthonormal`` (boolean): Simulation box orthogonal or not. Always True for CubicBox
    * ``cellGrid`` (ndarray[3]): Number of subcells for each dimension, [nx, ny, nz]
    * ``nearby`` (ndarray[27,nx*ny*nz]): Neigbours array


**Methods**
===========

.. autosummary::
    :nosignatures:

    chiripa.CubicBox.CubicBox._define_regions
    chiripa.CubicBox.CubicBox.center_of_box 
    chiripa.CubicBox.CubicBox.linked_cell

**API**
=======

.. autoclass:: chiripa.CubicBox.CubicBox
    :members: 
    :show-inheritance:
    :private-members:
    :special-members: __init__
    

