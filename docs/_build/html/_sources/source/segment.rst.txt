Segment Object
--------------
**Overview**
============

Segment object is used to represent a molecular structure

    
A PET monomer represented as a Segment:

.. image:: ../imgs/PET_segment.png
    :width: 600 px
    :align: center
    :alt: Principle of angle bending

The attributes of the Segment object are summarized in the figure.

**Example**
===========

.. code-block::

    # This generates the example of the figure
    s1 = Segment(filecoord="../data/pet_1mon_aa.pdb")
    # This creates a empty segment
    s1 = Segment()
    # Creating a segment for an ethylene monomer 
    s1 = Segment(natoms=2, xlist=[0.0, -0.002], ylist=[0.765, -0.777], zlist=[0.0, 0.0], elementlist=['C','C'])



**Attributes**
==============

    _filecoord: str
        Path to the coordinate file used as input.
    _filetop: str
        Path to the topology file used as input. If is None, the topology is guessed.
    _filetypeatoms: str
        Path to a file containing the type atoms. (see above the format of this file)
    _logger: logger
        LOGGER TO DO
    _natoms: int
        Number of atoms.
    _coords: ndarray[natoms,3] type=float64
        Coordinates of the atoms
    _topology: topology
        Topology object
    _typeelements: ndarray[natoms] type=str
        Type of elements read from self._filetypeatoms
    _elements: ndarray[natoms] type=str
        Name of the elements
    _charge: int
        Net charge of the segment
    _isBOassigned: boolean
        If it is ``True`` the order bonds have been assigened using the algorithm
        reported in https://doi.org/10.1186/s13321-019-0340-0
    _dummy_head_atom: integer
        Index of the atom acting as dummy head atom to mimic polymer chain. The value is -1 when is not used.
    _dummy_tail_atom: integer
        Index of the atom acting as dummy tail atom to mimic polymer chain. The value is -1 when is not used.

**Methods**
===========

.. autosummary::
    :nosignatures:

    chiripa.Segment.Segment.assign_bond_orders 
    chiripa.Segment.Segment.calc_vdw_volume_VABC
    chiripa.Segment.Segment.center_of_geom
    chiripa.Segment.Segment.center_of_mass
    chiripa.Segment.Segment.check_parameter_consistence
    chiripa.Segment.Segment.euler_orientation
    chiripa.Segment.Segment.get_coords
    chiripa.Segment.Segment.load_from_disk
    chiripa.Segment.Segment.read_gro_from_scratch
    chiripa.Segment.Segment.read_pdb_from_scratch
    chiripa.Segment.Segment.read_sdf_coordtopo_from_scratch
    chiripa.Segment.Segment.read_topology_from_pdb
    chiripa.Segment.Segment.read_topology_from_xyz
    chiripa.Segment.Segment.read_xyz_from_scratch
    chiripa.Segment.Segment.set_charge
    chiripa.Segment.Segment.set_topology_from_disk
    chiripa.Segment.Segment.set_typeatoms
    chiripa.Segment.Segment.translate_vector
    chiripa.Segment.Segment.write_xyz
    
**API**
=======

.. automodule:: chiripa.Segment
    :synopsis: Definition of a Segment as a representation of molecules.
    :members:
    :undoc-members:
    :show-inheritance:
    :private-members:
    :special-members: __init__
    

