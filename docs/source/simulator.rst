Simulator Object
----------------

**Overview**
============

Simulator object 


**Example**
===========

.. code-block::

    # This generates the example of the figure


**Attributes**
==============

    molecules: list of Segments
        The list of segments in the simulator. The number of segments in this list is unlimited. 
        However, to calculate :math:`{\chi}` parameter, only two segments must be given.     
    natoms: int
        Total number of atoms in the simulator
    box: Box
        A box object configurating the simulation box

**Methods**
===========

.. autosummary::
    :nosignatures:

    chiripa.Simulator.Simulator.add_molecule 
    chiripa.Simulator.Simulator.get_coordinates
    chiripa.Simulator.Simulator.get_elements
    chiripa.Simulator.Simulator.get_ith_molecule
    chiripa.Simulator.Simulator.get_natoms
    chiripa.Simulator.Simulator.get_nmolecules
    chiripa.Simulator.Simulator.get_unique_elements
    chiripa.Simulator.Simulator.remove_molecule
    chiripa.Simulator.Simulator.write_DFTBopt_simulator
    chiripa.Simulator.Simulator.write_PDB_simulator
    chiripa.Simulator.Simulator.write_gaussian_simulator
    chiripa.Simulator.Simulator.write_nwchem_simulator
    
**API**
=======

.. autoclass:: chiripa.Simulator.Simulator
    :special-members: __init__
    


