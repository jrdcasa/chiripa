Molecular Graph
---------------
**Overview**
============

Molecular Graph object is used to represent a molecular structure

**Example**
===========
 
**Attributes**
==============
    _natoms: int
        Number of atoms (= number of nodes in the graph).
    _nmols: int
        Number of molecules (= number of non-conected graphs).
    _bonds: list of sets
        Bonds of the molecule (edges of the graph) [{0,1},...{10,11}]
    _cycles: list
        XXXXXX
    _graphdict: dict
        Dictionary to represent the graph. Example: 0:[1,2] means that node 0 is linked to 1 and 2.
    _iatch: ndarray, np.shape(natoms,), type=int
        Molecule (or chain) for each atom.
    _undirected: bool
        Graph undirected or not.

**Methods**
===========


**API**
=======

.. autoclass:: chiripa.MolecularGraph.MolecularGraph
    :special-members: __init__

    

