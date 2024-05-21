================
Chiripa keywords
================

A number of keywords must be defined to run a ``chiripa`` simulation.
These keywords are stored in a dictionary and then passed to 
the object `Chi_Universe <chi_universe.html>`_ to run the :math:`{\chi}_{ij}` calculation.

The following keys are allowed in the dictionaty:

``"names"``: 
    (`list <https://docs.python.org/3.8/tutorial/datastructures.html#more-on-lists>`_, mandatory)

    A list containing the name of the two segments. 

    Example::

        d["names"] = ["n-hexane", "nitrobenzene"]

``"filecoords"``:
    (`list <https://docs.python.org/3.8/tutorial/datastructures.html#more-on-lists>`_ mandatory)  

    A 2-elements list containing the path of the coordinate files for each segment.
    Allowed formats are: 

        * `pdb <http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html>`_
        * `xyz <https://en.wikipedia.org/wiki/XYZ_file_format>`_
        * `gro <http://manual.gromacs.org/documentation/2020/reference-manual/file-formats.html?highlight=gro#gro>`_ 
        * `sdf <https://en.wikipedia.org/wiki/Chemical_table_file>`_ 

    Example::

        d["filecoords"] = ["./n-hexane.pdb", "./nitrobenzene.pdb"]

``"filetop"``
    (`list <https://docs.python.org/3.8/tutorial/datastructures.html#more-on-lists>`_, mandatory)

    A 2-elements list containing the path of the topology files for each segment. In the case of xyz, gro and pdb without 
    CONECT section, the topology is guessed using the algoritm in [#]_ .

        * `pdb <http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html>`_
        * `xyz <https://en.wikipedia.org/wiki/XYZ_file_format>`_
        * `gro <http://manual.gromacs.org/documentation/2020/reference-manual/file-formats.html?highlight=gro#gro>`_ 
        * `sdf <https://en.wikipedia.org/wiki/Chemical_table_file>`_ 

    Example::

        "filetop" = ["./n-hexane.pbd", "./nitrobenzene.pdb"]

.. warning::

    The order in ``"name"``, ``"filecoords"`` and ``"filetop"`` lists must be consistent. That is the element 0(1) in 
    these lists must be correspond to the Segment 1(2) in the simulator 

``"coordination_numbers_Z"``
    (`bool <https://docs.python.org/3.8/c-api/bool.html>`_, optional, default=False)

    If ``true`` the coordination numbers (Z) of all pairs are calculated. 
    Keep in mind, that if the file **Z_results.log** exists in the working directory the calculation does not carried out.
    This avoid to calculate the Z number each time that you want to check the state of the QM calculations in the server.

    Example::

        "coordination_numbers_Z" = True


``"Z_parameters"``
    (`dict <https://docs.python.org/3.8/tutorial/datastructures.html#dictionaries>`_, optional)

    The components of this dictionary are:
        
        * ``Z_samples`` : (int, default=20) Number of samples to average.
        * ``Z_puttrialsmonomers`` : (int, default=1000) Number maximum of trials to put the segment 2 (screening segment)
          around the segment 1 (base segment)
        * ``Z_debug`` : (bool, default=False) Keep debug information or not
        * ``Z_nonbonded`` : (str, default="truhlar") van der Waals radii set to be used in the overlap algorithm. Allowed values are:

            * "truhlar" [#]_
            * "okuwaki_correction" [#]_

``"calculate_volume"``
    (`bool <https://docs.python.org/3.8/c-api/bool.html>`_, optional, default=False)

    Specifies if the segment(monomer) volume is to be calculated. The calculation of 
    the van der Waals volume using the method reported by Zhao et al. [#]_

        "Fast Calculation of van der Waals Volume as a Sum of Atomic and
        Bond Contributions and Its Application to Drug Compounds",  J. Org. Chem. 2003, 68, 7368-7373
        https://pubs.acs.org/doi/10.1021/jo034808o.

``"interaction_energy""``
    (`bool <https://docs.python.org/3.8/c-api/bool.html>`_, optional, default=False)

    Specifies if the interaction energy between segment pairs is to be calculated using a QM package.
    If True, the options defining the calculation should be given in ``"energy_parameters"`` dictionary.

``"energy_parameters"``
    (`dict <https://docs.python.org/3.8/tutorial/datastructures.html#dictionaries>`_, mandatory)
    
    The components of this dictionary are:

        * ``"qm_engine"`` : (str, default="gaussian"). QM package to perform the calculations of the interaction energies between pairs. Allowed values are "gaussian" and "nwchem"
        * ``"qm_path_exe"`` : (str, mandatory). Full Path to the executable file of the QM package.
        * ``"qm_scratch_dir"`` : (str, default="./"). Path to the scratch directory. If it does not exist the working directory ("./") is used.
        * ``"qm_charge"`` : (int, mandatory). Total charge of the system.
        * ``"qm_multiplicity"`` (int, mandatory). Spin multiplicity of the system.
        * ``"qm_basisset"`` (str, default="3-21g"). Basis set to be used
        * ``"qm_method"`` (str, default="hf"). Level of theory to be used.
        * ``"qm_task"`` (str, default="energy"). Task to perform. Allowed values are "energy" and "opt"
        * ``"number_configurations"`` (int, default=4). Number of configurations to be created for each pair (i,j)



.. [#] Zhang Q, Zhang W, Li Y, Wang J, Zhang L and  Hou T
       "A rule-based algorithm for automatic bond type perception"
       Journal of Cheminformatics 2012, 4:26
       https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-4-26
.. [#] Mantina et al. J. Phys. Chem. A, Vol. 113, No. 19, 2009, 5806-5812 (Table 12) (https://pubs.acs.org/doi/10.1021/jp8111556)               
.. [#] Okuwaki, K.; Mochizuki, Y.; Doi, H.; Ozawa, T. (Table 1 and 2) J. Phys. Chem. B 2018, 122 (1), 338–347. https://doi.org/10.1021/acs.jpcb.7b08461 
.. [#]  Zhao YH, Abraham MH, and Zissimos AM "Fast Calculation of van der Waals Volume as a Sum of Atomic and
        Bond Contributions and Its Application to Drug Compounds",  J. Org. Chem. 2003, 68, 7368-7373
        https://pubs.acs.org/doi/10.1021/jo034808o.

