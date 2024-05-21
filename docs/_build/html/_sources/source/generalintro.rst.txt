====================
General introduction
====================

Flory-Huggins (FH) interaction parameter (:math:`{\chi}`) is included in the FH theory to
quantify the thermodynamic interaction energy between distinct components. This parameter
is an important ingredient in mesoscale simulations of polymers, usign DPD [#]_ (Dissipative Particle Dynamics) or MesoDyn [#]_ techniques. 

Recently, some authors have proposed the use of ab-initio methods to estimate this parameter between segments [#]_ or a molecular cluster [#]_

CHIRIPA is a Python framework that allows one to calculate the interaction parameter (:math:`{\chi}_{ij}`)
by determining the interactions between two monomers (or segments) at a
molecular level using QM and DFT methods. 
The methodology implemented closely follows the algorithms proposed by 
Blanco-Fan [#]_ , [#]_ and Okuwaki et al. [3]_ , [#]_

CHIRIPA stands for CHI inteRactIon PArameter. The word chiripa is a Spanish word that means “by luck” making reference to the use of random numbers in the implemented algorithms.

**References**

.. [#] Fraaije, J. G. E. M.; van Vlimmeren, B. A. C.; Maurits, N.M.;  Postma,  M.;  Evers,  O.  A.;  Hoffmann,  C.;  Altevogt,  P.;Goldbeck-Wood, G.J. Chem. Phys.1997,106, 4260-4269. https://aip.scitation.org/doi/10.1063/1.473129
.. [#] Hoogerbrugge, P. J. & Koelman, J.M.V.A. Simulation microscopic hydrodynamic phenomena with dissipative particle dynamics. Europhys. Lett. 19, 155–160 (1992). https://iopscience.iop.org/article/10.1209/0295-5075/19/3/001
.. [#] Okuwaki, K.; Mochizuki, Y.; Doi, H.; Ozawa, T. Fragment Molecular Orbital Based Parametrization Procedure for Mesoscopic Structure Prediction of Polymeric Materials. J. Phys. Chem. B 2018, 122 (1), 338–347. https://doi.org/10.1021/acs.jpcb.7b08461
.. [#] Sepehr F, Paddison SJ, Dissipative Particle Dynamics interaction parameters from ab initio calculations, Chemical Physics Letters, 645, 2016, 20-26 https://www.sciencedirect.com/science/article/pii/S0009261415009653?via%3Dihub
.. [#] Blanco, M. Molecular Silverware. I. General Solutions to Excluded Volume Constrained Problems. J. Comput. Chem. 1991, 12 (2), 237–247. https://doi.org/10.1002/jcc.540120214
.. [#] Fan, C. F.; Olafson, B. D.; Blanco, M.; Hsu, S. L. Application of Molecular Simulation To Derive Phase Diagrams of Binary Mixtures. Macromolecules 1992, 25 (14), 3667–3676. https://doi.org/10.1021/ma00040a010
.. [#] Okuwaki, K.; Doi, H.; Mochizuki, Y. An Automated Framework to Evaluate Effective Interactionparameters for Dissipative Particle Dynamics Simulations Basedon the Fragment Molecular Orbital (FMO) Method J. Comput. Chem. Japan 2018, 17 (2), 102–109. https://doi.org/10.2477/jccj.2017-0048.



