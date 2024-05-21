Atomic_data
-----------

**Overview**
============

This module is like a database. It contains data about atomic information used 
in ``chiripa``

**Example**
===========

.. code-block::

    import chiripa as chi
    d1 = chi.atomic_mass["C"]
    print(d1)
    # print 12.01

    d3 = chi.element_cov_radius["N"]
    print(d3)
    # print 0.68 

    d2 = chi.element_vdw_truhlar_radius["C"]
    print(d2)
    # print 1.40

    # Calculation of the vdW volume 
    element_vdw_vmd_volume_bondi = dict()
    factor = 4.0 / 3.0 *  3.14159265359
    for key, val in element_vdw_vmd_radius_bondi.items():
        4.0 * Pi / 3.0 = 4.188790205
        element_vdw_vmd_volume_bondi[key] = round(factor * val*val*val,2)
    for item, v in element_vdw_vmd_volume_bondi.items():
        print("'{0:s}': {1:.02f},".format(item, v))

    v2 = chi.element_vdw_vmd_volume_bondi["C"]
    print(v2)
    # print 20.58

    d1 = chi.distances_revaluated_Okuwaki["c3-c3"]
    print(d1)
    # print 3.65


**API**
=======

.. automodule:: chiripa.atomic_data
    

