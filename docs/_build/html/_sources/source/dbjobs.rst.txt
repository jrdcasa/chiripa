DBJobs
------
**Overview**
============

This class is used to create a SQLite3 database. This database keeps track of all
QM jobs. 

.. image:: ../imgs/database.png
    :width: 600 px
    :align: center


**Example**
===========

.. code-block::

    # Create a database if it does not exist
    db1 = chi.DBjobs("out/test_02.db", logger=self.log)
    # Create a table named qmjobs
    success = db1.create_table_qmjobs()
    # Insert some data in the table
    sql_insert = "INSERT INTO qm_jobs ('id_job', 'name_job', 'cog') VALUES ({0:d}, '{1:s}', '{2:s}')". \
            format(db_index, inputname, str(cog_list)[0:-1])
    db1.insert_data(sql_insert)
    # Close connection 
    db1.close_connection()


**Methods**
===========

.. autosummary::
    :nosignatures:

    chiripa.DBJobs.DBjobs.check_for_table
    chiripa.DBJobs.DBjobs.close_connection
    chiripa.DBJobs.DBjobs.commit_db
    chiripa.DBJobs.DBjobs.create_table_qmjobs
    chiripa.DBJobs.DBjobs.insert_data
    chiripa.DBJobs.DBjobs.number_of_rows
    chiripa.DBJobs.DBjobs.remove_row
    chiripa.DBJobs.DBjobs.remove_table
    chiripa.DBJobs.DBjobs.sort_by_column
    chiripa.DBJobs.DBjobs.update_allcolumns_tosamevalue
    chiripa.DBJobs.DBjobs.update_data_row


**API**
=======

.. autoclass:: chiripa.DBJobs.DBjobs
    :members:
    :show-inheritance:
    :undoc-members:
    :special-members: __init__, __str__
    

