Server
------
**Overview**
============

This is an Abstract Class (ABC in python). It just provides an interface for Server subclasses. 
Server class represents an entry point to generate a subclass to send jobs for different
simulation packages (example: Gaussian, NWchem) and different type of servers (example: Slurm, localhost)

**Example**
===========

This class cannot be directly instantiated because is an abstract class

**Attributes**
==============
    * ``_nameserver`` (string) : Name of the server (example: localhost, trueno.csic.es).
    * ``_username`` (string) : User name in the server or None for localhost
    * ``_key_file`` (string) : Path to the ssh key file (more information about ssh keys can be found in this `web <https://serverpilot.io/docs/how-to-use-ssh-public-key-authentication/>`_ )
    * ``_client`` (SSHClient) : A SSHClient instance to handle the SSH connection (more in `paramiko <http://docs.paramiko.org/en/stable/api/client.html>`_)
    * ``_dbindex`` (int) : Index to know the last element added to the database
    * ``_db_base`` (DBJobs) : Instance to handle the database
      


**API**
=======

.. autoclass:: chiripa.Server.Server
    :members:
    :show-inheritance:
    :special-members: __init__, __str__
    

