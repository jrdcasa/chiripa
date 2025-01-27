from abc import ABC, abstractmethod
from chiripa import DBjobs

class Server(ABC):

    """Abstract class for a Server object

    """

    # ===========================================================================================
    def __init__(self, nameserver, dbname, username=None,
                 keyfile=None, logger=None):

        """
        Server abstract class

        ``Parameters``::
            * **nameserver** (type: string): Name of the server (example: localhost, trueno.csic.es)
            * **dbname**: (type: string) : File name of the database
            * **username** (type: string, default=None): User name in the server or None for localhost
            * **keyfile** (type: string, default=None): Path to the ssh key file (more information about ssh keys can be found in
                this `web <https://serverpilot.io/docs/how-to-use-ssh-public-key-authentication/>`_ )
            * **logger** (type: logger, default=None): Logger instance

        """

        self._nameserver = nameserver
        self._username = username
        self._key_file = keyfile
        self._client = None

        self._db_index = 0

        if dbname is None:
            self._db_base = DBjobs("automaticdb.db", logger=logger)
        else:
            self._db_base = DBjobs(dbname, logger=logger)

    # ===========================================================================================
    @abstractmethod
    def execute_cmd(self):
        """
            Virtual method to be implemented in all derived classes
            Implementation on the server to execute a command (from the WORKING directory)
        """
        pass

    # ===========================================================================================
    @abstractmethod
    def qm_calc(self):
        """
            Virtual method to be implemented in all derived classes
        """

        pass

    # ===========================================================================================
    @abstractmethod
    def qm_allsend(self):
        """
            Virtual method to be implemented in all derived classes
        """

        pass

    # ===========================================================================================
    @abstractmethod
    def get_energy_from_calculations(self):
        """
            Virtual method to be implemented in all derived classes
        """

        pass

    # ===========================================================================================
    @abstractmethod
    def extract_energy_calculations(self):
        """
            Virtual method to be implemented in all derived classes
        """

        pass

    # ===========================================================================================
    @abstractmethod
    def get_logs(self):
        """
            Virtual method to be implemented in all derived classes
        """

        pass