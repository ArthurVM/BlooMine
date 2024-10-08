import sys
from os import kill, path
from signal import alarm, signal, SIGALRM, SIGKILL
from subprocess import PIPE, Popen

from .utilities import *


class Error(Exception):
    """ Base class for exceptions in this module.
    """
    pass


class InputError(Error):
    """ Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message

        sys.stderr.write(f"INPUT ERROR : {self.expression}\n{self.message}")
        sys.exit(1)


class command():
    """ A class for handling shell executed commands
    TODO: clean up this class
    """
    def __init__(self, command, subprocessID):
        self.command = str(command)
        self.subprocessID = str(subprocessID)

    def run(self, timeout = -1):

        kill_tree = True

        class Alarm(Exception):
            pass

        def alarm_handler(signum, frame):
            raise Alarm

        p=Popen(self.command, shell = True, stdout = PIPE, stderr = PIPE, env = None)

        if timeout != -1:
            signal(SIGALRM, alarm_handler)
            alarm(timeout)

        try:
            stdout, stderr = p.communicate()
            if timeout != -1:
                alarm(0)

        except Alarm as a:
            pids = [p.pid]
            if kill_tree:
                pids.extend(self.get_process_children(p.pid))
            for pid in pids:
                # This is to avoid OSError: no such process in case process dies before getting to this line
                try:
                    kill(pid, SIGKILL)
                except OSError:
                    pass
            return -1, '', ''

        return p.returncode, stdout, stderr

    def run_comm(self, if_out_return, so=None, se=None, exit=True):
        """ Run the command with or without an exit code
        """
        print(f"COMMAND={self.command}")

        if so != None:
            print(f"{self.subprocessID}\nCOMMAND={self.command}\n", file=so, flush=True)
        if se != None:
            print(f"{self.subprocessID}\nCOMMAND={self.command}\n", file=se, flush=True)

        returncode, stdout, stderr = self.run(360000)

        ## capture stdout and stderr
        if so != None:
            print(f"{stdout.decode()}\n", file=so, flush=True)
        if se != None:
            print(f"{stderr.decode()}\n", file=se, flush=True)

        if returncode and stderr and exit:
            print(f"{self.command} failed. Please check the log files to work out what went wrong.")
            print(f"STDOUT :\n{stdout.decode()}\nSTDERR :\n{stderr.decode()}\n")
            sys.exit(1)

        if if_out_return:
            return returncode, stdout, stderr

    def run_comm_quiet(self, if_out_return, so=None, se=None, exit=True):
        """ Run the command quietly with or without an exit code
        """

        if so != None:
            print(f"{self.subprocessID}\nCOMMAND={self.command}\n", file=so, flush=True)
        if se != None:
            print(f"{self.subprocessID}\nCOMMAND={self.command}\n", file=se, flush=True)

        returncode, stdout, stderr = self.run(360000)

        ## capture stdout and stderr
        if so != None:
            print(f"{stdout.decode()}\n", file=so, flush=True)
        if se != None:
            print(f"{stderr.decode()}\n", file=se, flush=True)

        if returncode and stderr and exit:
            print(f"{self.command} failed. Please check the log files to work out what went wrong.")
            print(f"STDOUT :\n{stdout.decode()}\nSTDERR :\n{stderr.decode()}\n")
            sys.exit(1)

        if if_out_return:
            return returncode, stdout, stderr

    def get_process_children(self, pid):
        p = Popen(f"ps --no-headers -o pid --ppid {pid}", shell = True,
	         stdout = PIPE, stderr = PIPE)
        stdout, stderr = p.communicate()
        return [int(p) for p in stdout.split()]
