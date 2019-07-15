import shlex
import subprocess


def run_command(command):
    '''
    This function runs a command using subprocess
    :param command: the command to run
    :return: the output value
    '''

    process = subprocess.Popen(shlex.split(command))
    out, err = process.communicate()

    if err is not None:
        raise Exception(err)

    return out

