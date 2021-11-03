"""
Messages which will be written to the console in response to an inconsistency.

Stores warning messages used thought the system to alert the user of an input
error or failure to run to completion, and an explaination.

"""

import sys

def query_yes_no(question, default="yes"):
    """
    Ask a yes/no question via raw_input() and return their answer.

    Args:
        question (str): a string that is presented to the user.
        default (str): the presumed answer if the user just hits <Enter>.
                It must be "yes" (the default), "no" or None (meaning
                an answer is required of the user).
            
    Returns:
        True for "yes" or False for "no".

    """

    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y]/n "
    elif default == "no":
        prompt = " y/[N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")

def exiting(msg):
    """
    Writes out the exit message to the settings file, then exits EVolve.

    Args:
        msg (str): The exit message which is causing EVolve to crash.
    """

    with open('outputs/settings.txt', "a") as settings_f:
        settings_f.write("\n")
        settings_f.write(f"EXIT_MESSAGE: {msg}")
        settings_f.write("\n")
    
    exit(msg)

def invalid_option():
    """
    Can't run Evolve with hydrogen escape but without fastchem.

    Provides option to run with fastchem.
    """

    print('Evolve cannot be run with the hydrogen escape option but without FastChem.')
    answer = query_yes_no('Would you like to continue the run including the FastChem model?', default='yes')

    if answer == True:
        return True
    elif answer == False:
        exiting('Exiting.\nPlease rerun without H escape if you do not wish to use FastChem.')
