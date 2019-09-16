import os
import chimera
import re
from chimera import runCommand # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages
from chimera import selection
from chimera import dialogs
from DockPrep import prep
import time
import AddH
from DockPrep.prefs import prefs, defaults, INCOMPLETE_SC, MEMORIZED_SETTINGS

# Import related to reading from command line
from select import epoll, EPOLLIN
import subprocess

def read_with_timeout(fd, timeout__s):
    """ Reads from fd until there is no new data for at least timeout__s seconds.
    Requires: linux > 2.5.44
    """
    buf = []
    e = epoll()
    e.register(fd, EPOLLIN)
    while True:
        ret = e.poll(timeout__s)
        if not ret or ret[0][1] is not EPOLLIN:
            break
        buf.append(fd.read(1))
    return "".join(buf)

def read_reply_log():
    r = dialogs.find('reply')
    r_text = r.text.get('1.0', 'end')
    return (r, r_text)

def rc(command, subprocess = False, process = None, output = False, timeout = 0):
    if subprocess:
        process.stdin.write(command + "\n")
        process.stdin.flush()
        if output:
            response = read_with_timeout(process.stdout, timeout)
            return(response)
        else:
            return 0
    else:
        runCommand(command)
        if output:
            response = read_reply_log()
            return(response)
        else:
            return 0






def save_and_clear_reply_log(logfile_name):
    r, r_text = read_reply_log()
    with open(logfile_name, 'a') as logfile:
        logfile.write(r_text)
    r.Clear()

def clean_protein(proc = None):
    if proc == None:
        rc("select protein")
        rc("select invert")
        rc("delete selected")
        rc("select ligand")
        rc("delete selected")
        rc("select #0:@.B")
        rc("delete selected")
    else:
        # Clean out all non-peptide atoms and ligands
        rc("select protein", subprocess = True, process = proc)
        rc("select invert",  subprocess = True, process = proc)
        rc("delete selected",  subprocess = True, process = proc)
        rc("select ligand",  subprocess = True, process = proc)
        rc("delete selected",  subprocess = True, process = proc)

        # Delete all alternative locations for residues
        rc("select #0:@.B",  subprocess = True, process = proc)
        rc("delete selected",  subprocess = True, process = proc)

def generate_surface(proc = None, output = False, timeout = 0.5):
    if proc == None:
        rc("select protein", subprocess = False)
        rc("split", subprocess = False)
        rc("surface allComponents false", subprocess = False)
        return(0)

    else:
        rc("select protein",  subprocess = True, process = proc)
        rc("split",  subprocess = True, process = proc)

        # allComponents false excludes bubbles from surface computation
        response = rc("surface allComponents false",  
        subprocess = True, process = proc, output = output, timeout = timeout) 

        return(response)

def add_hydrogens_prep():
    # Add hydrogen atoms/charges, mutate MSE and other non-standard AAs to standard, and perform other basic DockPrep operations.
    models = chimera.openModels.list(modelTypes=[chimera.Molecule])
    prep(models, addHFunc=AddH.hbondAddHydrogens)


def extract_features(file_names, max_attempts = 5, 
    logfile_name = 'log__testChimeraFeatureExtraction.out'):

    # Start a Chimera prompt as a subprocess
    proc = subprocess.Popen(["chimera", "--nogui"],
        stdin = subprocess.PIPE,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE)

    for fn in file_names:
        attempt = 1

        while attempt <= max_attempts:
            try:
                print("Processing: " + fn)

                rc("open " + fn, subprocess = True, process = proc)
                rc("open " + fn, subprocess = False)
                #print(read_with_timeout(proc.stdout, 1))
                clean_protein(proc)
                clean_protein()

                # Clear the reply log before attempting to generate surface
                #save_and_clear_reply_log(logfile_name)

                # Try to generate surface, but skip this residue if it fails
                response = generate_surface(proc, True, 0.5)
                generate_surface()
                

                if 'connected surface components' not in response:
                    raise Exception("Surface computation failed")
                    

                # Generate surface through regular rc
                proc.stdin.write("~select\n")
                rc("~select")
                add_hydrogens_prep()




                print("Success: " + fn)
                proc.stdin.write("close all\n")

                # Save and clear the reply log before you continue
                #save_and_clear_reply_log(logfile_name)
                break

            except Exception as e:
                print("Exception: " + fn + " attempt: " \
                    + str(attempt) + " Message: " + str(e))
                attempt = attempt + 1
    return(proc.communicate()[0])

def extract_features_from_directory(directory_name, max_attempts = 5,
    logfile_name = 'log_test_ChimeraFeatureExtraction.out'):
    file_names = [directory_name + fn \
    for fn in os.listdir(directory_name) if fn.endswith(".pdb")]    

    extract_features(file_names, max_attempts, logfile_name)

extract_features_from_directory("/home/julian/")