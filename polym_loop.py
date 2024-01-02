#!/usr/bin/env python

################################################################################
#
# polym_loop.py
# This file is part of the Polymatic distribution.
#
#
# Description: Controls the simulated polymerization loop of the Polymatic
# algorithm. Polymerization steps are performed in cycles. After each bond is
# made, an energy minimization is performed in LAMMPS. After each cycle, a
# molecular dynamics step is performed in LAMMPS.
#
# Syntax:
#  ./polym_loop.py
#
# User parameters and file paths necessary for the polymerization should be
# specified at the beginning of the script.
#
################################################################################
#
# Polymatic: a general simulated polymerization algorithm
# Copyright (C) 2013, 2015 Lauren J. Abbott
#
# Polymatic is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# Polymatic is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License (COPYING)
# along with Polymatic. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

import sys
import os
import subprocess
import shutil
import glob

#
# User Setup
#

# Parameters
bonds = 0
bonds_tot = 16
bonds_cyc = 5
md_cyc = 3
md_max = 100
keep = 0

# Input Scripts
input_polym = 'polym.in'
input_min = 'min.in'
input_md0 = 'md0.in'
input_md1 = 'md1.in'
input_md2 = 'md2.in'

# Polymatic Scripts
script_step = 'polym.pl'
script_init = 'polym_init.pl'
script_final = 'polym_final.pl'

# LAMMPS
#lmps = 'lmp_serial -in /file/path/to/lammps-*/src/lmp_serial'
lmps = 'mpirun -np 8 /usr/local/bin/lmp_mpi'


# Main
#

def main():
    print_header()
    polym_loop()
    print_footer()

#
# Functions
#

def polym_loop():

    # Variables
    global bonds

    # Initialization
    polym_init()

    # Step 0
    setup_step()
    em()

    # Steps 1-N
    while (bonds < bonds_tot):

        # Polymerization step
        bonds += 1
        setup_step()
        code = polym_step()
        if (code == 1):
            break
        em()

        # Stop if finished
        if (bonds == bonds_tot):
            break

        # Molecular dynamics
        if (bonds % bonds_cyc == 0):
            setup_md(0)
            if ((bonds/bonds_cyc) % md_cyc == 0):
                md(2)
            else:
                md(1)

    # Finalization
    polym_final()

def polym_step():

    # Variables
    global bonds
    attempts = 1

    # Attempt until successful or max attempts
    while (1):

        # Polymerization step
        cmd = 'perl ../scripts/%s -i init.lmps -t ../types.txt -s ../scripts/%s -o data.lmps' \
            % (script_step, input_polym)
        sys.stdout.flush()
        code = subprocess.call(cmd.split())

        # Bond made
        if (code == 0):
            print ("  Attempts: %d" % attempts)
            if (not os.path.isfile('data.lmps')):
                err_exit("Polymerization step script did not complete properly.")
            return 0

        # No pair found
        elif (code == 3):

            # Stop if maximum attempts reached
            if (attempts > md_max):
                print ("  No pair was found within the maximum number of attempts.")
                bonds -= 1
                os.chdir('..')
                if (keep == 0):
                    for path in glob.glob('step_*'):
                        shutil.rmtree(path)
                return 1

            # Molecular dynamics
            setup_md(attempts)
            md(0)
            attempts += 1

        # Error
        else:
            err_exit("Polymerization step script did not complete properly.")

def polym_init():
    print( "Initialization:")
    if (script_init != 0):
        cmd = 'perl scripts/%s -i data.lmps -t types.txt -s scripts/%s -o temp.lmps' \
            % (script_init, input_polym)
        sys.stdout.flush()
        code = subprocess.call(cmd.split())
    else:
        shutil.copy('data.lmps', 'temp.lmps')
    if (code != 0 or not os.path.isfile('temp.lmps')):
        err_exit("Polymerization initialization script did not complete properly.")

def polym_final():
    print ("Finalization:")
    if (script_final != 0):
        cmd = 'perl scripts/%s -i temp.lmps -t types.txt -s scripts/%s -o final.lmps' \
            % (script_final, input_polym)
        sys.stdout.flush()
        code = subprocess.call(cmd.split())
    else:
        shutil.copy('temp.lmps', 'final.lmps')
    if (code != 0 or not os.path.isfile('final.lmps')):
        err_exit("Polymerization finalization script did not complete properly.")
    os.remove('temp.lmps')

def em():

    # Input script
    shutil.copy('../scripts/' + input_min, input_min)

    # LAMMPS EM
    cmd = '%s -i %s' % (lmps, input_min)
    code = subprocess.call(cmd.split(), stdout=open('out', 'w'))
    if (not os.path.isfile('min.lmps')):
        err_exit("LAMMPS energy minimization did not complete properly.")

    # Data file
    shutil.copy('min.lmps', '../temp.lmps')
    os.chdir('..')

    # Keep files?
    if (keep == 0):
        for path in glob.glob('step_*'):
            shutil.rmtree(path)

def md(num):

    # Input script
    if (num == 0):
        inp = input_md0
        shutil.copy('../../scripts/' + inp, inp)
    elif (num == 1):
        inp = input_md1
        shutil.copy('../scripts/' + inp, inp)
    elif (num == 2):
        inp = input_md2
        shutil.copy('../scripts/' + inp, inp)

    # LAMMPS MD
    cmd = '%s -i %s' % (lmps, inp)
    code = subprocess.call(cmd.split(), stdout=open('out', 'w'))
    if (not os.path.isfile('md.lmps')):
        err_exit("LAMMPS molecular dynamics did not complete properly.")

    # Data file
    if (num == 0):
        shutil.copy('md.lmps', '../init.lmps')
    else:
        shutil.copy('md.lmps', '../temp.lmps')
    os.chdir('..')

    # Keep files?
    if (keep == 0):
        if (num == 0):
            for path in glob.glob('md_*'):
                shutil.rmtree(path)
        else:
            for path in glob.glob('step_*'):
                shutil.rmtree(path)

def setup_step():

    # Directory
    print ("Step %d:" % bonds)
    directory = 'step_' + '{0:03d}'.format(bonds)
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        err_exit("Directory '%s' already exists." % (directory))
    os.chdir(directory)

    # Data file
    if (bonds == 0):
        shutil.copy('../temp.lmps', 'data.lmps')
    else:
        shutil.copy('../temp.lmps', 'init.lmps')

def setup_md(num):

    # Directory
    if (num == 0):
        directory = 'step_' + '{0:03d}'.format(bonds) + '_md'
    else:
        directory = 'md_' + '{0:03d}'.format(num)
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        err_exit("Directory '%s' already exists." % (directory))
    os.chdir(directory)

    # Data file
    if (num == 0):
        shutil.copy('../temp.lmps', 'data.lmps')
    else:
        shutil.copy('../init.lmps', 'data.lmps')

def print_header():
    print ("Polymatic Simulated Polymerization\n")
    print("Parameters\n----------")
    print ("Initial bonds:             %d" % bonds)
    print ("Total bonds:               %d" % bonds_tot)
    print ("Bonds per cycle:           %d" % bonds_cyc)
    print ("Frequency of MD type 2:    %d" % md_cyc)
    print ("Maximum bond attempts:     %d\n" % md_max)
    print ("Polymerization Loop\n-------------------")

def print_footer():
    print ("\nSummary\n-------")
    print ("Bonds made:                %d" % bonds)
    print ("Completion percentage:     %.f%%" % (bonds*100.0/bonds_tot))

def err_exit(error):
    print ("Error: %s" % error)
    sys.exit(1)

if __name__ == '__main__':
    main()

