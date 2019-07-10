#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import re
import os
import shutil
import uuid
import subprocess

from psi4 import core
from psi4 import extras
from psi4.driver.p4util.exceptions import *

class XDMDispersion(object):
    """Class for tXDM dispersion calculations

    """

    def __init__(self, a1=-1.0, a2=0.0, vol=""):
        self.a1 = a1
        self.a2 = a2
        self.vol = vol.lower()
        if (not extras.addons("postg")):
            raise XDMError("Cannot find the postg executable for the XDM dispersion correction")

    def run_postg(self,wfn=None,derint=0):
        """Calls postg and returns energies and derivatives"""

        # grab xdm parameters
        (a1,a2,vol) = (self.a1,self.a2,self.vol)

        # Validate arguments
        if wfn is None:
            raise ValidationError("Call to run_postg without a wavefunction.")
        if (a1 < 0 or len(vol) == 0):
            raise ValidationError("""Call to run_postg with incorrect parameters a1 = %.4f, a2 = %.4f, vol = %s.""" % (a1,a2,vol))

        # check that the volume token is either one of the postg tokens or a float
        allowed_vols = ['blyp','b3lyp','bhandhlyp','bhandh','bhah','bhahlyp','camb3lyp','cam-b3lyp','pbe','pbe0','lcwpbe','lc-wpbe','pw86','pw86pbe','b971','b97-1','hf']
        if vol not in allowed_vols:
            try:
                float(vol)
            except:
                raise XDMError("""Call to run_postg with invalid volume token: %s""" % vol)

        # Find environment by merging PSIPATH and PATH environment variables
        lenv = {
            'PATH': ':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':') if x != '']) + \
            ':' + os.environ.get('PATH'),
            'LD_LIBRARY_PATH': os.environ.get('LD_LIBRARY_PATH')
        }
        #   Filter out None values as subprocess will fault on them
        lenv = {k: v for k, v in lenv.items() if v is not None}
 
        # Setup unique scratch directory and move in (assume we are in psi4)
        current_directory = os.getcwd()
        psioh = core.IOManager.shared_object()
        psio = core.IO.shared_object()
        os.chdir(psioh.get_default_path())
        postg_tmpdir = 'psi.' + str(os.getpid()) + '.' + psio.get_default_namespace() + '.postg.' + str(uuid.uuid4())[:8]

        if os.path.exists(postg_tmpdir) is False:
            os.mkdir(postg_tmpdir)
            os.chdir(postg_tmpdir)
   
        # Write molden file, with no virtual orbitals
        try:
            occa = wfn.occupation_a()
            occb = wfn.occupation_b()
        except AttributeError as err:
            os.chdir(current_directory)
            raise XDMError("The wavefunction passed to run_postg does not have occupation numbers.") from err
        moldenfile = './tmp.molden'
        mw = core.MoldenWriter(wfn)
        mw.write(moldenfile, wfn.Ca(), wfn.Cb(), wfn.epsilon_a(), wfn.epsilon_b(), occa, occb, False)

        # Call postg program
        command = ['postg', '%.4f' % a1, '%.4f' % a2, moldenfile, vol]
        try:
            child = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=lenv)
        except OSError as e:
            os.chdir(current_directory)
            raise XDMError('Error executing postg: %s' % e)
        out, err = child.communicate()
        if child.returncode != 0:
            os.chdir(current_directory)
            raise XDMError("""Error running postg. Please check the temporary files in directory: %s\n--- Some info about the error from postg follows, maybe ---\n %s""" % (postg_tmpdir, err.decode('utf-8')))

        # Parse output
        for line in out.splitlines():
            line = line.decode('utf-8')
            if re.match('dispersion energy',line):
                sline = line.split()
                exdm = float(sline[2])
                
        # Clean up files and remove scratch directory
        os.chdir('..')
        try:
            shutil.rmtree(postg_tmpdir)
        except OSError as err:
            raise OSError("""Unable to remove postg temporary directory: %s""" % postg_tmpdir) from err
        os.chdir(current_directory)

        return exdm

## """Module with functions that interface with postg."""
## import re
## import uuid
## import socket
## import subprocess
## 
## from .util import parse_dertype
## from .molecule import Molecule
## 
## def run_gcp(self, func=None, dertype=None, verbose=False):  # dashlvl=None, dashparam=None
##
##     # Parse grad output
##     if derint != 0:
##         derivfile = './gcp_gradient'
##         dfile = open(derivfile, 'r')
##         dashdderiv = []
##         for line in geom.splitlines():
##             lline = line.split()
##             if len(lline) != 4:
##                 continue
##             if lline[0] == 'Gh':
##                 dashdderiv.append([0.0, 0.0, 0.0])
##             else:
##                 dashdderiv.append([float(x.replace('D', 'E')) for x in dfile.readline().split()])
##         dfile.close()
## 
##         if len(dashdderiv) != self.natom():
##             raise ValidationError('Program gcp gradient file has %d atoms- %d expected.' % \
##                 (len(dashdderiv), self.natom()))
## 
##     # Prepare results for Psi4
##     if isP4regime and derint != 0:
##         core.set_variable('GCP CORRECTION ENERGY', dashd)
##         psi_dashdderiv = core.Matrix.from_list(dashdderiv)
## 
##     # Print program output to file if verbose
##     if not verbose and isP4regime:
##         verbose = True if core.get_option('SCF', 'PRINT') >= 3 else False
##     if verbose:
## 
##         text = '\n  ==> GCP Output <==\n'
##         text += out.decode('utf-8')
##         if derint != 0:
##             with open(derivfile, 'r') as handle:
##                 text += handle.read().replace('D', 'E')
##             text += '\n'
##         if isP4regime:
##             core.print_out(text)
##         else:
##             print(text)
## 
##     # return -D & d(-D)/dx
##     if derint == -1:
##         return dashd, dashdderiv
##     elif derint == 0:
##         return dashd
##     elif derint == 1:
##         return psi_dashdderiv
## 
## try:
##     # Attach method to libmints psi4.Molecule class
##     core.Molecule.run_gcp = run_gcp
## except (NameError, AttributeError):
##     # But don't worry if that doesn't work b/c
##     #   it'll get attached to qcdb.Molecule class
##     pass
## 
