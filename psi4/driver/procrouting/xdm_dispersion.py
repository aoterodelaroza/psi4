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

from psi4 import core
from psi4 import extras
from psi4.driver.p4util.exceptions import ValidationError, XDMError

class XDMDispersion(object):
    """Class for tXDM dispersion calculations

    """

    def __init__(self, a1=-1.0, a2=0.0, vol=""):
        self.a1 = a1
        self.a2 = a2
        self.vol = vol
        if (not extras.addons("postg")):
            raise XDMError("Cannot find the postg executable for the XDM dispersion correction")

    def run_postg(self,wfn=None,derint=0):
        """Calls postg and returns energies and derivatives"""

        # Validate arguments
        if wfn is None:
            raise ValidationError("Call to run_postg without a wavefunction.")
        if (self.a1 < 0 or len(self.vol) == 0):
            raise ValidationError("""Call to run_postg with incorrect parameters a1 = %.4f, a2 = %.4f, vol = %s.""" % (self.a1,self.a2,self.vol))

        filename = "blah.molden"
        mw = core.MoldenWriter(wfn)
        mw.write(filename, wfn.Ca(), wfn.Cb(), wfn.epsilon_a(), wfn.epsilon_b(), wfn.occupation_a(), wfn.occupation_b(), False)

        exdm = 0.0
        return exdm

## """Module with functions that interface with postg."""
## import os
## import re
## import uuid
## import shutil
## import socket
## import subprocess
## 
## try:
##     from psi4.driver.p4util.exceptions import *
##     from psi4 import core
##     isP4regime = True
## except ImportError:
##     from .exceptions import *
##     isP4regime = False
## from .util import parse_dertype
## from .molecule import Molecule
## 
## 
## def run_gcp(self, func=None, dertype=None, verbose=False):  # dashlvl=None, dashparam=None
##
##     # TODO temp until figure out paramfile
##     allowed_funcs = ['HF/MINIS', 'DFT/MINIS', 'HF/MINIX', 'DFT/MINIX',
##         'HF/SV', 'DFT/SV', 'HF/def2-SV(P)', 'DFT/def2-SV(P)', 'HF/def2-SVP',
##         'DFT/def2-SVP', 'HF/DZP', 'DFT/DZP', 'HF/def-TZVP', 'DFT/def-TZVP',
##         'HF/def2-TZVP', 'DFT/def2-TZVP', 'HF/631Gd', 'DFT/631Gd',
##         'HF/def2-TZVP', 'DFT/def2-TZVP', 'HF/cc-pVDZ', 'DFT/cc-pVDZ',
##         'HF/aug-cc-pVDZ', 'DFT/aug-cc-pVDZ', 'DFT/SV(P/h,c)', 'DFT/LANL',
##         'DFT/pobTZVP', 'TPSS/def2-SVP', 'PW6B95/def2-SVP',
##         # specials
##         'hf3c', 'pbeh3c']
##     allowed_funcs = [f.lower() for f in allowed_funcs]
##     if func.lower() not in allowed_funcs:
##         raise Dftd3Error("""bad gCP func: %s. need one of: %r""" % (func, allowed_funcs))
## 
##     # Move ~/.dftd3par.<hostname> out of the way so it won't interfere
##     defaultfile = os.path.expanduser('~') + '/.dftd3par.' + socket.gethostname()
##     defmoved = False
##     if os.path.isfile(defaultfile):
##         os.rename(defaultfile, defaultfile + '_hide')
##         defmoved = True
## 
##     # Find environment by merging PSIPATH and PATH environment variables
##     lenv = {
##         'PATH': ':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':') if x != '']) + \
##                 ':' + os.environ.get('PATH'),
##         'LD_LIBRARY_PATH': os.environ.get('LD_LIBRARY_PATH')
##         }
##     #   Filter out None values as subprocess will fault on them
##     lenv = {k: v for k, v in lenv.items() if v is not None}
## 
##     # Find out if running from Psi4 for scratch details and such
##     try:
##         import psi4
##     except ImportError as err:
##         isP4regime = False
##     else:
##         isP4regime = True
## 
##     # Setup unique scratch directory and move in
##     current_directory = os.getcwd()
##     if isP4regime:
##         psioh = core.IOManager.shared_object()
##         psio = core.IO.shared_object()
##         os.chdir(psioh.get_default_path())
##         gcp_tmpdir = 'psi.' + str(os.getpid()) + '.' + psio.get_default_namespace() + \
##             '.gcp.' + str(uuid.uuid4())[:8]
##     else:
##         gcp_tmpdir = os.path.expanduser('~') + os.sep + 'gcp_' + str(uuid.uuid4())[:8]
##     if os.path.exists(gcp_tmpdir) is False:
##         os.mkdir(gcp_tmpdir)
##     os.chdir(gcp_tmpdir)
## 
##     # Write gcp_parameters file that governs cp correction
## #    paramcontents = gcp_server(func, dashlvl, 'dftd3')
## #    paramfile1 = 'dftd3_parameters'  # older patched name
## #    with open(paramfile1, 'w') as handle:
## #        handle.write(paramcontents)
## #    paramfile2 = '.gcppar'
## #    with open(paramfile2, 'w') as handle:
## #        handle.write(paramcontents)
## 
## ###Two kinds of parameter files can be read in: A short and an extended version. Both are read from
## ###$HOME/.gcppar.$HOSTNAME by default. If the option -local is specified the file is read in from
## ###the current working directory: .gcppar
## ###The short version reads in: basis-keywo
## 
##     # Write dftd3_geometry file that supplies geometry to dispersion calc
##     numAtoms = self.natom()
##     geom = self.save_string_xyz()
##     reals = []
##     for line in geom.splitlines():
##         lline = line.split()
##         if len(lline) != 4:
##             continue
##         if lline[0] == 'Gh':
##             numAtoms -= 1
##         else:
##             reals.append(line)
## 
##     geomtext = str(numAtoms) + '\n\n'
##     for line in reals:
##         geomtext += line.strip() + '\n'
##     geomfile = './gcp_geometry.xyz'
##     with open(geomfile, 'w') as handle:
##         handle.write(geomtext)
##     # TODO somehow the variations on save_string_xyz and
##     #   whether natom and chgmult does or doesn't get written
##     #   have gotten all tangled. I fear this doesn't work
##     #   the same btwn libmints and qcdb or for ghosts
## 
##     # Call gcp program
##     command = ['gcp', geomfile]
##     command.extend(['-level', func])
##     if derint != 0:
##         command.append('-grad')
##     try:
##         #print('command', command)
##         dashout = subprocess.Popen(command, stdout=subprocess.PIPE, env=lenv)
##     except OSError as e:
##         raise ValidationError('Program gcp not found in path. %s' % e)
##     out, err = dashout.communicate()
## 
##     # Parse output
##     success = False
##     for line in out.splitlines():
##         line = line.decode('utf-8')
##         if re.match('  Egcp:', line):
##             sline = line.split()
##             dashd = float(sline[1])
##         if re.match('     normal termination of gCP', line):
##             success = True
## 
##     if not success:
##         os.chdir(current_directory)
##         raise Dftd3Error("""Unsuccessful gCP run.""")
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
## #    # Clean up files and remove scratch directory
## #    os.unlink(paramfile1)
## #    os.unlink(paramfile2)
## #    os.unlink(geomfile)
## #    if derint != 0:
## #        os.unlink(derivfile)
## #    if defmoved is True:
## #        os.rename(defaultfile + '_hide', defaultfile)
## 
##     # clean up files and remove scratch directory
##     os.chdir('..')
##     try:
##         shutil.rmtree(gcp_tmpdir)
##     except OSError as err:
##         raise OSError('Unable to remove gcp temporary directory: {}'.format(gcp_tmpdir)) from err
##     os.chdir(current_directory)
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
