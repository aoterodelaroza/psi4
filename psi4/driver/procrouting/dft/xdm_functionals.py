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
"""
List of XDM-corrected functionals
"""

funcs = []

funcs.append({
    "name": "BLYP-XDM",
    "x_functionals": {
        "GGA_X_B88": {}
    },
    "c_functionals": {
        "GGA_C_LYP": {}
    },
    "citation":
    '    A. D. Becke, Phys. Rev. A 38, 3098 (1988).\n' + \
    '    C. Lee, W. Yang, and R. G. Parr, Phys. Rev. B 37, 785 (1988).\n' + \
    '    A. Otero-de-la Roza, E. R. Johnson, J. Chem. Phys. 138, 204109 (2013).\n',
    "description":
    '    BLYP GGA Exchange-Correlation Functional plus XDM dispersion.\n',
    "dispersion": {
        "type": "xdm",
        "params": {
            's6': 1.25,
            'alpha6': 20.0,
            'sr6': 1.1
        }
    }
})

functional_list = {}
for functional in funcs:
    functional_list[functional["name"].lower()] = functional

