/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/**********************************************************
 * xdm_dispersion.h: XDM dispersion correction.
 * Alberto Otero de la Roza <aoterodelaroza@gmail.com> and Joseph Weatherby <jweatherby@upei.ca>
 * June 27th, 2019
 * Based on the dispersion code by Robert Parrish <robparrish@gmail.com>
 ***********************************************************/

#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "xdm_dispersion.h"

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <memory>

namespace psi {

  XDMDispersion::XDMDispersion() {}

  XDMDispersion::~XDMDispersion() {}

  std::shared_ptr<XDMDispersion> XDMDispersion::build(double a1, double a2) {
    auto disp = std::make_shared<XDMDispersion>();
    disp->a1_ = a1;
    disp->a2_ = a2;
    return disp;
  }

  void XDMDispersion::print(std::string out, int level) const {
    if (level < 1) return;
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    printer->Printf("   => XDM Dispersion <=\n\n");

    printer->Printf("XDM dispersion");
    printer->Printf("\n");

    printer->Printf("XDM citation");
    printer->Printf("\n");

    printer->Printf("    a1 = %12.4f   a2 = %12.4f   \n",a1_,a2_);
    printer->Printf("\n");
  }

  std::string XDMDispersion::print_energy(std::shared_ptr<Molecule> m) {
    double e = compute_energy(m);
    std::stringstream s;
    s.setf(std::ios::scientific);
    s.precision(11);

    s << "   " << " XDMDispersion Energy: " << e << " [Eh]" << std::endl;

    return s.str();
  }

  std::string XDMDispersion::print_gradient(std::shared_ptr<Molecule> m) {
    SharedMatrix G = compute_gradient(m);
    double *g = G->pointer()[0];
    std::stringstream s;
    s.setf(std::ios::scientific);
    s.precision(11);

    s << "   " << " XDMDispersion Gradient ([a.u.]): " << std::endl << std::endl;
    s << "    Atom #:       E_x                E_y                 E_z" << std::endl;
    s << "   -----------------------------------------------------------------" << std::endl;

    for (int k = 1; k <= m->natom(); k++) {
      // clang-format off
      s << "  " << std::setw(5) << k <<
        std::setw(20) << g[(k - 1) * 3 + 0] <<
        std::setw(20) << g[(k - 1) * 3 + 1] <<
        std::setw(20) << g[(k - 1) * 3 + 2] << std::endl;
      // clang-format on
    }
    return s.str();
  }

  std::string XDMDispersion::print_hessian(std::shared_ptr<Molecule> m) {
    SharedMatrix H = compute_hessian(m);
    double **h = H->pointer();

    std::stringstream s;
    s.setf(std::ios::scientific);
    s.precision(11);

    s << "   " << " XDMDispersion Hessian ([a.u.]): " << std::endl << std::endl;
    for (int k = 1; k <= m->natom(); k++) {
        for (int j = 1; j <= m->natom(); j++) {
            // clang-format off
            s << "    Atom Pair A = " << k << " B = " << j << ":" << std::endl << std::endl;
            s << "                   xB                 yB                  zB" << std::endl;
            s << "   -----------------------------------------------------------------" << std::endl;

            s << "  " << std::setw(5) << "xA" <<
              std::setw(20) << h[(k - 1) * 3 + 0][(j - 1) * 3 + 0] <<
              std::setw(20) << h[(k - 1) * 3 + 0][(j - 1) * 3 + 1] <<
              std::setw(20) << h[(k - 1) * 3 + 0][(j - 1) * 3 + 2] << std::endl;
            s << "  " << std::setw(5) << "yA" <<
              std::setw(20) << h[(k - 1) * 3 + 1][(j - 1) * 3 + 0] <<
              std::setw(20) << h[(k - 1) * 3 + 1][(j - 1) * 3 + 1] <<
              std::setw(20) << h[(k - 1) * 3 + 1][(j - 1) * 3 + 2] << std::endl;
            s << "  " << std::setw(5) << "zA" <<
              std::setw(20) << h[(k - 1) * 3 + 2][(j - 1) * 3 + 0] <<
              std::setw(20) << h[(k - 1) * 3 + 2][(j - 1) * 3 + 1] <<
              std::setw(20) << h[(k - 1) * 3 + 2][(j - 1) * 3 + 2] << std::endl;
            s << std::endl;
            // clang-format on
        }
    }
    return s.str();
  }

  double XDMDispersion::compute_energy(std::shared_ptr<Molecule> m) {
    double E = -1.0;

    printf("xxxx XDM's compute_energy has been called! xxxx\n");
    printf("xxxx with a1 = %.4f\n",a1_);
    printf("xxxx with a2 = %.4f\n",a2_);
    printf("xxxx for now, returning an energy of 1 Hartree\n");

    return E;
  }

  SharedMatrix XDMDispersion::compute_gradient(std::shared_ptr<Molecule> m) {
    auto G = std::make_shared<Matrix>("XDMDispersion Gradient", m->natom(), 3);
    double **Gp = G->pointer();

    for (int i = 0; i < m->natom(); i++) {
        for (int j = 0; j < 3; j++) {
            Gp[i][j] = 0.0;
        }
    }

    return G;
  }

  SharedMatrix XDMDispersion::compute_hessian(std::shared_ptr<Molecule> m) {
    throw PSIEXCEPTION("XDMDispersion: Hessians not implemented");
  }

}  // end namespace
