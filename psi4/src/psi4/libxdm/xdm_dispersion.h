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

#ifndef xdm_dispersion_h
#define xdm_dispersion_h

/**********************************************************
 * xdm_dispersion.h: XDM dispersion correction.
 * Alberto Otero de la Roza <aoterodelaroza@gmail.com> and Joseph Weatherby <jweatherby@upei.ca>
 * June 27th, 2019
 * Based on the dispersion code by Robert Parrish <robparrish@gmail.com>
 ***********************************************************/
#include "psi4/psi4-dec.h"
#include <string>

namespace psi {

  class Molecule;

  class XDMDispersion {
  protected:
    double a1_;
    double a2_;

  public:
    XDMDispersion();
    virtual ~XDMDispersion();

    static std::shared_ptr<XDMDispersion> build(double a1 = 0.0, double a2 = 0.0);

    double get_a1() const { return a1_; }
    double get_a2() const { return a2_; }

    void set_a1(double a1) { a1_ = a1; }
    void set_a2(double a2) { a2_ = a2; }

    std::string print_energy(std::shared_ptr<Molecule> m);
    std::string print_gradient(std::shared_ptr<Molecule> m);
    std::string print_hessian(std::shared_ptr<Molecule> m);

    virtual double compute_energy(std::shared_ptr<Molecule> m);
    virtual SharedMatrix compute_gradient(std::shared_ptr<Molecule> m);
    virtual SharedMatrix compute_hessian(std::shared_ptr<Molecule> m);

    virtual void print(std::string out_fname = "outfile", int level = 1) const;
    void py_print() const { print("outfile", 1); }
  };
}

#endif
