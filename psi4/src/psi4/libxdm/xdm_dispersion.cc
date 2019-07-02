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
#include "psi4/libfock/cubature.h"
#include "psi4/libfock/points.h"
#include "psi4/libfock/v.h"
#include "xdm_dispersion.h"

#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {

  XDMDispersion::XDMDispersion() {}

  XDMDispersion::~XDMDispersion() {}

  std::shared_ptr<XDMDispersion> XDMDispersion::build(double a1, double a2) {
    auto disp = std::make_shared<XDMDispersion>();
    outfile->Printf("\nxxxx in routine build\n");
    disp->a1_ = a1;
    disp->a2_ = a2;
    return disp;
  }

  void XDMDispersion::print(std::string out, int level) const {
    if (level < 1) return;
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    outfile->Printf("\nxxxx in routine print\n");
    printer->Printf("   => XDM Dispersion <=\n\n");

    printer->Printf("XDM dispersion");
    printer->Printf("\n");

    printer->Printf("XDM citation");
    printer->Printf("\n");

    printer->Printf("    a1 = %12.4f   a2 = %12.4f   \n",a1_,a2_);
    printer->Printf("\n");
  }

  std::string XDMDispersion::print_energy(std::shared_ptr<Molecule> m) {
    // double e = compute_energy(m);
    double e = 0.0;
    std::stringstream s;
    s.setf(std::ios::scientific);
    s.precision(11);

    s << "   " << " XDMDispersion Energy: " << e << " [Eh]" << std::endl;

    outfile->Printf("\nxxxx in routine print_energy\n");
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
    outfile->Printf("\nxxxx in routine print_gradient\n");
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
    outfile->Printf("\nxxxx in routine print_hessian\n");
    return s.str();
  }

  double XDMDispersion::compute_energy(std::shared_ptr<scf::HF> hf) {
    double E = -100000.0;

    std::map<std::string, std::string> opt_map;
    opt_map["DFT_PRUNING_SCHEME"] = "FLAT";

    std::map<std::string, int> opt_int_map;
    opt_int_map["DFT_RADIAL_POINTS"] = 100; // options_.get_int("DFT_VV10_RADIAL_POINTS");
    opt_int_map["DFT_SPHERICAL_POINTS"] = 302; // options_.get_int("DFT_VV10_SPHERICAL_POINTS");
    DFTGrid xdmgrid = DFTGrid(hf->molecule(), hf->V_potential()->basis(), opt_int_map, opt_map, hf->options());

    // hf->basis()
    // hf->molecule()
    // hf->V_potential()
    // hf->V_potential()->basis()
    // hf->options()
    // hf->functional()

    // hf->V_potential()->grid()
    //   std::shared_ptr<DFTGrid> grid_;
    // hf->V_potential()->functional()
    //   std::shared_ptr<SuperFunctional> functional_;
    // hf->V_potential()->properties()
    //   std::vector<std::shared_ptr<PointFunctions>> point_workers_;
    // std::vector<std::shared_ptr<SuperFunctional>> functional_workers_;

    // std::shared_ptr<BlockOPoints> get_block(int block);
    // size_t nblocks();
    // std::map<std::string, double>& quadrature_values() { return quad_values_; }
    // std::vector<std::shared_ptr<PointFunctions>> pw = hf->V_potential()->properties();
    // std::map<std::string, SharedVector>& pv = pw[0].point_values()

    // std::map<std::string, SharedVector>& point_values() { return point_values_; }
    // point_values_["RHO_A"] = std::make_shared<Vector>("RHO_A", max_points_);

    // see prepare_vv10_cache in v.cc for the following on  how to prepare the workers
    int rank = 0;

    const int max_points = xdmgrid.max_points();
    const int max_functions = xdmgrid.max_functions();
    std::vector<std::shared_ptr<PointFunctions>> xdm_point_workers;

    int num_threads_ = 1;
#ifdef _OPENMP
    num_threads_ = omp_get_max_threads();
#endif

    for (size_t i = 0; i < num_threads_; i++) {
      auto point_tmp = std::make_shared<RKSFunctions>(hf->V_potential()->basis(), max_points, max_functions);
      point_tmp->set_ansatz(2); // to get up to the laplacian
      point_tmp->set_pointers(hf->V_potential()->Dao()[0]);
      xdm_point_workers.push_back(point_tmp);
    }

    std::vector<std::map<std::string, SharedVector>> xdm_tmp_cache;
    xdm_tmp_cache.resize(xdmgrid.blocks().size());

    double rhosum = 0.0;
#pragma omp parallel for private(rank) schedule(guided) num_threads(num_threads_)
    for (size_t Q = 0; Q < xdmgrid.blocks().size(); Q++) {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        
        // Get workers and compute data
        // std::shared_ptr<SuperFunctional> fworker = functional_workers_[rank];
        std::shared_ptr<PointFunctions> pworker = xdm_point_workers[rank];
        std::shared_ptr<BlockOPoints> block = xdmgrid.blocks()[Q];

        pworker->compute_points(block,false);
        xdm_tmp_cache[Q] = pworker->point_values();

        int npoints = block->npoints();
        double* x = block->x();
        double* y = block->y();
        double* z = block->z();
        double* w = block->w();
        double* rho_a = pworker->point_value("RHO_A")->pointer();
        double* lap_a = pworker->point_value("LAPL_RHO_A")->pointer();
        double* grho2_a = pworker->point_value("GAMMA_AA")->pointer();
        double* tau_a = pworker->point_value("TAU_A")->pointer();

        // printf("I am in block %zu with %zu points\n",Q,npoints);
        // for (int i=0; i<npoints; i++){
        //   printf("%.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
        //          x[i],y[i],z[i],rho_a[i],lap_a[i],grho2_a[i],tau_a[i]);
        //   rhosum += w[i] * rho_a[i];
        // }
    }
    // printf("rho_sum = %.10f\n",rhosum);

    // printf("molecule\n");
    // for (int i=0; i<hf->molecule()->natom(); i++){
    //   printf("%s %.10f %.10f %.10f\n",hf->molecule()->symbol(i).c_str(),
    //          hf->molecule()->x(i),hf->molecule()->y(i),hf->molecule()->z(i));
    // }

    // printf("done!\n");
    // exit(0);

    // // how is the exc calculated? (v.cc)
    // for (size_t Q = 0; Q < xdmgrid.blocks().size(); Q++) {
    //   std::shared_ptr<BlockOPoints> block = xdmgrid.blocks()[Q];
    //   // std::shared_ptr<PointFunctions> pworker = point_workers_[rank];

    // //   // calculate rho, etc.
    // //   pworker->compute_points(block, false);
    // //   // do the quadrature
    // //   std::vector<double> qvals = dft_integrators::rks_quadrature_integrate(block, fworker, pworker);
    // //   functionalq[rank] += qvals[0];
    // //   dft_integrators::rks_integrator(block, fworker, pworker, V_local[rank]);
    // }
    // // quad_values_["FUNCTIONAL"] = std::accumulate(functionalq.begin(), functionalq.end(), 0.0);
    // printf("here\n");
    // exit(0);

    outfile->Printf("\nxxxx in routine compute_energy\n");
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

    outfile->Printf("\nxxxx in routine compute_gradient\n");
    return G;
  }

  SharedMatrix XDMDispersion::compute_hessian(std::shared_ptr<Molecule> m) {
    outfile->Printf("\nxxxx in routine compute_hessian\n");
    throw PSIEXCEPTION("XDMDispersion: Hessians not implemented");
  }

}  // end namespace
