/* ***************************************************************************
 *
 *  Copyright (C) 2013-2022 University of Dundee
 *  All rights reserved.
 *
 *  This file is part of SAMoS (Soft Active Matter on Surfaces) program.
 *
 *  SAMoS is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  SAMoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ****************************************************************************/

/*!
 * \file pair_bare_coulomb_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-Jun-2022
 * \brief Implementation of PairBareCoulombPotential class
 */

#include "pair_bare_coulomb_potential.hpp"

void PairBareCoulombPotential::compute(double dt)
{
    int N = m_system->size();
    double alpha = m_alpha;

    if (m_system->compute_per_particle_energy())
    {
        for (int i = 0; i < N; i++)
        {
            Particle &p = m_system->get_particle(i);
            p.set_pot_energy("barecoulomb", 0.0);
        }
    }
    vector<double> fx(N, 0.0);
    vector<double> fy(N, 0.0);
    vector<double> fz(N, 0.0);
    vector<vector<double>> fx_loc;
    vector<vector<double>> fy_loc;
    vector<vector<double>> fz_loc;

    int tid, numth;

#pragma omp parallel 
    {
        #pragma omp master
        numth = omp_get_num_threads();
    }
    fx_loc.resize(numth);
    fy_loc.resize(numth);
    fz_loc.resize(numth);
    for (int i = 0; i < numth; i++)
    {
        fx_loc[i].resize(N);
        fy_loc[i].resize(N);
        fz_loc[i].resize(N);
    }
    m_potential_energy = 0.0;
#pragma omp parallel private(tid)
    {
        tid = omp_get_thread_num();
        for (int i = 0; i < N; i++)
        {
            fx_loc[tid][i] = 0;
            fy_loc[tid][i] = 0;
            fz_loc[tid][i] = 0;
        }
#pragma omp for schedule(auto) nowait reduction(+ \
                                                : m_potential_energy)
        for (int i = 0; i < N; i++)
        {
            Particle &pi = m_system->get_particle(i);
            for (int j = i + 1; j < N; j++)
            {
                Particle &pj = m_system->get_particle(j);
                if (m_has_pair_params)
                {
                    int pi_t = pi.get_type() - 1, pj_t = pj.get_type() - 1;
                    alpha = m_pair_params[pi_t][pj_t].alpha;
                }
                double dx = pj.x - pi.x, dy = pj.y - pi.y, dz = pj.z - pi.z;
                // m_system->apply_periodic(dx, dy, dz);
                double r_sq = dx * dx + dy * dy + dz * dz;
                double r = sqrt(r_sq);
                // Handle potential
                double potential_energy = alpha / r;
                m_potential_energy += potential_energy;
                // Handle force
                double r_3 = r * r_sq;
                double force_factor = alpha / r_3;
                fx_loc[tid][i] -= force_factor * dx;
                fy_loc[tid][i] -= force_factor * dy;
                fz_loc[tid][i] -= force_factor * dz;

                // Use 3d Newton's law
                fx_loc[tid][j] += force_factor * dx;
                fy_loc[tid][j] += force_factor * dy;
                fz_loc[tid][j] += force_factor * dz;
            }
        }
#pragma omp critical
        for (int i = 0; i < N; i++)
        {
            fx[i] += fx_loc[tid][i];
            fy[i] += fy_loc[tid][i];
            fz[i] += fz_loc[tid][i];
        }
    }
    for (int i = 0; i < N; i++)
    {
        Particle &pi = m_system->get_particle(i);
        pi.fx += fx[i];
        pi.fy += fy[i];
        pi.fz += fz[i];
    }
}
