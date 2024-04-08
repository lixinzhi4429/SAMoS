/* ***************************************************************************
 *
 *  Copyright (C) 2013-2016 University of Dundee
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
 * \file integrator_pairdiss_pos.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 29-Dec-2015
 * \author Silke Henkes shenkes@lorentz.leidenuniv.nl
 * \date 8 April 2024
 * \brief Declaration of IntegratorPairdissPos class
 */ 

#ifndef __INTEGRATOR_PAIRDISS_POS_H__
#define __INTEGRATOR_PAIRDISS_POS_H__

#include <cmath>

#include "integrator.hpp"
#include "rng.hpp"
#include <Eigen/Sparse>
#include "value.hpp"
//#include "neighbour_list.hpp"

using std::sqrt;

/*! IntegratorBrownian class implements Brownian dynamics for the particle position.
    Particle director will not be integrated. 
*/
class IntegratorPairdissPos : public Integrator
{
public:
  
  //! Constructor
  //! \param sys Pointer to a System object containing all particles
  //! \param msg Internal message handler
  //! \param pot Pairwise and external interaction handler
  //! \param align Pairwise and external alignment handler
  //! \param nlist Neighbour list object
  //! \param cons Enforces constraints to the manifold surface
  //! \param temp Temperature control object
  //! \param param Contains information about all parameters 
  IntegratorPairdissPos(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align, NeighbourListPtr nlist,  ConstrainerPtr cons, ValuePtr temp, pairs_type& param) : Integrator(sys, msg, pot, align, nlist, cons, temp, param)
  { 
    m_known_params.push_back("zeta");
    m_known_params.push_back("xipair");
    m_known_params.push_back("seed");
    m_known_params.push_back("r_int");
    m_known_params.push_back("use_particle_radii");
    m_known_params.push_back("phase_in");
    m_known_params.push_back("tphase");
    string param_test = this->params_ok(param);
    if (param_test != "")
    {
      m_msg->msg(Messenger::ERROR,"Parameter \""+param_test+"\" is not a valid parameter for brownian_pos integrator.");
      throw runtime_error("Unknown parameter \""+param_test+"\" in pairdiss_pos integrator.");
    }
    if (param.find("zeta") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Pairdiss dynamics integrator for particle position. Single particle friction not set. Using default value 1.");
      m_zeta = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Pairdiss dynamics integrator for particle position. Setting single particle friction to "+param["zeta"]+".");
      m_zeta = lexical_cast<double>(param["zeta"]);
    }
    if (param.find("xipair") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Pairdiss dynamics integrator for particle position. Pair particle friction not set. Using default value 1.");
      m_xipair = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Pairdiss dynamics integrator for particle position. Setting pair particle friction to "+param["xipair"]+".");
      m_xipair = lexical_cast<double>(param["xipair"]);
    }
    m_msg->write_config("integrator.pairdiss_pos.xipair",lexical_cast<string>(m_xipair));
    if (param.find("r_int") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No potential range (r_int) specified for Pairdiss dynamics integrator. Setting it to 1.3.");
      m_r_int = 1.3;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global interaction range (r_int) for Pairdiss dynamics integrator is set to "+param["r_int"]+".");
      m_r_int = lexical_cast<double>(param["r_int"]);
    }
    if (param.find("seed") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Pairdiss dynamics integrator for particle position. No random number generator seed specified. Using default 0.");
      m_rng = make_shared<RNG>(0);
      m_msg->write_config("integrator.pairdiss_pos.seed",lexical_cast<string>(0));
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Pairdiss dynamics integrator for particle position. Setting random number generator seed to "+param["seed"]+".");
      m_rng = make_shared<RNG>(lexical_cast<int>(param["seed"]));
      m_msg->write_config("integrator.pairdiss_pos.seed",param["seed"]);
    }
    if (param.find("use_particle_radii") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Pairdiss dynamics integrator is set to use particle radii to control its range. Parameter r_int will be multiplicative.");
      m_use_particle_radii = true;
      m_msg->write_config("integrator.pairdiss_pos.use_radii","true");
    }
    if (param.find("phase_in") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Pairdiss dynamics integrator. Gradually phasing in the pair dissipation for new particles.");
      m_phase_in = true;
      m_msg->write_config("integrator.pairdiss_pos.phase_in","true");
    }  
    if (param.find("tphase") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No phase in time specified for Pairdiss dynamics integrator. Setting it to 1000*dt.");
      m_tphase = 1000*m_dt;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Phase in time Pairdiss dynamics integrator is set to "+param["tphase"]+".");
      m_tphase = lexical_cast<double>(param["tphase"]);
    } 
    m_msg->msg(Messenger::WARNING,"Pairdiss dynamics integrator. FDT compliant fluctuations are not yet implemented. Defaults to single particle with fluctuation scaling mu = 1/(zeta+0.5*<z>*xipair). Do not use as thermal integrator.");
    
    
  }
  
  
  //! Propagate system for a time step
  void integrate();
  
  
private:
  
  RNGPtr  m_rng;          //!< Random number generator 
  double  m_zeta;           //!< single particle friction
  double  m_xipair;           //!< pair particle friction

  double m_r_int;           //!< interaction range (multiplicative of radii)
  bool m_use_particle_radii; //!< Use particle radii to determine contacts
  bool m_phase_in;          //!< Phase in pair part of dissipation for consistency
  double m_tphase;          //!< Phase in time scale


  // sparse matrix solver from Eigen
  Eigen::SparseLU<Eigen::SparseMatrix<double,Eigen::RowMajor>> sol;
  
  
};

typedef shared_ptr<IntegratorPairdissPos> IntegratorPairdissPosPtr;

#endif
