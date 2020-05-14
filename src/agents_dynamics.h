/*  AgentsDynamics
    defines movement rules/differential equations of agents in
    SwarmDynamics(swarmdyn.h, swarmdyn.cpp)
    v0.1, 13.5.2020

    (C) 2020 Pascal Klamser

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef agents_dynamics_H
#define agents_dynamics_H
// OWN MODULES:
#include "common_defines.h"
#include "agents.h"
#include "agents_operation.h"
// #include "mathtools.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <math.h>
#include <random>       // std::default_random_engine

void ParticleBurstCoast(particle &a, params * ptrSP, gsl_rng *r);
template<class agent>
void Boundary(agent &a, double sizeL,  int BC);             // calculate boundary conditions
void MovePredator(predator &pred, std::vector<particle> &a, params *ptrSP, gsl_rng *r);
void CreatePredator(std::vector<particle> &a, params *ptrSP,
                         predator &pred, gsl_rng *r);
void CreateFishNet(std::vector<particle> &a, params *ptrSP,
                   std::vector<predator> &preds, gsl_rng *r);
void MoveFishNet(std::vector<predator> &preds, params *ptrSP);
std::vector<double> predictX(std::vector<double> & x,
                             std::vector<double> & v,
                             std::vector<double> & force,
                             double t_burst, double up_rate,
                             double friction, bool alongForce=false);
// computes force needed to set v=0 in time=t_burst with friction
double force2stop(double v, double t_burst, double friction);
void FishNetKill(std::vector<particle> &a, std::vector<predator> &preds,
                 params *ptrSP);
void PredKill(std::vector<particle> &a, predator &pred, params *ptrSP,
              gsl_rng *r);
#endif
