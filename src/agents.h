/*  agents
    defines structures/classes of agents used in SwarmDynamics(swarmdyn.h, swarmdyn.cpp)
    v0.1, 13.5.2020

    (C) 2020 Pascal Klamser, Pawel Romanczuk

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

#ifndef agents_H
#define agents_H
#include "common_defines.h"
#include <gsl/gsl_rng.h>
#include <H5Cpp.h>      // for hdf5 output
#include <set>
#include <vector>

struct particle{
    std::vector<double> x;            // agent position vector
    std::vector<double> v;            // agent velocity vector
    std::vector<double> vold;         // agent velocity vector from step before
    std::vector<double> u;            // agent direction unit vector
    double phi;             // agent direction polar angle [0,2Pi]
    double vproj;           // vel along heading direction
    std::vector<int>    cell;         // cell index for spatial sub-division (linked list algorithm)
    std::vector<double> force;        // total social force vector
    std::vector<double> force_att;    // attraction vector
    std::vector<double> force_rep;    // repulsion vector
    std::vector<double> force_sor;    // social repulsion vector
    std::vector<double> force_alg;    // alignment vector
    std::vector<double> force_flee;   // flee vector
    double fitness;         // current fitness of particle
    double fit_decrease;    // how much fitness got decreased in current time_step
    bool dead;              // if killed by predator
    unsigned int bin_step;
    unsigned int steps_till_burst;
    unsigned int id;        // ID of each agent (needed to split and merge agents correctly)
    std::vector<unsigned int> NN;    // vector containing all NN

    // counters of interaction partners - important only for local metric coupling (not global)
    int counter_rep;        // counter repulsion partners 
    int counter_sor;        // counter social repulsion partners 
    int counter_alg;        // counter allignment partners 
    int counter_att;        // counter attraction partners 
    int counter_flee;       // counter flee
    std::vector<double> out(void);  // for output
};
typedef struct particle particle;


struct predator{
    unsigned int id;        // ID of each agent (needed to split and merge agents correctly)
    std::vector<double> x;            // agent position vector
    std::vector<double> v;            // agent velocity vector
    std::vector<double> u;            // agent direction unit vector
    double phi;             // agent direction polar angle
    double phi_start;       // phi at creation of predator
    double vproj;           // vel along heading direction
    std::vector<int>    cell;         // cell index for spatial sub-division (linked list algorithm)
    std::vector<double> force;        // predation force vector
    unsigned int kills; // counts how many prey got killed (only if pred_kill == 1)
    unsigned int state;     // state of predator:0=approach, 1=hunt
    std::vector<unsigned int> NN;    // vector containing all F seen by P
    std::set<unsigned int> NNset;    // set of all F seeing P
    std::set<unsigned int> NN2set;   // set of all F seeing NNset which are not in NNset
    std::set<unsigned int> NN3set;   // set of all F seeing NN2set  which are not in NNset or NN2set
    std::vector<double> out(void);
};
typedef struct predator predator;

// data structure for system parameters, and auxiliary variables
struct params{

    std::string location;   // path of output-files
    std::string fileID;     // name of output-files
    unsigned int N;                  // number of agents
    int Npred;              // 1 if predator, 0 if not
    int Ndead;              // # of dead/killed agents
    double sizeL;           // system size
    double Dphi;            // angular noise intensity  (= variance/2)
    double noisep;          // auxiliary var. - std. deviation of the angular noise

    double soc_strength;    // social strength
    double env_strength;    // environmental strength
    double burst_rate;    // alignment strength

    double prob_social;    // social rep. strength

    double rep_range;       // repulsion range
    double alg_range;       // alignment range
    double att_range;       // attraction range

    double alphaTurn;
    double flee_range;      // flee range

    double pred_time;       // time at which the predator comes in
    double pred_speed0;     // how fast the predator moves

    double sim_time;        // simulation time in natural time units
    double dt;              // time step of integration
    double output;          // output time in natural units

    int BC;                 // switch for boundary conditions
    int IC;                 // switch for initial conditions
    int N_confu;// determines if the predator 1: strictly moves relative to COM, 2: follows the COM, 3: adjust start position to hit com with straight move
    int pred_kill;         // if predator kills prey or flies through
    int pred_move;         // if predator moves randomly (=0) OR follows closest (=1)

    double trans_time;      // transient time before output starts

    int sim_steps;          // auxiliary var. - total number of integration steps int(sim_time/dt)
    int step_output;        // auxiliary var. - integrations steps between outputs int(output/dt)

    int output_mode;        // switch for output data (full,mean)
    bool out_extend;        // derived from output_mode
    bool out_mean;          // derived from output_mode
    bool out_particle;      // derived from output_mode
    int out_h5;             // switch for ouput data format (txt, HDF5)
    unsigned int outstep;       // current output step (needed for hdf5)
    unsigned int outstep_pred;  // current output step (needed for hdf5)
    unsigned int total_outstep; // total Nr of output steps
    unsigned int total_outstep_pred; // total Nr of output steps

    double kill_rate;        // kills per second if prob_selected = 1 and prob_catch = 1
    double kill_range;       // distance between pred and prey at which prob_kill > 0

    // output-arrays
    std::vector< std::vector<double> > dataOutMean;
    std::vector< std::vector<double> > dataOutSwarm;
    std::vector< std::vector<double> > dataOutSwarmPred;
    
    double burst_duration;
    unsigned int burst_steps;      // steps which prey stays in bin_mode
    double beta;            // relaxation rate of velocity along heading

    gsl_rng *r;
};
typedef struct params params;

#endif
