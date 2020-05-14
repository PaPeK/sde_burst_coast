/*  SocialForces
    defines social force models which describe the interactions
    between agents in SwarmDynamics(swarmdyn.h, swarmdyn.cpp)
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
#include "social_forces.h"

// Social-Force-Models------------------------------------------------
void SFM_4Zone(params *ptrSP, double u_ji[2],
               double dist_interaction, double v_ji[2],
               unsigned int c[3], double f0[2],
               double f1[2], double f2[2]){
    /* computes the 4 forces assuming input of c = {0, 0, 0, 0}
     * and fx = {0, 0}
     */
    // FORCE0 (repulsion)
    if(dist_interaction <= ptrSP->rep_range){
        c[0]++;
        f0[0] = u_ji[0];
        f0[1] = u_ji[1];
    }
    // FORCE1 (alignment)
    else if((dist_interaction > ptrSP->rep_range) && 
            (dist_interaction <= ptrSP->alg_range)){
        c[1]++;
        f1[0] = v_ji[0];
        f1[1] = v_ji[1];
    }
    // FORCE2 (attraction)
    else if((dist_interaction > ptrSP->alg_range) &&
            (dist_interaction <= ptrSP->att_range)){
        c[2]++;
        f2[0] = u_ji[0];
        f2[1] = u_ji[1];
    }
}
