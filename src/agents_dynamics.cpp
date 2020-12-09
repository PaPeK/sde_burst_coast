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
#include "agents_dynamics.h"

void draw_social_or_environmental_force(particle &a, params *ptrSP,
										gsl_rng *r, std::vector<double> &force,
										std::vector<double> &hvec, double &force_mag)
{
		double helper = 0;
		double dt = ptrSP->dt;
		double lphi = 0.0;
		double sizeL = ptrSP->sizeL;
		double beta = ptrSP->beta;
		double prob_social = ptrSP->prob_social;
		double draw_social = gsl_rng_uniform(r);
		
        if(draw_social <= prob_social)
        {
            force_mag = ptrSP->soc_strength;
            
            // implementation corresponds to couzin-model:
            if(a.counter_rep > 0)
            {
                vec_add221(force, a.force_rep);
            }
            else
            {
                // if both alg AND attraction: weighted average!
                if(a.force_alg[0] != 0 || a.force_alg[1] != 0)
                {
                    // hvec = vec_set_mag(a.force_alg, 1);
                    hvec = vec_set_mag(a.force_alg, a.counter_alg);
                    vec_add221(force, hvec);
                }
                if(a.force_att[0] != 0 || a.force_att[1] != 0)
                {
                    // hvec = vec_set_mag(a.force_att, 1);
                    hvec = vec_set_mag(a.force_att, a.counter_att);
                    vec_add221(force, hvec);
                }
            }
            
            // normalize:
            if(force[0] != 0 || force[1] != 0)
            {
                force = vec_set_mag(force, 1);
            }
            // if no social-force -> swim straight
            else
            {
                force[0] = cos(a.phi);
                force[1] = sin(a.phi);
            }
        }
        // if no force_flee -> draw random environmental cue
        else
        {
            if(a.counter_flee > 0)
            {
				force = a.force_flee;
                force = vec_set_mag(force, 1);
            }
            else
            {
                force_mag = ptrSP->env_strength;
                lphi = 2 * M_PI * gsl_rng_uniform(r);
                // random direction only forwards:
                // lphi = a.phi + (M_PI * gsl_rng_uniform(r) - M_PI / 2);
                force[0] = cos(lphi);
                force[1] = sin(lphi);
            }
        }
        
        // WALL: only if not periodic boundary conditions
        if(ptrSP->BC >= 5)
        {
            force = vec_set_mag(force, force_mag);
            a.force = force;
            
            //------------------ WALL collision avoidance ------------------
            // USE this if increased turning enabled (quick direction change)
            // estimate position at future burst
            // if position is outside of tank: adopt force such that agent is 1BL away from wall 
            // 1: Estimate future position x_fut  
            // 2: if x_fut outside of Tank use the force closest to the intended force 
            //      which ensures that the agent is inside the tank at the next burst
            std::vector<double> x_fut = predictXatNextBurst(a.x, a.v, a.force, 
                                                            dt * ptrSP->burst_steps,
                                                            dt * a.steps_till_burst,
                                                            beta);
            double r_fut = vec_length(x_fut);
            
            if(r_fut > sizeL - 2)
            {
                helper = closestForceDirection(a, ptrSP);
                force[0] = force_mag * cos(helper);
                force[1] = force_mag * sin(helper);
            }
        }
        
        force = vec_set_mag(force, force_mag);
        a.force = force;
}

bool overshoot_check(particle &a, std::vector<double> &force, double &force_mag, double &lphi) 
{ 
	std::vector<double> new_u {cos(lphi), sin(lphi)};
    
    double angForceV0 = acos(vec_dot(a.u, force) / force_mag);
    double angForceV1 = acos(vec_dot(new_u, force) / force_mag);
    double angV0V1 = acos(vec_dot(new_u, a.u));
    
    bool OvershootI = (angForceV0 < angForceV1);
    bool OvershootII = not OvershootI and (angV0V1 > angForceV0);
    bool OvershootIII = fabs(lphi - a.phi) > M_PI;
    bool overshoot=fabs(lphi - a.phi) > 0.01 and (OvershootI or OvershootII or OvershootIII); 
    
    return overshoot;
}

void consider_boundary(particle &a, params *ptrSP)
{
	int BC = ptrSP->BC;
    double sizeL = ptrSP->sizeL;
    
	if(BC!=-1)
    {
        Boundary(a, sizeL, BC);
    }
}

void ParticleBurstCoast(particle &a, params *ptrSP, gsl_rng *r)
{
	std::vector<double> force(2);
	std::vector<double> hvec(2);
	
    double force_mag = ptrSP->soc_strength;
    
    bool first_burst = (a.bin_step == ptrSP->burst_steps);
    bool bursting = (a.bin_step > 0);
   
    if(first_burst)
    {
		draw_social_or_environmental_force(a, ptrSP, r, force, hvec, force_mag);
    }
    else if(bursting)
    {
		// burst-mode: keep initial force
        force = a.force;
    }
    else
    {
		// coast-mode: no force
        a.force[0] = force[0] = 0;
        a.force[1] = force[1] = 0;
    }
    
    
    if(bursting)
    {
        a.bin_step -= 1;
    }
    if(a.steps_till_burst > 0)
    {
        a.steps_till_burst -= 1;
    }
    else
    {
        a.steps_till_burst = 0;
    }
    
    
    // Calculate polar angle
    double lphi = a.phi;
    double vproj = 0.0;
	vproj = a.vproj;     // to use correct time-step
    
    // speed adjustment
    double forcev = force[0] * cos(lphi) + force[1] * sin(lphi);
    
	double dt = ptrSP->dt;
	double beta = ptrSP->beta;
    // a.vproj += (-beta * vproj * vproj * vproj + forcev) * dt;
    a.vproj += (-beta * vproj + forcev) * dt;
    
    // a.vproj += rnv;
    // prevents F of swimming back
    if (a.vproj < 0)
    {
      a.vproj = 0.001;
      lphi += M_PI / 2;
    }
    
    // normal turn:
    double forcep= -force[0] * sin(lphi) + force[1] * cos(lphi);
    double alphaTurn = ptrSP->alphaTurn;
    lphi += alphaTurn * forcep * dt / vproj;
   
   
   // due to vproj-dependence extremely large values might occur
    if(overshoot_check(a, force, force_mag, lphi))
    {
		//set direction to force-direction
        lphi = atan2(force[1], force[0]);
    }
    
    
    lphi = fmod(lphi, 2*M_PI);
    a.phi = lphi;
    a.u[0] = cos(lphi);
    a.u[1] = sin(lphi);
    
    // Move particles with speed in units of [vel.al. range. / time]
    a.v[0] = a.vproj*a.u[0];
    a.v[1] = a.vproj*a.u[1];
    a.x[0] += a.v[0]*dt;
    a.x[1] += a.v[1]*dt;
    
    // Reset all forces
    a.force_rep[0] = a.force_rep[1] = 0.0;
    a.force_att[0] = a.force_att[1] = 0.0;
    a.force_alg[0] = a.force_alg[1] = 0.0;
    a.force_flee[0] = a.force_flee[1]=0.0;
    a.counter_rep = 0;
    a.counter_alg = 0;
    a.counter_att = 0;
    a.counter_flee = 0;
    
    consider_boundary(a, ptrSP);
    
    
     if ( a.steps_till_burst == 0 )
     {
        a.bin_step = ptrSP->burst_steps;
        double draw_update;
        unsigned int steps = 1;
        double burst_rate = ptrSP->burst_rate;
        unsigned int max_steps = 5 / (burst_rate * dt); // 1 burst takes on average 1/burst_rate times which are 1/ (burst_rate * dt) steps
        while(a.steps_till_burst == 0)
        {
            draw_update = gsl_rng_uniform(r);
            if (draw_update <= burst_rate * dt)
            {
                a.steps_till_burst = steps;
            }
            if (steps > max_steps )// break condition (no infinite loops)
            { 
                a.steps_till_burst = steps;   
            }
            steps += 1;
        }
    }
    
    consider_boundary(a, ptrSP);
}


template<class agent>
void Boundary(agent &a, double sizeL,  int BC)
{
// Function for calculating boundary conditions
// and update agent position and velocity accordingly
double tx;
double dphi=0.1;
double tmpphi;
double dx=0.0;
double dist2cen;
double diff;

std::vector<double> hv(2);
std::vector<double> wall_normal = a.x;

// -1 is Open boundary condition
switch (BC)
{
    case 0:
        // Periodic boundary condition
        a.x[0] = fmod(a.x[0] + sizeL, sizeL);
        a.x[1] = fmod(a.x[1] + sizeL, sizeL);
        break;
    // TODO: any case which changes the velocity should also change u, phi
    case 1:
        // Inelastic box boundary condition
        if(a.x[0]>sizeL)
        {
            a.x[0]=0.9999*sizeL;
            if(a.v[1]>0.)
                tmpphi=(0.5+dphi)*M_PI;
            else
                tmpphi=-(0.5+dphi)*M_PI;

            a.v[0]=cos(tmpphi);
            a.v[1]=sin(tmpphi);
        }
        else if(a.x[0]<0)
        {
            a.x[0]=0.0001;
            
            if(a.v[1]>0.)
            {
                tmpphi=(0.5-dphi)*M_PI;
            }
            else
            {
                tmpphi=-(0.5-dphi)*M_PI;
            }

            a.v[0]=cos(tmpphi);
            a.v[1]=sin(tmpphi);
        }
        if(a.x[1]>sizeL)
        {
            a.x[1]=0.9999*sizeL;
            if(a.v[0]>0.)
                tmpphi=-dphi;
            else
                tmpphi=M_PI+dphi;

            a.v[0]=cos(tmpphi);
            a.v[1]=sin(tmpphi);

        }
        else if(a.x[1]<0)
        {
            a.x[1]=0.0001;
            if(a.v[0]>0.)
                tmpphi=dphi;
            else
                tmpphi=M_PI-dphi;

            a.v[0]=cos(tmpphi);
            a.v[1]=sin(tmpphi);
        }

        break;
        
    case 2:
        // Elastic box boundary condition
        if(a.x[0]>sizeL)
        {
            dx=2.*(a.x[0]-sizeL);
            a.x[0]-=dx;
            a.v[0]*=-1.;
        }
        else if(a.x[0]<0)
        {
            dx=2.*a.x[0];
            a.x[0]-=dx;
            a.v[0]*=-1.;
        }
        if(a.x[1]>sizeL)
        {
            dx=2.*(a.x[1]-sizeL);
            a.x[1]-=dx;
            a.v[1]*=-1.;
        }
        else if(a.x[1]<0)
        {
            dx=2.*a.x[1];
            a.x[1]-=dx;
            a.v[1]*=-1.;
        }
        break;

    case 3:
        // Periodic boundary condition in x
        // Elastic  boundary condition in y
        tx= fmod(a.x[0], sizeL);
        if(tx < 0.0f)
            tx += sizeL;
        a.x[0]=tx;

        if(a.x[1]>sizeL)
        {
            dx=2.*(a.x[1]-sizeL);
            a.x[1]-=dx;
            a.v[1]*=-1.;
        }
        else if(a.x[1]<0)
        {
            dx=2.*a.x[1];
            a.x[1]-=dx;
            a.v[1]*=-1.;
        }
        break;
    case 4:
        // Periodic boundary condition in x
        // Inelastic  boundary condition in y
        tx= fmod(a.x[0], sizeL);
        if(tx < 0.0f)
        {
            tx += sizeL;
        }
        a.x[0]=tx;

        if(a.x[1]>sizeL)
        {
            dx=(a.x[1]-sizeL);
            a.x[1]-=dx+0.0001;
            a.v[1]=0.0;
        }
        else if(a.x[1]<0)
        {
            dx=a.x[1];
            a.x[1]-=dx-0.0001;
            a.v[1]=0.0;
        }
        break;
        
    case 5:
        // Elastic circle boundary condition
        dist2cen = vec_length(a.x);
        diff = dist2cen - sizeL;
        
        if (diff > 0)
        {
            // 1. mirror the position at the circular wall
            vec_mul221(wall_normal, 1 / dist2cen); // wall_normal = 1 * a.x / |a.x|
            hv = vec_mul(wall_normal, - 2 * diff); // hv = wall_normal * 2 * diff
            vec_add221(a.x, hv);
            // 2. mirror the velocity at the circular wall
            diff = vec_dot(wall_normal, a.v);
            hv = vec_mul(wall_normal, - 2 * diff);
            vec_add221(a.v, hv);
            // 3. update rest of agent properties
            a.u = vec_set_mag(a.v, 1);
            a.phi = atan2(a.u[1], a.u[0]);
        }
        break;
        
    case 6:
        // half-elastic circle boundary condition
        dist2cen = vec_length(a.x);
        diff = dist2cen - sizeL;
        
        if(diff > 0)
        {
            // 1. mirror the position at the circular wall
            vec_mul221(wall_normal, 1 / dist2cen); // wall_normal = 1 * a.x / |a.x|
            hv = vec_mul(wall_normal, - 2 * diff); // hv = wall_normal * 2 * diff
            vec_add221(a.x, hv);
            // 2. set velocity component normal to wall to 0
            diff = vec_dot(wall_normal, a.v);
            hv = vec_mul(wall_normal, -diff);
            vec_add221(a.v, hv);
            // vec_div221(a.v, 4);
            // 3. update rest of agent properties
            a.u = vec_set_mag(a.v, 1);
            a.phi = atan2(a.u[1], a.u[0]);
        }
        break;
    }
}

template
void Boundary(particle &a, double sizeL,  int BC);

template
void Boundary(predator &a, double sizeL,  int BC);


// Predators represent fishNet: 
//  -line of length < sizeL/2 (otherwise measures are not working)
//  -1.random direction selected
//  -2. identify com
//  -3. create fishNet sizeL/2 away from net in designated direction
void CreateFishNet(std::vector<particle> &a, params *ptrSP,
                   std::vector<predator> &preds, gsl_rng *r)
{
    std::vector<double> com(2);
    std::vector<double> hv(2);
    std::vector<int> empty(0);
   
    GetCenterOfMass(a, ptrSP, empty, com);
    
    double netLength = ptrSP->sizeL / 4;
    double phi = 2 * M_PI * gsl_rng_uniform(r);
    
    std::vector<double> v_net_move {cos(phi), sin(phi)};
    std::vector<double> v_net_perp = vec_perp(v_net_move);
    
    double dist2com  = 0.25;
    
    // dist2com = 0.25 * gsl_rng_uniform(r);
    std::vector<double> xpred_start = vec_mul(v_net_move, - ptrSP->sizeL * dist2com);
    vec_add221(xpred_start, com);   // thats the position L/2 away from COM
    
    double var_noise = M_PI / 4; 
    double noise = var_noise * gsl_rng_uniform(r) - var_noise / 2;
    
    phi += noise;
    v_net_move[0] = cos(phi);
    v_net_move[1] = sin(phi);
    v_net_perp = vec_perp(v_net_move);
    hv = vec_mul(v_net_perp, -netLength/2 );
    vec_add221( xpred_start, hv );
    
    double netDist = netLength / ( preds.size() - 1 );
    
    for(unsigned int i=0; i<preds.size(); i++)
    {
        hv = vec_mul(v_net_perp, i * netDist);
        preds[i].x = vec_add( xpred_start, hv );
        preds[i].phi = phi;
        preds[i].u = v_net_move;
        preds[i].v[0] = ptrSP->pred_speed0 * preds[i].u[0];
        preds[i].v[1] = ptrSP->pred_speed0 * preds[i].u[1];
        Boundary(preds[i], ptrSP->sizeL, ptrSP->BC);
    }
}


// Predators represent fishing-rot: 
//  -start maximum distance from COM
//  -random start angle
void CreatePredator(std::vector<particle> &a, params *ptrSP,
                    predator &pred, gsl_rng *r)
{
    // random direction:
    double phi = gsl_rng_uniform(r) * 2 * M_PI;
    
    pred.phi = phi;
    pred.u[0] = cos(phi);
    pred.u[1] = sin(phi);
    
    // sizeL/2 (approx max. distance) distance from com:
    std::vector<double> com(2);
    std::vector<double> hv = pred.u;
    std::vector<int> empty(0);
    
    GetCenterOfMass(a, ptrSP, empty, com);
    vec_mul221(hv, - ptrSP->sizeL/2);
    pred.x = vec_add(com, hv);
    pred.v[0] = ptrSP->pred_speed0 * pred.u[0];
    pred.v[1] = ptrSP->pred_speed0 * pred.u[1];
    Boundary(pred, ptrSP->sizeL, ptrSP->BC);
}


// Blind predator is supposed to have at different speeds the same persistence length:
//  -compare any speed v with v0=1
//      -agent with speed v needs dt=1m/v to travel 1m
//      -it is supposed to have the same angle deviation during 1m as a particle moving with v0
//          -Wiener-process= sqrt(dt) * Dphi
//              - sqrt(1/v0) * Dphi = sqrt(1/v) * \hat{Dphi}
//              - sqrt(1/1) *  Dphi = sqrt(1/v) * \hat{Dphi}
//              -> \hat{Dphi} = sqrt(v) * Dphi
//  -thus must the standard deviation for the angular noise increase with sqrt(v)
void MovePredator(predator &pred, std::vector<particle> &a, params *ptrSP, gsl_rng *r)
{
    double lphi;
    
    if(ptrSP->pred_move == 0)
    {
        double rnp = ptrSP->noisep * gsl_ran_gaussian(r, 1.0); //  ptrSP->noisep * gsl_ran_gaussian(r, 1.0);
        lphi = pred.phi;
        lphi +=  rnp * sqrt(ptrSP->pred_speed0); // to keep persistence length the same
        lphi = fmod(lphi, 2*M_PI);
    }
    // directly going for closest prey (no special dynamics)
    else
    {
        std::vector<double> hv(2);
        std::vector<double> minvec = CalcDistVec(pred.x, a[0].x, ptrSP->BC, ptrSP->sizeL); // pointing to a
        
        double mindist = vec_length(minvec);
        
        for(int i=1; i<a.size(); i++)
        {
            hv = CalcDistVec(pred.x, a[i].x, ptrSP->BC, ptrSP->sizeL); // pointing to a
            double dist = vec_length(hv);
            
            if (dist < mindist)
            {
                mindist = dist;
                minvec = hv;
            }
        }
        lphi = atan2(minvec[1], minvec[0]);
        // std::vector<double> com(2);
        // std::vector<int> empty(0);
        // GetCenterOfMass(a, ptrSP, empty, com);
        // vec_sub221(com, pred.x);
        // lphi = atan2(com[1], com[0]);
    }
    pred.phi = lphi;
    pred.u[0] = cos(lphi);
    pred.u[1] = sin(lphi);
    pred.v[0] = ptrSP->pred_speed0 * pred.u[0];
    pred.v[1] = ptrSP->pred_speed0 * pred.u[1];
    pred.x[0] = pred.x[0] + pred.v[0] * ptrSP->dt;
    pred.x[1] = pred.x[1] + pred.v[1] * ptrSP->dt;
    Boundary(pred, ptrSP->sizeL, ptrSP->BC);
    // move to com:
}


void MoveFishNet(std::vector<predator> &preds, params *ptrSP)
{
    for(unsigned int i=0; i<preds.size(); i++)
    {
        preds[i].x[0] = preds[i].x[0] + preds[i].v[0] * ptrSP->dt;
        preds[i].x[1] = preds[i].x[1] + preds[i].v[1] * ptrSP->dt;
        Boundary(preds[i], ptrSP->sizeL, ptrSP->BC);
    }
}


double predictV(double v, double force,
                double time, double friction)
{
    // predicts the velocity along a selected position, assuming a constant force.
    // it is derived from the equations of motion used in this model:
    // dv/dt = -friction * v + F
    // The solution of this differential equation is:
    // v(t) = (v0 - F/friction) * exp(-friction * t) + F/friction
    // with v0 as velocity at time t=0 and F as the force.
    double h = force / friction;
    double v_pred = (v - h) * exp(-friction * time) + h;
    return v_pred;
}


double predictPos(double v, double force,
                double time, double friction)
{
    // predicts the position along a selected position, assuming a constant force.
    // it is dierived from the differential equation:
    // dx/dt = v
    // with v as in the function "predictV" derived:
    // v(t) = (v0 - F/friction) * exp(-friction * t) + F/friction
    // the solution of this DEQ is:
    // r(t) = (v0/friction * F/frictin^2) * ( 1-exp(-friction * t) ) + F/friction * t
    double h = force / friction;
    double x_pred = (v/friction - h/friction) * ( 1 - exp(-friction * time) ) + h * time;
    
    return x_pred;
}


std::vector<double> predictXatNextBurst(std::vector<double> & x,
                                        std::vector<double> & v,
                                        std::vector<double> & force,
                                        double t_burst, double t_tnb,
                                        double friction)
{
    // the function name should actually be "predict position a little after the next burst
    // assuming the next burst would be delayed a little".
    // Obviously too long.
    // This delay, the little time after, is necessary because if the agent is 
    // at the next burst exactly at, or very close to the tank wall,
    // it could not adjust its force to avoid the collision.
    // Therefore I predict a little longer.
    // The little longer is exaclty the time the agent normally bursts.
    // If an agent would still be inside, if it would coast instead of bursting,
    // than the burst force at the next burst is for sure sufficient to avoid collision.
    unsigned int dim = x.size();
    
    std::vector<double> x_anb(dim, 0); // anb = at next burst
    
    double tb = fmin(t_burst, t_tnb); // burst time
    double tc = t_burst; // coast always t_burst longer see explanation above 
    
    if(t_tnb > t_burst)
    {
        tc += ( t_tnb - tb);
    }
    
    for(unsigned int i=0; i<dim; i++)
    {
        x_anb[i] = predictPos(v[i], force[i], tb, friction);
        double vb = predictV(v[i], force[i], tb, friction); // velocity after burst and before coast
        x_anb[i] += predictPos(vb, 0., t_tnb - tb, friction);
        x_anb[i] += x[i];
    }
    
    return x_anb;
}


double forceChange2NotCollide(particle &a, params * ptrSP,
                              double dphi)
{
    // change in force direction results in an inside position
    // 1. It starts with an change of dphi
    // 2. it repeatedly changes the force direction by dphi until it 
    //      finds the first position inside the boundary
    // 3. Now at each next step dphi is halfed 
    // 4. if an an inside position is found it changes the force by -dphi
    //      i.e. it makes the change smaller
    // 5. if an an outside position is found it changes the force by +dphi
    // 6. 2 break conditions exists: final resolution reached OR change = PI
    double t_burst = ptrSP->dt * ptrSP->burst_steps;
    double t_tnb = ptrSP->dt * a.steps_till_burst;
    double friction = ptrSP->beta;
    double force_mag = vec_length(a.force);
    double forceAngle = atan2(a.force[1], a.force[0]);
    double insideChange = M_PI;
    
    std::vector<double> force(2);
    std::vector<double> x_anb(2);
    
    bool outside = true;
    bool foundInsideAngle = false;
   
    double forceChange = dphi;
    
    while (fabs(forceChange) < M_PI and outside)
    {
        force[0] = force_mag * cos(forceAngle + forceChange);
        force[1] = force_mag * sin(forceAngle + forceChange);
        x_anb = predictXatNextBurst(a.x, a.v, force, t_burst, t_tnb, friction);
        
        double r_x = vec_length(x_anb);
        
        if(r_x < ptrSP->sizeL)
        {
            insideChange = forceChange;
            foundInsideAngle = true;
            outside = false;
        }
        if (foundInsideAngle)
        {
            dphi /= 2;
        }
        if(outside)
        {
            forceChange += dphi;
        }
        else if(fabs(dphi) > M_PI/64)
        {
            forceChange -= dphi;
            outside = true;
        }
    }
    
    return insideChange;
}


double closestForceDirection(particle &a, params * ptrSP)
{
    double forceAngle = atan2(a.force[1], a.force[0]);
    double angularChangeCcw = forceChange2NotCollide(a, ptrSP, M_PI/4);
    double angularChangeCw = forceChange2NotCollide(a, ptrSP, -M_PI/4);
    double angle = angularChangeCcw;
    
    if(fabs(angularChangeCw) < fabs(angle))
    {
        angle = angularChangeCw;
    }
    
    return angle + forceAngle;
}


// killing of individuals if in kill-range
//      to speed up computation: instead of NxN_p distance computations 
//                               compute if agent in front of net and closer than kill_range
//      assumptions: -sizeL/2 < NetLength
//                   -net is moving perpendicular to its elongation
void FishNetKill(std::vector<particle> &a, std::vector<predator> &preds,
                 params *ptrSP)
{
    std::vector<double> v_net = CalcDistVec(preds[0].x, preds[1].x, 
                                            ptrSP->BC, ptrSP->sizeL);
                                            
    double dist = vec_length(v_net);
    double NetLength = ( preds.size() -1 ) * dist;
    
    vec_div(v_net, dist);
    
    std::vector<double> v_net_move = preds[0].u;
    
    // DEBUGGING
    double zero = fabs(vec_dot( v_net, v_net_move ));
    if ( zero > 0.0001 )
    {
        std::cout<< "zero: " << zero << std::endl;
    }
    if(ptrSP->sizeL / 2 < NetLength)
    {
        std::cout<< "ptrSP->sizeL / 2 < NetLength: FishNetKill not working" << std::endl;
    }
    // compute projection of distance 
    std::vector<double> r_jp(2);
    
    double front, side; 
    
    for(unsigned int i=0; i<a.size(); i++)
    {
        r_jp = CalcDistVec(preds[0].x, a[i].x, ptrSP->BC, ptrSP->sizeL); // r_jp = a.x - pred.x -> pointing to a
        front = vec_dot(v_net_move, r_jp);
        side = vec_dot(v_net, r_jp);
        
        if(front > 0 && front < ptrSP->kill_range && side > 0 && side <= NetLength)
        {
            a[i].dead = true;
            ptrSP->Ndead++;
        }
    }
}


// killing of individuals if in kill-range (kill_mode=1)
//                        if in kill-range/sqrt(N_s) (kill_mode=2, confusion)
//                              with N_s as # of agents sensed (r_ip < 2*kill_range)
void PredKill(std::vector<particle> &a, predator &pred,
              params *ptrSP, gsl_rng *r)
{
    double dist;
    double r_sense = 4 * ptrSP->kill_range;
    
    if(ptrSP->pred_kill == 1)
    {
        for(unsigned int i=0; i<a.size(); i++)
        {
            dist = CalcDist(pred.x, a[i].x, ptrSP->BC, ptrSP->sizeL);
            
            if(dist <= ptrSP->kill_range)
            {
                a[i].dead = true;
                ptrSP->Ndead++;
            }
        }
    }
    else if(ptrSP->pred_kill > 1)
    {
        std::vector<double> catch_probs;
        
        catch_probs.reserve(a.size());
        
        int N_sensed = 0;
        
        std::vector<unsigned int> possible_kill_id;
        
        possible_kill_id.reserve(a.size());
        
        double prob_catched;
        double total_catch_prob = 0;
        
        for(unsigned int i=0; i<a.size(); i++)
        {
            dist = CalcDist(pred.x, a[i].x, ptrSP->BC, ptrSP->sizeL);
            
            if(dist <= r_sense)
            {
                N_sensed++;
                if(dist < ptrSP->kill_range)
                {
                    possible_kill_id.push_back(i);
                    prob_catched = (ptrSP->kill_range - dist) / ptrSP->kill_range; // is always > 0 due to if-condition above
                    catch_probs.push_back( prob_catched );
                    total_catch_prob += prob_catched;
                }
            }
        }
        
        double effective_kill_rate;
        
        // confusion effect assumed to be 1/N_sensed
        if(N_sensed > 0)
        {
            effective_kill_rate = ptrSP->kill_rate * SigmoidFunction( N_sensed, ptrSP->N_confu, -1);
        }
        double prob_killed = 0;
        
        for(unsigned int i=0; i<possible_kill_id.size(); i++)
        {
			// probabilistic 
            if (ptrSP->pred_kill == 2)
            {
                prob_killed = ptrSP->kill_rate * ptrSP->dt;
            }
            // probabilistic + confusion
            else if (ptrSP->pred_kill == 3)
            {
                prob_killed = effective_kill_rate * ptrSP->dt;
            }
            // probabilistic + confusion + selection 
            else if (ptrSP->pred_kill == 4)
            {
                double prob_selected = catch_probs[i] / total_catch_prob;
                prob_killed = effective_kill_rate * ptrSP->dt * prob_selected;
            }
            
            double luck = gsl_rng_uniform(r); 
            
            if(prob_killed > luck)
            {
                unsigned int ii = possible_kill_id[i];
                a[ii].dead = true;
                ptrSP->Ndead++;
            }
        }
    }
}

double force2stop(double v, double t_burst, double friction)
{
    // this function computes the force to decelerate the agent along the
    // selected direction of motion to zero.
    // It is derived by setting v(t) to zero:
    // v(t) = 0 = (v0 - F) * exp(-friction * t) + F/friction
    // which results in:
    // F = - ( v0 * friction * exp(-friction * t) ) / ( 1 - exp(-friction * t) ) 
    double h = exp(- friction * t_burst);
    double force = - friction * v * h / ( 1 - h );
    
    return force;
}
