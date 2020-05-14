/*  SwarmDynamics 
    Stochastic differential equation model of agents interacting via attraction/repulsion/alignment.
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
#include "swarmdyn.h"

// GLOBAL STRUCT DECLARATION


// GLOBAL VARIABLE DECLARATION
long seed;          // seed rng

// gsl random number generator
gsl_rng *r; // global generator
const gsl_rng_type * T;

//Runs with ./swarmdyn
int main(int argc, char **argv){

    int s;              // variable for current step number eg: currenttime/dt

    params SysPara;   // system

    // Init and output the system paramters
    SetCoreParameters(&SysPara);
    ParseParameters(argc, argv, &SysPara);
    InitSystemParameters(&SysPara);
    OutputParameters(SysPara);

    // initialize agents and set initial conditions
    std::vector<particle> agent(SysPara.N);    // particles or prey
    std::vector<particle> agent_dead(0);
    InitSystem(agent, SysPara);
    InitRNG();
    ResetSystem(agent, &SysPara, false, r);

    double dt = SysPara.dt;
    SysPara.r = r;

    std::vector<predator>  preds(SysPara.Npred);
    InitPredator(preds);

    int sstart = 0;
    double t1 = clock(), t2 = 0.; // time variables for measuring comp. time
    std::cout<< "Go";
    // Perform numerical integrate
    t1 = clock();
    for(s=sstart; s < SysPara.sim_steps; s++){
        // define some basic time-flags
        bool time_pred = (s >= static_cast<int>(SysPara.pred_time/dt));
        bool time_output = (s >= static_cast<int>(SysPara.trans_time/dt));
        // Perform a single step
        // first split: output handles agents who are dead but
        //  NN of non-dead agents (also imporant for NN2)
        split_dead(agent, agent_dead, preds);
        if (agent.size() == 0){
            Output(s, agent, SysPara, preds, true);
            break;
        }
        Step(s, agent, &SysPara, preds);
        // Data output
        if(s%SysPara.step_output==0 && time_output)
        {
            // UNCOMMENT if simulation per output step are interesting:
            // t2= clock();
            // double tdiff=(t2-t1)/CLOCKS_PER_SEC;
            // printf("s=%d; time / out. step = %.4f\n",s,tdiff);
            // t1=t2;
            Output(s, agent, SysPara, preds);
            SysPara.outstep += 1;
            if (time_pred)
                SysPara.outstep_pred += 1;
        }
    }
    // if minimum output generated -> assumes equilibration run
    // -> save final positions velocities
    merge_dead(agent, agent_dead);
    if (SysPara.outstep == 1)
        WritePosVel(agent, &SysPara, "final_posvel_" + SysPara.fileID, false);
    return 0;
}


long unsigned int getseed(int const K)
{

    typedef std::chrono::high_resolution_clock hiclock;

    auto gett= [](std::chrono::time_point<hiclock> t0)
    {
        auto tn = hiclock::now();
        return static_cast<long unsigned int>(std::chrono::duration_cast<std::chrono::microseconds>(tn-t0).count());
    };

    long unsigned int diffs[10];
    diffs[0] = gett(hiclock::now());
    for(int i=1; i!=10; i++)
    {
        auto last = hiclock::now();
        for(int k=K; k!=0; k--)
        {
            diffs[i]= gett(last);
        }
    }

    return *std::max_element(&diffs[1],&diffs[9]);
}


void InitRNG(){
    // Initialize random number generator
    // time_t  t1;                     // Get system time for random number seed
    // time(&t1);
    // seed=t1;
    // std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    // seed=std::chrono::time_cast<long int>(t1);
    seed = static_cast<unsigned long>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
#if PRINT_PARAMS
    printf("Random number seed: %d\n",static_cast<int>(seed));
#endif
    gsl_rng_env_setup();
        T=gsl_rng_default;
        r=gsl_rng_alloc(T);
    gsl_rng_set(r, seed);
}

void Step(int s, std::vector<particle> &a, params* ptrSP, std::vector<predator> &preds)
{
    // function for performing a single (Euler) integration step
    int i = 0;
    double rnp = 0.0;

    double dt = ptrSP->dt;
    int N = a.size();
    unsigned int j;
    int ii;

    // CREATE PREDATOR
    if (s>=ptrSP->pred_time/dt && s-1<ptrSP->pred_time/dt){
        if ( preds.size() == 1 )
            CreatePredator(a, ptrSP, preds[0], r);
        else
            CreateFishNet(a, ptrSP, preds, r);
    }
    // reCreate Fishnet: evertime the net has passed
    if ( preds.size() > 1 ){
        int steps_since_pred = s - static_cast<int>(ptrSP->pred_time/ptrSP->dt);
        int steps_needed;
        if (ptrSP->pred_speed0 > 0)
            steps_needed = static_cast<int>( ptrSP->sizeL / 2 /
                                            (ptrSP->pred_speed0 * ptrSP->dt) );
        else
            steps_needed = static_cast<int>( (ptrSP->sim_time - ptrSP->trans_time) / 4 / ptrSP->dt);
        if ( steps_since_pred % steps_needed == 0 )
            CreateFishNet(a, ptrSP, preds, r);
    }
    // Reset simulation-step specific values to default
    for (i=0; i<N; i++)
        a[i].NN.resize(0);
    for (i=0; i<preds.size(); i++){
        preds[i].NNset.clear();
        preds[i].NN.resize(0);
    }
    // INTERACTION:
    if (s < ptrSP->pred_time/dt)
        InteractionVoronoiF2F(a, ptrSP);
    else
        InteractionVoronoiF2FP(a, ptrSP, preds); // has build in pred->voronoinn computation

    // Update all agents
    for(i=0;i<N;i++)
    {
        // Generate noise
        rnp = ptrSP->noisep * gsl_ran_gaussian(r, 1.0);
        ParticleBurstCoast(a[i], ptrSP, r);
    }
    // PREDATOR RELATED STUFF(P-move, .... )
    if (s>=ptrSP->pred_time/dt){
        if ( preds.size() > 1 ){
            MoveFishNet(preds, ptrSP);
            if ( ptrSP->pred_kill != 0 ) // predator kills if in kill_range (kill_range/sqrt(N_det))
                FishNetKill(a, preds, ptrSP);
        }
        else{
            MovePredator(preds[0], a, ptrSP, r);
            if ( ptrSP->pred_kill != 0 ) // predator kills if in kill_range (kill_range/sqrt(N_det))
                PredKill(a, preds[0], ptrSP, r);
        }
    }
}


void Output(int s, std::vector<particle> &a, params &SP,
            std::vector<predator> &preds, bool forceSave){
    std::vector<double> out;

    if (s < SP.pred_time / SP.dt){
        if (SP.out_mean){
            out = Out_swarm(a, SP);
            DataCreateSaveWrite(SP.dataOutMean, out, SP,
                                "swarm", forceSave);
        }
        if (SP.out_particle)
            WriteParticles<particle>(a, SP, "part", SP.outstep);
    }
    else{
        if (SP.out_mean){
            out = Out_swarm(a, SP);
            DataCreateSaveWrite(SP.dataOutSwarm, out, SP,
                                "swarm", forceSave);
            out = Out_swarm_fishNet(a, preds, SP);
            DataCreateSaveWrite(SP.dataOutSwarmPred, out, SP, 
                                "swarm_fishNet", forceSave);
        }
        if (SP.out_particle){
            WriteParticles<particle>(a, SP, "part", SP.outstep);
            WriteParticles<predator>(preds, SP, "pred", SP.outstep_pred);
        }
    }
}


std::vector<double> Out_swarm(std::vector<particle> &a, params &SP){
    int N = a.size();
    double dist = 0;    // distance between predator and single prey
    std::vector<double> avg_x(2);   // average position vector
    // helper
    double hd = 0;
    std::vector<double> hv(2, 1);
    // to compute the aspect ratio:
    double max_IID = 0;
    std::vector<double> max_IID_vec(2, 1);
    // compute averages
    std::vector<int> allprey(N);
    std::iota (std::begin(allprey), std::end(allprey), 0); //Fill with 0, 1,...N
    std::vector<double> basic_avgs = basic_particle_averages(a, allprey);
    std::vector<double> avg_v {basic_avgs[0], basic_avgs[1]};
    std::vector<double> avg_u {basic_avgs[2], basic_avgs[3]};
    double avg_vsquare = basic_avgs[5];
    double avg_speed = basic_avgs[4];

    double IID = 0;     // Inter Individual Distance
    double NND = 0;     // Nearest Neighbor Distance
    double ND = 0;      // Neighbor Distance
    double nd = 0;
    for(int i=0; i<N; i++){
        vec_add221(avg_x, a[i].x);
        // NND:
        nd = 0;
        for (auto it=a[i].NN.begin(); it!=a[i].NN.end(); ++it){
            dist = CalcDist(a[i].x, a[*it].x, SP.BC, SP.sizeL);
            nd += dist;
            hd = fmin(hd, dist);
        }

        nd /= a[i].NN.size();
        ND += nd;

        // NND: (not in NN-loop because async-update do not has always NN)
        hd = SP.N * SP.alg_range;  // arbitrary large value
        for(int j=0; j<i-1; j++){
            dist = CalcDist(a[i].x, a[j].x, SP.BC, SP.sizeL);
            hd = fmin(hd, dist);
        }
        // IID:
        for(int j=i+1; j<N; j++){
            hv = CalcDistVec(a[i].x, a[j].x, SP.BC, SP.sizeL);
            dist = vec_length(hv);
            hd = fmin(hd, dist); // for NND
            IID += dist;
            if(dist > max_IID){ // maxIID and vector for aspect_ratio
                max_IID = dist;
                max_IID_vec = hv; 
            }
        }
        NND += hd;

    }
    vec_div221(avg_x, N);
    NND /= N;
    ND /= N;
    IID /= (N-1) * N;
    double avgvel = vec_length(avg_v);
    double pol_order = vec_length(avg_u);
    // double avgdirection = atan2(avg_v[1], avg_v[0]);

    // Calc normalized angular momentum (normalized by radial distance)
    // milling OP:
    double L_norm = 0;
    for(int i=0; i<N; i++){
        hv[0] = a[i].x[0] - avg_x[0];
        hv[1] = a[i].x[1] - avg_x[1];
        L_norm += (hv[0] * a[i].v[1] - hv[1] * a[i].v[0]) / (vec_length(hv));   // L/r=(\vec{r} x \vec{v})/r
    }
    L_norm = fabs(L_norm) / N;

    double Area_ConHull = AreaConvexHull(a, allprey);
    // computes elongation and aspect ratio:
    double elongation = get_elongation(a, avg_u, allprey);
    double aspect_ratio = get_elongation(a, max_IID_vec, allprey);
    double a1, a2, a_com_maxIID;
    a1 = (avg_v[0] * max_IID_vec[0] +
          avg_v[1] * max_IID_vec[1]) /
         (avgvel * vec_length(max_IID_vec));
    a2 = acos(-a1);
    a1 = acos(a1);
    a_com_maxIID = fmin(a1, a2); 

    // generate output
    std::vector<double> out_vec;
    out_vec.reserve(17);
    out_vec.push_back( N); 
    out_vec.push_back( pol_order); 
    out_vec.push_back( L_norm); 
    out_vec.push_back( avg_x[0]); 
    out_vec.push_back( avg_x[1]); 
    out_vec.push_back( avg_v[0]); 
    out_vec.push_back( avg_v[1]); 
    out_vec.push_back( avg_speed); 
    out_vec.push_back( avg_vsquare); 
    out_vec.push_back( elongation); 
    out_vec.push_back( aspect_ratio); 
    out_vec.push_back( a_com_maxIID); 
    out_vec.push_back( Area_ConHull); 
    out_vec.push_back( NND); 
    out_vec.push_back( IID); 
    out_vec.push_back( ND); 

    return out_vec;
}




// output for preds.size() > 1: predator resemble a fishing net
// Most important: 
// -count how many in front and in kill-range
// -prey distance to net = min(prey distance to predators)
// -FRONT prey distance to net = min(front prey distance to predators)
//      - front-prey difficult for periodic BC -> Fitness also useless
//      - alternativ: fish-net swimming always up, fish only endangered if in front and close-> killing range
//         pred.NNset = prey seeing predator -> have normal repulsion 
//                                              have fitness decrease if 
std::vector<double> Out_swarm_fishNet(std::vector<particle> &a,
                                      std::vector<predator> &preds, params &SP)
{
    // ASSUMING: - fish-net on right-side (pred.x[0] >= sizeL/2)
    //           - fish-net moving up (pred.v[0] = 0)
    //           - fish-net is horizontal line (preds[0].x[1] = preds[1].x[1] = preds[3].x[1])
    // # of prey in front and in kill_range
    int NfrontClose = 0;
    int Nfront = 0;
    double yfront = preds[0].x[1];
    double NetLength = 0;
    double dist;
    std::vector<double> r_jp(2);
    if (preds.size() > 1){
        std::vector<double> v_net_move = preds[0].u;
        std::vector<double> v_net = CalcDistVec(preds[0].x, preds[1].x, 
                                                SP.BC, SP.sizeL);
        dist = vec_length(v_net);
        NetLength = ( preds.size() -1 ) * dist;
        double front, side; 
        for (unsigned int i=0; i<a.size(); i++){
            r_jp = CalcDistVec(preds[0].x, a[i].x, SP.BC, SP.sizeL); // r_jp = a.x - pred.x -> pointing to a
            front = vec_dot(v_net_move, r_jp);
            side = vec_dot(v_net, r_jp);
            if ( side > 0 && side <= NetLength ){
                Nfront++;
                if ( front > 0 && front < SP.kill_range )
                    NfrontClose++;
            }
        }
    }
    // prey-distance to net (very general)
    double distMin = SP.sizeL * 2;
    std::vector<double> hvec(2);
    double distNet = 0;
    for (int i=0; i<a.size(); i++){
        for (int j=0; j<preds.size(); j++){
            hvec = CalcDistVec(preds[j].x, a[i].x, SP.BC, SP.sizeL);
            dist = vec_length(hvec);
            distMin = fmin(distMin, dist);
        }
        distNet += distMin;
        distMin = SP.sizeL * 2; // reset maximum distance
    }
    distNet /= a.size();
    // generate output
    std::vector<double> out_vec;
    out_vec.reserve(4);
    out_vec.push_back( NfrontClose);
    out_vec.push_back( distNet);
    out_vec.push_back( Nfront);
    out_vec.push_back( SP.Ndead);
    return out_vec;
}

std::vector<double> Out_swarm_pred(std::vector<particle> &a,
                                   predator &pred, params &SP){
    std::vector<double> basic_avgs;
    std::vector<int> allprey(a.size());
    std::iota (std::begin(allprey), std::end(allprey), 0); //Fill with 0, 1,...N
    // CLU  ############################################################
    int N_clu = a.size();
    // averaging
    basic_avgs = basic_particle_averages(a, allprey);
    std::vector<double> clu_avg_v {basic_avgs[0], basic_avgs[1]};
    double clu_avg_s = basic_avgs[4];
    double clu_avg_vsquare = basic_avgs[5];
    double clu_dpi = get_avg_pred_dist(a, allprey, pred, SP);
    // generate output
    std::vector<double> out_vec;
    out_vec.reserve(5);
    out_vec.push_back( N_clu);
    out_vec.push_back( clu_avg_s);
    out_vec.push_back( clu_avg_vsquare);
    out_vec.push_back( clu_dpi);
    out_vec.push_back( SP.Ndead);
    return out_vec;
}


void DataCreateSaveWrite(std::vector< std::vector<double> > &data,
                     std::vector<double> &out, params &SP,
                     std::string name, bool forceSave){
    unsigned int n = data.size();
    unsigned int m = out.size();
    // CREATE (initialize data-vector):
    if (n == 0){
        unsigned int left_outsteps = SP.total_outstep - SP.outstep;
        std::vector< std::vector<double> > vec2d(left_outsteps, std::vector<double>(m));
        data = vec2d;
        n = left_outsteps;
    }
    // SAVE current values to data-array
    unsigned int current = SP.outstep - (SP.total_outstep - n);
    data[current] = out;
    // WRITE TO FILE (last output-step)
    if ( (SP.outstep == SP.total_outstep - 1) || forceSave ){
        if (SP.out_h5){
            std::string f_h5out = SP.location + "out_" + SP.fileID + ".h5";
            h5CreateWriteDset(f_h5out, SP.fileID, name, data, SP.out_extend);
        }
        else
            WriteVector2d(SP.location + name + "_" + SP.fileID + ".dat",
                          data, false);
    }
}


/*Write_out: takes a std::vector and writes its content to 
* a hdf5-dataset (h5dset) or to a file (name) depending
* on the output mode (SP.out_h5)
*/
void Write_out(std::vector<double> &out, params &SP,
           std::vector<hsize_t> & vec_dim,
           H5::DataSet *h5dset, std::string name){
    if (SP.out_h5){
        // create offset of data 
        std::vector<hsize_t> offset(vec_dim.size());  // offset of hyperslab  
        if (SP.out_extend)
            offset[0] = vec_dim[0] - 1;  // to not overwrite preceding run
        offset[offset.size() - 2] = SP.outstep;      // corresponds to time offset
        h5WriteDouble(h5dset, out, offset);
    }
    else
        WriteVector(SP.location + name + "_" + SP.fileID + ".dat", out, false);
}


template<class T>
std::vector<double> basic_particle_averages(std::vector<particle> &a,
                                            std::vector<T> &nodes){
    int N = nodes.size();
    //initialize
    unsigned int ii = 0;
    std::vector<double> avg_v(2), avg_u(2);
    double vsquare, avg_s, avg_vsquare;
    avg_s = avg_vsquare = 0;
    // compute output
    for(int i=0; i<N; i++){
        ii = nodes[i];
        vec_add221(avg_v, a[ii].v);
        vec_add221(avg_u, a[ii].u);
        vsquare = a[ii].v[0] * a[ii].v[0] + a[ii].v[1] * a[ii].v[1];
        avg_s += sqrt(vsquare);
        avg_vsquare += vsquare;
    }
    vec_div221(avg_v, N);
    vec_div221(avg_u, N);
    avg_s /= N;
    avg_vsquare /= N;
    // return output
    std::vector<double> out {avg_v[0], avg_v[1], avg_u[0], avg_u[1], avg_s, avg_vsquare};
    return out;
}
template
std::vector<double> basic_particle_averages(std::vector<particle> &a,
                                            std::vector<int> &nodes);
template
std::vector<double> basic_particle_averages(std::vector<particle> &a,
                                            std::vector<unsigned int> &nodes);

template<class T>
double get_avg_pred_dist(std::vector<particle> &a,
                         std::vector<T> &nodes, predator &pred,
                         params &SP){
    T N = nodes.size();
    unsigned int ii;
    double dist = 0;
    double dpi;
    std::vector<double> r_pi(2);
    for(T i=0; i<N; i++){
        ii = nodes[i];
        r_pi = CalcDistVec(pred.x, a[ii].x, SP.BC, SP.sizeL);
        dpi = vec_length(r_pi);
        dist += dpi;
    }
    dist /= N;
    return dist; 
}
template
double get_avg_pred_dist(std::vector<particle> &a,
                         std::vector<int> &nodes, predator &pred,
                         params &SP);
template
double get_avg_pred_dist(std::vector<particle> &a,
                         std::vector<unsigned int> &nodes, predator &pred,
                         params &SP);


template<class part>
void WriteParticles(std::vector<part> &a, params &SP, 
                    std::string name, double outstep){
    if (a.size() == 0)
        return;
    std::vector<double> out = a[0].out();
    if (SP.out_h5){
        H5::DataSet *h5dset;
        // load/create h5-file
        std::string f_h5out = SP.location + "out_" + SP.fileID + ".h5";
        H5::H5File *H5out;
        if ( exists(f_h5out) )
            H5out = new H5::H5File(f_h5out.c_str(), H5F_ACC_RDWR);
        else
            H5out = new H5::H5File(f_h5out.c_str(), H5F_ACC_TRUNC);
        std::vector<hsize_t> dim(3, 0);
        // if outstep == 0: create OR extend-dataset + read-dimension
        //      extend: easy... just extend 0 dimension
        //      no-extend: creat with dim(time=total_outstep-SP.outstep, N, out.size())
        std::string n_dset;
        if (SP.fileID == "xx")
            n_dset = "/" + name;
        else
            n_dset = "/" + SP.fileID + "/" + name;
        if ( outstep == 0){
            if ( SP.out_extend ){    // multiple Runs -> open and extend existing files
                h5dset = new H5::DataSet(H5out->openDataSet(n_dset.c_str()));
                h5read_extend_dim(h5dset, dim);
            }
            else{
                dim[0] = SP.total_outstep - SP.outstep; // is less/equal for particleD
                dim[1] = a.size();  // could be problematic if agents already killed -> BUG
                dim[2] = out.size();
                if ( (dim[1] < SP.N) && (dim[1] != SP.Npred) ) // assume that only prey get killed -> must be prey
                    dim[1] = SP.N;
                h5dset = new H5::DataSet(h5CreateDSet(H5out, dim,
                                                      n_dset.c_str(), "double"));
            }
        }
        // else: load dataSet and read dimension
        else{
            h5dset = new H5::DataSet(H5out->openDataSet(n_dset.c_str()));
            h5readDimension(h5dset, dim);
        }
        // now create-offset
        std::vector<hsize_t> offset(dim.size());  // offset of hyperslab  
        if (SP.out_extend)
            offset[0] = dim[0]-1;  // to not overwrite preceding run
        offset[dim.size()-3] = outstep;  // time-offset
        // write data to offset
        for(int i=0; i<a.size(); i++){
            out = a[i].out();
            offset[offset.size()-2] = a[i].id;  // corresponds to particle offset
            h5WriteDouble(h5dset, out, offset);
        }
        delete h5dset;
        delete H5out;
        // h5dset->close();
        // H5out->close();
    }
    else {
        std::ofstream outFile((SP.location + name + "_" + SP.fileID
                               + ".dat").c_str(), std::ios::app);
        std::vector<double> out_default(out.size(), 0);
        unsigned int id = 0;
        for(int i=0; i<a.size(); i++){
            out = a[i].out();
            while (id < a[i].id){
                for (int j=0; j<out_default.size(); j++)
                    outFile << out_default[j] << " ";
                outFile << std::endl;
                id++;
            }
            for (int j=0; j<out.size(); j++)
                outFile << out[j] << " ";
            outFile << std::endl;
            id++;
        }
    }
}

template
void WriteParticles(std::vector<particle> &a, params &SP, 
                    std::string name, double outstep);
template
void WriteParticles(std::vector<predator> &a, params &SP, 
                    std::string name, double outstep);
