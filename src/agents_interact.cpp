/*  AgentsInteract
    computes interactions between Agents based on metric or
    voronoi for SwarmDynamics(swarmdyn.h, swarmdyn.cpp)
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

#include "agents_interact.h"


void InteractionVoronoiF2F(std::vector<particle> &a, params *ptrSP)
{
    // calculates local voronoi interactions
    typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
    typedef CGAL::Triangulation_vertex_base_with_info_2<int, K>         Vb;
    typedef CGAL::Triangulation_data_structure_2<Vb>                    Tds;
    typedef CGAL::Delaunay_triangulation_2<K, Tds>                      Delaunay;
    typedef Delaunay::Vertex_handle                                     Vertex_handle;
    typedef Delaunay::Edge_iterator                                     Edge_iterator;
    typedef Delaunay::Point                                             Point;
    typedef std::pair<Point, int>                                       PPoint;

    int N = a.size();

    // create the delaunay triangulation network
    std::vector< std::pair<Point, int> > Vr;   // stores the point location and index
    std::vector< std::pair< std::vector<double>, int > > posId;
    Vr.reserve(4 * N);
    posId.reserve(4 * N);
    for(int i=0; i<N; i++){
        Point p1(a[i].x[0], a[i].x[1]);
        PPoint p2 = std::make_pair(p1, i);
        Vr.push_back(p2);
        makePairAndPushBack(posId, a[i].x, i);
    }
    // produce replicate prey for periodic BC
    if (!ptrSP->BC){
        std::vector< std::pair< std::vector<double>, int > > newPosId = 
            GetCopies4PeriodicBC( posId, ptrSP->sizeL );
        for(int i=0; i<newPosId.size(); i++){
            Point p1(newPosId[i].first[0], newPosId[i].first[1]);
            PPoint p2 = std::make_pair(p1, newPosId[i].second);
            Vr.push_back(p2);
        }
    }
    // do delauney-triangulation
    Delaunay t;
    t.insert(Vr.begin(),Vr.end());

    //iterates over all finite edges and apply interaction 
    //  (infinite edges connect the infinite vertex with the vertices of the complex hull)
    for(Edge_iterator ei=t.finite_edges_begin(); ei!=t.finite_edges_end(); ei++){
      // Get a vertex from the edge, edge is stored as pair of the neighboring face and the vertex opposite to it
      Delaunay::Face& f = *(ei->first);
      int i = ei->second;
      Vertex_handle vi = f.vertex(f.cw(i));     // cw = clockwise rotation in face starting at vertex i
      Vertex_handle vj = f.vertex(f.ccw(i));    // ccw = counter clockwise .....
      IntCalcPrey(a, vi->info(), vj->info(), ptrSP, false);    // info returns the index
      IntCalcPrey(a, vj->info(), vi->info(), ptrSP, false);    // 2) symm=false because vij=vj for ASYNC_UPDATE

    }
}

void InteractionVoronoiF2FP(std::vector<particle> &a, params *ptrSP, std::vector<predator> &preds)
{
    // calculates local voronoi interactions
    typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
    typedef CGAL::Triangulation_vertex_base_with_info_2<int, K>         Vb;
    typedef CGAL::Triangulation_data_structure_2<Vb>                    Tds;
    typedef CGAL::Delaunay_triangulation_2<K, Tds>                      Delaunay;
    typedef Delaunay::Vertex_handle                                     Vertex_handle;
    typedef Delaunay::Edge_iterator                                     Edge_iterator;
    typedef Delaunay::Point                                             Point;
    typedef std::pair<Point, int>                                       PPoint;

    int N = a.size();
    int predId = -1;

    // create the delaunay triangulation network
    std::vector< std::pair<Point, int> > Vr;   // stores the point location and index
    std::vector< std::pair< std::vector<double>, int > > posId;
    Vr.reserve(4 * ( N + preds.size() ) );
    posId.reserve(4 * ( N + preds.size() ) );
    for(int i=0; i<N; i++){
        Point p1(a[i].x[0], a[i].x[1]);
        PPoint p2 = std::make_pair(p1, i);  // prey labeled with corresponding index
        Vr.push_back(p2);
        makePairAndPushBack(posId, a[i].x, i);
    }
    for (int i=0; i<preds.size(); i++){
        Point p1(preds[i].x[0], preds[i].x[1]);
        PPoint p2 = std::make_pair(p1, predId);     // predator labeled with negative index
        Vr.push_back(p2);
        makePairAndPushBack(posId, preds[i].x, predId);
        predId--;
    }
    // produce replicate prey/predator for periodic BC
    if (!ptrSP->BC){
        std::vector< std::pair< std::vector<double>, int > > newPosId = 
            GetCopies4PeriodicBC( posId, ptrSP->sizeL );
        for(int i=0; i<newPosId.size(); i++){
            Point p1(newPosId[i].first[0], newPosId[i].first[1]);
            PPoint p2 = std::make_pair(p1, newPosId[i].second);
            Vr.push_back(p2);
        }
    }

    Delaunay t;
    t.insert(Vr.begin(),Vr.end());

    //iterates over all finite edges and apply interaction 
    //  (infinite edges connect the infinite vertex with the vertices of the complex hull)
    for(Edge_iterator ei=t.finite_edges_begin(); ei!=t.finite_edges_end(); ei++){
        // Get a vertex from the edge, edge is stored as pair of the neighboring face and the vertex opposite to it
        Delaunay::Face& f = *(ei->first);
        int i = ei->second;
        Vertex_handle vi = f.vertex(f.cw(i));             // cw = clockwise rotation in face starting at vertex i
        Vertex_handle vj = f.vertex(f.ccw(i));            // ccw = counter clockwise .....
        if ((vi->info() >= 0) && (vj->info() >= 0)){      // both F
          IntCalcPrey(a, vi->info(), vj->info(), ptrSP, true);
        }
    }
    for (int i=0; i<a.size(); i++)
        for (int j=0; j<preds.size(); j++)
            IntCalcPred(a, i, preds[j], ptrSP);
}


void InteractionGlobal(std::vector<particle> &a, params *ptrSP)
{
    // Simple brute force algorithm for global interactions
    // checking all the N*(N-1)/2 combinations
    int i,j;
    int N = a.size();

    for(i=0;i<N;i++){
        for(j=i+1;j<N;j++){
            IntCalcPrey(a, i, j, ptrSP, false);
            IntCalcPrey(a, j, i, ptrSP, false);    // info returns the index
        }
    }
}

void InteractionPredGlobal(std::vector<particle> &a, params *ptrSP, std::vector<predator> &preds)
{
    // Simple brute force algorithm for global interactions /w predator only
    // checking all the N combinations
    int N = a.size();
    int i, j;

    for(i=0;i<N;i++){
        for(j=0;j<N;j++)
            IntCalcPred(a, i, preds[j], ptrSP);
    }
}

void IntCalcPrey(std::vector<particle> &a, int i, int j, params *ptrSP, bool symm)
{
    // Function updating social forces for a pair of interacting agents
    // Please Note that for local metric interaction the cutoff distance
    // is set to 1.2 x interaction range
    // only compute interaction if necessary (at start of burst)
    if ( a[i].bin_step != ptrSP->burst_steps )
        return;
    symm = false; // ASYNC_UPDATE
    // check if interaction already computed (only relevant for periodic BC)
    if (!ptrSP->BC){    // BC=0: periodic BC
        int there;
        there = where_val_in_vector<unsigned int>(a[i].NN,
                                                  static_cast<unsigned int>(j));
        if (there != a[i].NN.size())   // if interaction already computed
            return;
    }
    std::vector<double> r_ji(2);
    double u_ji[2];
    double dist_interaction;
    unsigned int c[3] = {0, 0, 0};
    double f0[2], f1[2], f2[2];
    f0[0] = f1[0] =  f2[0] = f0[1] = f1[1] =  f2[1] = 0;
    u_ji[0] = u_ji[1] = 0.0;
    // Calc relative distance vector and corresponding unit vector
    r_ji = CalcDistVec(a[i].x, a[j].x, ptrSP->BC, ptrSP->sizeL);  // vec i->j
    dist_interaction = vec_length(r_ji);
    if(dist_interaction > 0.0)
    {
        u_ji[0]=r_ji[0]/dist_interaction;
        u_ji[1]=r_ji[1]/dist_interaction;
    }
    double v_ji[2];     // needed for allignment and selective att,rep
    // v_ji[0] = a[j].v[0] - a[i].v[0];
    // v_ji[1] = a[j].v[1] - a[i].v[1];
    v_ji[0] = a[j].v[0];
    v_ji[1] = a[j].v[1];
    SFM_4Zone(ptrSP, u_ji, dist_interaction, v_ji,
              c, f0, f1, f2);
    if (c[0]){
        a[i].force_rep[0] -= f0[0];
        a[i].force_rep[1] -= f0[1];
        a[i].counter_rep++;
    }
    if (c[1]){
        a[i].force_alg[0] += f1[0];
        a[i].force_alg[1] += f1[1];
        a[i].counter_alg++;
    }
    if (c[2]){
        a[i].force_att[0] += f2[0];
        a[i].force_att[1] += f2[1];
        a[i].counter_att++;
    }
    // global, voronoi have symmetric interactions
    if (symm){
        if (c[0]){
            a[j].force_rep[0] += f0[0];
            a[j].force_rep[1] += f0[1];
            a[j].counter_rep++;
        }
        if (c[1]){
            a[j].force_alg[0] -= f1[0];
            a[j].force_alg[1] -= f1[1];
            a[j].counter_alg++;
        }
        if (c[2]){
            a[j].force_att[0] -= f2[0];
            a[j].force_att[1] -= f2[1];
            a[j].counter_att++;
        }
    }

    if (c[0] + c[1] + c[2] > 0){
        a[i].NN.push_back(j);
            if (symm)
                a[j].NN.push_back(i);
    }
}

void IntCalcPred(std::vector<particle> &a, int i, predator &pred, params *ptrSP)
{
    // Function updating social forces for predator and a single prey
    //////////////////
    // check if interaction already computed (only relevant for periodic BC)
    if (!ptrSP->BC){    // BC=0: periodic BC
        bool contains;
        contains = pred.NNset.count(static_cast<unsigned int>(i));
        if (contains) // interaction already computed
            return;
    }
    std::vector<double> r_ip(2);
    double u_ip[2];
    double ang;
    double dist_interaction;
    u_ip[0]=u_ip[1]=0.0;

    // Calc relative distance vector and corresponding unit vector
    r_ip = CalcDistVec(a[i].x, pred.x, ptrSP->BC, ptrSP->sizeL);
    ang = atan2(r_ip[1], r_ip[0]);
    dist_interaction = vec_length(r_ip);
    if(dist_interaction>0.0)
    {
        u_ip[0] = cos(ang);
        u_ip[1] = sin(ang);
    }
    double fstrength=0.0;
    // Calculate flee force
    // here fstrength represents the probability to react to predator as environmental cue
    // after flee_range the probability is below 1 to detect the predator:
    fstrength = 2 / ( 1 + dist_interaction / ptrSP->flee_range );
    double random = gsl_rng_uniform(ptrSP->r);
    // if(fstrength > 0.0){
    if(random < fstrength){
        a[i].counter_flee++;
        pred.NNset.insert(i);
        // ALTERNATIVE:
        a[i].force_flee[0] -= 1 * u_ip[0];
        a[i].force_flee[1] -= 1 * u_ip[1];
    }
}
