/*  AgentsOperation
    defines operations on Agents and is part of SwarmDynamics(swarmdyn.h, swarmdyn.cpp)
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
#include "agents_operation.h"

void GetCenterOfMass(std::vector<particle> &a, params *ptrSP, std::vector<int> &cluster, 
                     std::vector<double> &out, bool revise, unsigned int rev_time, double quantile)
{
    // Calculates center of mass

    // if no cluster is specified, get com of all particles
    if (cluster.size() == 0){
      cluster.resize(a.size(), 0);
      std::iota (std::begin(cluster), std::end(cluster), 0); //Fill with 0, 1,...N
    }
    int NN = cluster.size();
    int i, ii;
    if (!ptrSP->BC)
    {
        // Calculates periodic boundary condition COM, see Wikipedia
        double avcosx = 0;
        double avsinx = 0;
        double avcosy = 0;
        double avsiny = 0;
        double scaling = 2*M_PI/ptrSP->sizeL;
        for (i = 0; i < NN; i++)
        {
            ii = cluster[i];
            avcosx += cos(scaling*a[ii].x[0]);
            avsiny += sin(scaling*a[ii].x[1]);
            avcosy += cos(scaling*a[ii].x[1]);
            avsinx += sin(scaling*a[ii].x[0]);
        }
        avcosx /= NN;
        avsiny /= NN;
        avcosy /= NN;
        avsinx /= NN;
        out[0] = (atan2(-avsinx,-avcosx) + M_PI)/scaling;
        out[1] = (atan2(-avsiny,-avcosy) + M_PI)/scaling;
        // Adjusts to take out outliers
        if(revise){
            std::vector<double> dist(NN);
            std::vector<double> distsort(NN);
            std::vector<double> hv(2);
            int j, i, ii = 0;
            int k;
            int counter;
            // Adjusts "rev_time" times
            for (k = 0; k < rev_time; k++)
            {
                counter = 0; // counts number of prey which are not outliers
                for (i = 0; i < NN; i++)
                    ii = cluster[i];
                    hv = CalcDistVec(a[ii].x, out, ptrSP->BC, ptrSP->sizeL);
                    dist[i] = vec_length(hv);
                distsort = dist;
                std::sort(distsort.begin(), distsort.end());
                // Makes the array dist[] contain 1 if not an outlier, 0 if it is an outlier
                for (i = 0; i < NN; i++)
                {
                    for (j = 0; j < NN; j++)
                        if (dist[i] == distsort[j])
                            break;
                    // 1 if true, 0 if false
                    dist[i] = (j < quantile*NN);
                    counter += dist[i];
                }
                // Recalculate COM
                avcosx = 0;
                avsinx = 0;
                avcosy = 0;
                avsiny = 0;
                scaling = 2*M_PI/ptrSP->sizeL;
                for (i = 0; i < NN; i++)
                {
                    if (dist[i])
                    {
                        ii = cluster[i];
                        avcosx += cos(scaling*a[ii].x[0]);
                        avsiny += sin(scaling*a[ii].x[1]);
                        avcosy += cos(scaling*a[ii].x[1]);
                        avsinx += sin(scaling*a[ii].x[0]);
                    }
                }
                avcosx /= counter;
                avsiny /= counter;
                avcosy /= counter;
                avsinx /= counter;
                out[0] = (atan2(-avsinx,-avcosx) + M_PI)/scaling;
                out[1] = (atan2(-avsiny,-avcosy) + M_PI)/scaling;
            }
        }
    }
    else
    {
        // Regular COM
        unsigned int counter = 0;
        out[1] = 0;
        out[0] = 0;
        for (i = 0; i < NN; i++)
        {
            ii = cluster[i];
            counter++;
            out[0] += a[ii].x[0];
            out[1] += a[ii].x[1];
        }
        out[0]/= counter;
        out[1]/= counter;
    }
}


double AreaConvexHull(std::vector<particle> &a, std::vector<int> &nodes){
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_2 Point_2;
    typedef CGAL::Polygon_2<K> Polygon_2;
    int ii;
    Point_2 points[nodes.size()];
    Point_2 convhull[nodes.size()];
    for (unsigned int i=0; i<nodes.size(); i++){
            ii = nodes[i];
            points[i] = Point_2(a[ii].x[0], a[ii].x[1]);
    }
    Point_2 *ptr = CGAL::convex_hull_2( points, points+nodes.size(), convhull);

    // create a polygon and put some points in it
    Polygon_2 p;
    for (int i=0; i<ptr-convhull; i++)
      p.push_back(convhull[i]);

    double Area = p.area();
    return Area;
}


double get_elongation(std::vector<particle> &a, std::vector<double> &dir, std::vector<int> &nodes){
    std::vector<double> p_dir(2, 0); // defines perpendicular direction
    double len = vec_length(dir);
    std::vector<double> cdir = dir;
    cdir[0] /= len;
    cdir[1] /= len;
    p_dir[0] = cdir[1];
    p_dir[1] = -cdir[0];
    double min, max, p_min, p_max, dist, p_dist;
    int ii = 0;
    min = p_min = std::numeric_limits<double>::max();
    max = p_max = std::numeric_limits<double>::lowest();
    dist = p_dist = 0;
    for(int i=0; i<nodes.size(); i++){
        ii = nodes[i];
        dist = a[ii].x[0] * cdir[0] + a[ii].x[1] * cdir[1];
        p_dist = a[ii].x[0] * p_dir[0] + a[ii].x[1] * p_dir[1];
        min = fmin(min, dist);
        max = fmax(max, dist);
        p_min = fmin(p_min, p_dist);
        p_max = fmax(p_max, p_dist);
    }
    double elongation = (max - min) / (p_max - p_min);
    return fabs(elongation);
}

template <class O, class I>
std::vector<O> GetPredFrontPrey(std::vector<particle> &a, params *ptrSP, predator *pred, std::vector<I> &nodes)
{
    // returns vector of indices of prey in fron of pred 
    // only prey considered whose index is in "nodes" 
    unsigned int i;
    I ii;
    std::vector<double> r_pi(2);
    double front;   // if positive -> prey is in front 
    std::vector<O> results;
    for(i=0; i<nodes.size(); i++){
        ii = nodes[i];
        r_pi = CalcDistVec(pred->x, a[ii].x, ptrSP->BC, ptrSP->sizeL);
        front = pred->u[0] * r_pi[0] + pred->u[1] * r_pi[1];
        if (front > 0.0)
            results.push_back(ii);
    }
    return results;
}
template
std::vector<int> GetPredFrontPrey(std::vector<particle> &a, params *ptrSP, predator *pred, std::vector<int> &nodes);
template
std::vector<unsigned int> GetPredFrontPrey(std::vector<particle> &a, params *ptrSP, predator *pred, std::vector<unsigned int> &nodes);
template
std::vector<int> GetPredFrontPrey(std::vector<particle> &a, params *ptrSP, predator *pred, std::vector<unsigned int> &nodes);
template
std::vector<unsigned int> GetPredFrontPrey(std::vector<particle> &a, params *ptrSP, predator *pred, std::vector<int> &nodes);



void split_dead(std::vector<particle> &a, std::vector<particle> &d,
                std::vector<predator> &preds){
    for (int i=0; i<a.size(); i++){
        if (a[i].dead){
            d.push_back(a[i]);
            a.erase(a.begin() + i);
            i--;
        }
    }
}

void merge_dead(std::vector<particle> &a,
                                 std::vector<particle> &d){
    std::vector<particle> all(a.size() + d.size());
    for (int i=0; i<a.size(); i++)
        all[a[i].id] = a[i];
    for (int i=0; i<d.size(); i++)
        all[d[i].id] = d[i];
    d.resize(0);
    a = all;
}


void makePairAndPushBack(std::vector< std::pair< std::vector<double>, int > > &vecpair,
                         std::vector<double> &vec, int id){
    std::pair< std::vector<double>, int > newPair;
    newPair = std::make_pair( vec, id );
    vecpair.push_back( newPair );
}


std::vector< std::pair< std::vector<double>, int > > GetCopies4PeriodicBC(
        std::vector< std::pair< std::vector<double>, int > > &posId, double L)
{
    bool lower;
    bool left;
    double x, y;
    int id;
    std::vector<double> pos;
    std::vector<double> right {L, 0};
    std::vector<double> up {0, L};
    std::vector<double> newPos;
    std::vector< std::pair< std::vector<double>, int > > newPosId;
    newPosId.reserve(3 * posId.size()); // shifted copies of original pos (SAME ID)
    for (int i; i < posId.size(); i++){
        pos = posId[i].first;
        id = posId[i].second;
        if (pos[0] < L/2)
            left = true;
        if (pos[1] < L/2)
            lower = true;
        if (left && lower){
            newPos = vec_add(pos, up);
            makePairAndPushBack(newPosId, newPos, id);
            vec_add221(newPos, right);
            makePairAndPushBack(newPosId, newPos, id);
            newPos = vec_add(pos, right);
            makePairAndPushBack(newPosId, newPos, id);
        }
        else if (left && !lower){
            newPos = vec_sub(pos, up);
            makePairAndPushBack(newPosId, newPos, id);
            vec_add221(newPos, right);
            makePairAndPushBack(newPosId, newPos, id);
            newPos = vec_add(pos, right);
            makePairAndPushBack(newPosId, newPos, id);
        }
        else if (!left && !lower){
            newPos = vec_sub(pos, up);
            makePairAndPushBack(newPosId, newPos, id);
            vec_sub221(newPos, right);
            makePairAndPushBack(newPosId, newPos, id);
            newPos = vec_sub(pos, right);
            makePairAndPushBack(newPosId, newPos, id);
        }
        else{ // (!left && lower)
            newPos = vec_add(pos, up);
            makePairAndPushBack(newPosId, newPos, id);
            vec_sub221(newPos, right);
            makePairAndPushBack(newPosId, newPos, id);
            newPos = vec_sub(pos, right);
            makePairAndPushBack(newPosId, newPos, id);
        }
    }
    return newPosId;
}
