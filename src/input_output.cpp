/*  InputOutput
    defines basic in- and output operations for 
    SwarmDynamics(swarmdyn.h, swarmdyn.cpp)
    v0.1, 13.5.2020

    (C) 2020 Pascal Klamser, Pawel Rormanczuk

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
#include "input_output.h"


void WritePosVelPlus(std::vector<particle> &a, std::vector<double> plus,
                     params* ptrSP, std::string name, bool append){
    std::ios_base::openmode mode;
    if (name == "")
        name = "posvelplus";
    if (append)
        mode = std::ios_base::app;
    else
        mode = std::ios_base::out;
    std::ofstream outFile((ptrSP->location + name + ".dat").c_str(),
                          mode);
    unsigned int id = 0;
    for(unsigned int i=0; i < a.size(); i++){
        while (id < a[i].id){
            outFile << 0 << " " << 0 << " "
                    << 0 << " " << 0 << " "
                    << 0 << std::endl;
            id++;
        }
        outFile << a[i].x[0] << " " << a[i].x[1] << " "
                << a[i].v[0] << " " << a[i].v[1] << " " 
                << plus[i] << std::endl;
        id++;
    }
    // to assure the output has always same length
    for(unsigned int i=id; i < ptrSP->N; i++){
        outFile << 0 << " " << 0 << " "
                << 0 << " " << 0 << " "
                << 0 << std::endl;
    }

    outFile << std::endl;
}

void WritePosVel(std::vector<particle> &a, params* ptrSP,
                 std::string name, bool append){
    std::ios_base::openmode mode;
    if (name == "")
        name = "posvel";
    if (append)
        mode = std::ios_base::app;
    else
        mode = std::ios_base::out;
    std::ofstream outFile((ptrSP->location + name + ".dat").c_str(),
                          mode);
    for(unsigned int i=0; i < a.size(); i++)
        outFile << a[i].x[0] << " " << a[i].x[1] << " "
                << a[i].v[0] << " " << a[i].v[1] << std::endl;
    outFile << std::endl;
}


void WritePosVelDead(std::vector<particle> &a, params &SP, std::string name, predator &pred){
    if (name == "")
        name = "pos_vel_dead_";
    std::ofstream outFile((SP.location + name + SP.fileID
                           + ".dat").c_str(), std::ios::app);
    unsigned int id = 0;
    for(int i=0; i<a.size(); i++){
        while (id < a[i].id){
            outFile << 0 << " " << 0 << " "
                    << 0 << " " << 0 << " "
                    << 0 << std::endl;
            id++;
        }
        outFile << a[i].x[0] << " " << a[i].x[1] << " "
                << a[i].v[0] << " " << a[i].v[1] << " "
                << int(a[i].dead) << std::endl;
        id++;
    }
    // to ensure that there are always same Nr of rows
    for(int i=id; i<SP.N; i++){
        outFile << 0 << " " << 0 << " "
                << 0 << " " << 0 << " "
                << 0 << std::endl;
    }
    // write predator values:
    outFile << pred.x[0] << " " << pred.x[1] << " "
            << pred.v[0] << " " << pred.v[1] << " "
            << 0 << std::endl;
    outFile << std::endl;
}

const char* getCmdOption(char ** begin, char ** end, const std::string & option){
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end){
        return *itr;
    }
    return "0";
}

void ParseParameters(int argc, char **argv, params *SysParams)
{
    // Function for parsing the parameters from command line
    // free short option arguments: y, Z
    SysParams->sizeL = atof(getCmdOption(argv, argv+argc, "-L"));
    SysParams->location = getCmdOption(argv, argv+argc, "-l");
    SysParams->N = atoi(getCmdOption(argv, argv+argc, "-N"));
    SysParams->Npred = atoi(getCmdOption(argv, argv+argc, "-n"));
    SysParams->flee_range = atof(getCmdOption(argv, argv+argc, "-f"));
    SysParams->pred_time = atof(getCmdOption(argv, argv+argc, "-e"));
    SysParams->pred_speed0 = atof(getCmdOption(argv, argv+argc, "-S"));
    SysParams->Dphi = atof(getCmdOption(argv, argv+argc, "-D"));
    SysParams->beta = atof(getCmdOption(argv, argv+argc, "-b"));
    SysParams->dt = atof(getCmdOption(argv, argv+argc, "-d"));
    SysParams->sim_time = atof(getCmdOption(argv, argv+argc, "-t"));
    SysParams->trans_time = atof(getCmdOption(argv, argv+argc, "-T"));
    SysParams->rep_range = atof(getCmdOption(argv, argv+argc, "-h"));
    SysParams->alg_range = atof(getCmdOption(argv, argv+argc, "-a"));
    SysParams->att_range = atof(getCmdOption(argv, argv+argc, "-r"));
    SysParams->soc_strength = atof(getCmdOption(argv, argv+argc, "-H"));
    SysParams->burst_rate = atof(getCmdOption(argv, argv+argc, "-A"));
    SysParams->env_strength = atof(getCmdOption(argv, argv+argc, "-R"));
    SysParams->output = atof(getCmdOption(argv, argv+argc, "-o"));
    SysParams->BC = atoi(getCmdOption(argv, argv+argc, "-B"));
    SysParams->IC = atoi(getCmdOption(argv, argv+argc, "-I"));
    SysParams->N_confu = atoi(getCmdOption(argv, argv+argc, "-c"));
    SysParams->output_mode = atoi(getCmdOption(argv, argv+argc, "-m"));
    SysParams->kill_range = atof(getCmdOption(argv, argv+argc, "-G"));
    SysParams->pred_kill = atoi(getCmdOption(argv, argv+argc, "-x"));
    SysParams->pred_move = atoi(getCmdOption(argv, argv+argc, "-X"));
    SysParams->fileID = getCmdOption(argv, argv+argc, "-E");
    SysParams->out_h5 = atoi(getCmdOption(argv, argv+argc, "-J"));
    SysParams->prob_social = atof(getCmdOption(argv, argv+argc, "-Q"));
    SysParams->burst_duration = atof(getCmdOption(argv, argv+argc, "-u"));
    SysParams->alphaTurn = atof(getCmdOption(argv, argv+argc, "-Y"));
    SysParams->kill_rate = atof(getCmdOption(argv, argv+argc, "-O"));
    // Sets auxillary variables
    if(SysParams->output>SysParams->dt)
        SysParams->step_output=(int) (SysParams->output/SysParams->dt);
    else
        SysParams->step_output=1;
    SysParams->sim_steps =(int) (SysParams->sim_time/SysParams->dt);
}


void OutputParameters(params SysParams)
{
    FILE *fp;
    // Write parameters to file
    fp=fopen((SysParams.location + "parameters.dat").c_str(), "w");
    fprintf(fp,"[Parameters]\n");
    fprintf(fp,"size:               \t%g\n",SysParams.sizeL);
    fprintf(fp,"particles:          \t%d\n",SysParams.N);
    fprintf(fp,"dt:                 \t%g\n",SysParams.dt);
    fprintf(fp,"time:               \t%g\n",SysParams.sim_time);
    fprintf(fp,"trans_time:         \t%g\n",SysParams.trans_time);
    fprintf(fp,"Dp:                 \t%g\n",SysParams.Dphi);
    fprintf(fp,"rep_range:          \t%g\n",SysParams.rep_range);
    fprintf(fp,"alg_range:          \t%g\n",SysParams.alg_range);
    fprintf(fp,"att_range:          \t%g\n",SysParams.att_range);
    fprintf(fp,"soc_strength:       \t%g\n",SysParams.soc_strength);
    fprintf(fp,"prob_social:       \t%g\n",SysParams.prob_social);
    fprintf(fp,"burst_rate:         \t%g\n",SysParams.burst_rate);
    fprintf(fp,"env_strength:       \t%g\n",SysParams.env_strength);
    fprintf(fp,"beta:               \t%g\n",SysParams.beta);

    fprintf(fp,"Npred:              \t%d\n",SysParams.Npred);
    fprintf(fp,"flee_range:         \t%g\n",SysParams.flee_range);
    fprintf(fp,"pred_time:          \t%g\n",SysParams.pred_time);
    fprintf(fp,"pred_speed0:        \t%g\n",SysParams.pred_speed0);

    fprintf(fp,"output_mode:        \t%d\n",SysParams.output_mode);
    fprintf(fp,"N_confu:            \t%d\n",SysParams.N_confu);
    fprintf(fp,"output:             \t%g\n",SysParams.output);
    fprintf(fp,"BC:                 \t%d\n",SysParams.BC);
    fprintf(fp,"IC:                 \t%d\n",SysParams.IC);
    fprintf(fp,"predator_circle_rad:\t%g\n",SysParams.kill_range);
    fprintf(fp,"pred_kill:          \t%d\n",SysParams.pred_kill);
    fprintf(fp,"pred_move:          \t%d\n",SysParams.pred_move);
    fprintf(fp,"out_h5:             \t%d\n",SysParams.out_h5);
    fprintf(fp,"burst_duration:     \t%g\n",SysParams.burst_duration);
    fprintf(fp,"alphaTurn:          \t%g\n",SysParams.alphaTurn);
    fprintf(fp,"kill_rate:          \t%g\n",SysParams.kill_rate);

    fclose(fp);
}

void LoadCoordinatesCPP(params * ptrSP, std::string name, std::vector<particle> &a)
{
    if (name == "")
        name = "foo";
    std::ifstream inputFile((ptrSP->location + name + ".in").c_str());
    unsigned int ii = 0;
    double phi = 0;
    if (inputFile){
        double value;
        while ( inputFile >> value ) {
            a[ii].x[0] = value;
            if ( inputFile >> value ) a[ii].x[1] = value;;
            if ( inputFile >> value ) a[ii].v[0] = value;;
            if ( inputFile >> value ) a[ii].v[1] = value;;
            phi = atan2(a[ii].v[1], a[ii].v[0]);
            a[ii].phi = phi;
            a[ii].u[0] = cos(phi);
            a[ii].u[0] = sin(phi);
            a[ii].vproj = sqrt(a[ii].v[0]*a[ii].v[0] + a[ii].v[1]*a[ii].v[1]);
            ii++;
        }
    }
    std::cout<< ii << "CoordsLoaded ";
}

void LoadCoordinates(params * ptrSP, const char *fn, std::vector<particle> &a, int N, double sizeL)
{
    // Function for loading coordinates as initial conditions
    // fn - input file name, containing 4 columns: X,Y,VX,VY (same format as final_coords.dat)

    FILE *input_ptr;

    double init_mv[2];
    char line[128];
    int i;
    int line_count=0;
    int idx_max=0;
    input_ptr = fopen((ptrSP->location + fn).c_str(), "r");
    while(fgets(line, 128, input_ptr)!=NULL)
    {
        line_count++;
        }
    printf("line count=%d\n", line_count);
    rewind(input_ptr);
    fclose(input_ptr);

    float *init_coord = new float(4*line_count*sizeof(float));;
    input_ptr = fopen((ptrSP->location + fn).c_str(), "r");
    i=0;
    init_mv[0]=init_mv[1]=0.0;  // initial mean velocity
    while(fgets(line, 128, input_ptr)!=NULL)
    {
        sscanf(line,"%g\t%g\t%g\t%g",&init_coord[4*i+0],&init_coord[4*i+1],&init_coord[4*i+2],&init_coord[4*i+3]);
        init_mv[0]+=init_coord[4*i+2];
        init_mv[1]+=init_coord[4*i+3];
        i++;
        }
    rewind(input_ptr);
    fclose(input_ptr);

    init_mv[0]/=line_count;
    init_mv[1]/=line_count;
    printf("initial vel. mvx=%.3f, mvy=%.3f\n",init_mv[0],init_mv[1]);

    if(N<line_count)
        idx_max=N;
    else{
        idx_max=line_count;
        for(i=line_count;i<N;i++)
        {
            a[i].x[0] = double(init_coord[4*i-line_count+0]) + 0.01;
            a[i].x[1] = double(init_coord[4*i-line_count+1]) + 0.01;
            a[i].v[0] = init_mv[0];
            a[i].v[1] = init_mv[1];
        }
    }
    for(i=0;i<idx_max;i++)
    {
        a[i].x[0]=(double) init_coord[4*i+0];
        a[i].x[1]=(double) init_coord[4*i+1];
        a[i].v[0]=(double) init_coord[4*i+2];
        a[i].v[1]=(double) init_coord[4*i+3];
    }
    delete init_coord;
}

void LoadVector(params * ptrSP, std::string name, std::vector<double> &vec)
{
    vec.resize(0);
    if (name == "")
        name = "foo";
    std::ifstream inputFile((ptrSP->location + name + ".dat").c_str());
    if (inputFile){
        double value;
        while ( inputFile >> value )
            vec.push_back(value);
    }
}

template<class T>
void WriteVector(std::string file, std::vector<T> &vec, bool append){
    std::ios_base::openmode openmode;
    if(append)
        openmode = std::ios_base::app;
    else
        openmode = std::ios_base::trunc;
    std::ofstream outFile(file.c_str(), openmode);
    for (int i=0; i<vec.size(); i++)
        outFile << vec[i] << " ";
    outFile << std::endl;
}
template
void WriteVector(std::string file, std::vector<double> &vec, bool append);
template
void WriteVector(std::string file, std::vector<int> &vec, bool append);
template
void WriteVector(std::string file, std::vector<unsigned int> &vec, bool append);

template<class T>
void WriteVector2d(std::string file, std::vector< std::vector<T> > &vec, bool append){
    std::ios_base::openmode openmode;
    if(append)
        openmode = std::ios_base::app;
    else
        openmode = std::ios_base::trunc;
    int n = vec.size();
    int m = vec[0].size();
    if (n > 0){
        std::ofstream outFile(file.c_str(), openmode);
        for (int i=0; i<vec.size(); i++){
            for (int j=0; j<m; j++)
                outFile << vec[i][j] << " ";
            outFile << std::endl;
        }
        outFile << std::endl;
    }
}
template
void WriteVector2d(std::string file, std::vector< std::vector<double> > &vec, bool append);
template
void WriteVector2d(std::string file, std::vector< std::vector<int> > &vec, bool append);
template
void WriteVector2d(std::string file, std::vector< std::vector<unsigned int> > &vec, bool append);
