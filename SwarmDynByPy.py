'''
    SwarmDynByPy
    functions which define the parameters and provide the appropriate c-call
    to start the simulation in C++ via swarmdyn
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
'''
import numpy as np

def selectionLineParameter(dic, SL):
    '''
    set the parameter in which the selection lines differ:
    -probability to update social
    -burst rate
    -burst duration
    '''
    assert SL in [0, 1, 2], 'Unvalid parameter for selection line: {}'.format(SL)
    lineParas = {'burst_duration' : [0.089, 0.091, 0.097],
                 'burst_rate' : [5.7, 6.3, 4.7],
                 'prob_social' : [0.39, 0.47, 0.6],
                 'rep_range' : [3.92, 4.44, 4.71],
                 'att_range' : [15.84, 18.45, 19.5]}
    for k in lineParas.keys():
        dic[k] = lineParas[k][SL]
    # below: no orientation zone 
    dic['alg_range'] = lineParas['rep_range'][SL]


def get_base_params(pred_time, record_time, mode=None, trans_time=None):
    '''
    sets base parameters for different scenarios
    e.g. "burst-coast" scenario has a circular area where the agents 
        avaoid the wall

    '''
    if mode is None:
        mode = 'pred_prey'
    if trans_time is None:
        trans_time = 0
    assert mode in possible_modes, 'mode {} not known'.format(mode)
    params = dict()

    params['path'] = "./"
    params["fileID"] = 'xx'  # not passed to the code, only there
    params["out_h5"] = 1 # output-format 0:txt 1:hdf5
    params["trans_time"] = pred_time + trans_time    # time till output starts (transient)
    params["time"] = pred_time + record_time + trans_time    # total time of simulation
    params['dt'] = 0.001
    params['output'] = 0.03
    params["output_mode"] = 0 # 0:mean, 1:full, 2:no-output AND no-dummies (only pava_out)
    # IC    0:v-dis,x-dis-box, 1:v-order, x-dis, 2:v-order, x-dis-box
    #       3:v-order, x-dis-circ  4: v-milling, x-cricle  5:v-dis, x-circle
    #       6:v-dis, x-dis 99:from file
    params["IC"] = 3     # recommended:2 global, 3 voronoi
    # BC    -1:no, 0: periodic, 1:inelastic, 2:elastic, 3:x peri, y ela
    #       4:x peri, y inela 5 elastic circle, 6 inelastic circle 
    params["BC"] = -1

    params["N"] = 8
    params["Npred"] = 0

    # #################### F behavior
    params["Dphi"] = 0.02 # 0.02, only relevant for single-fishing agents: influences persistence length!
    params['beta'] = 2.51 # friction coefficient
    params["rep_range"] = 0.5
    params["alg_range"] = 16.2
    params["att_range"] = 30  # 30.0
    params['soc_strength'] = 110.82    # strength of social force (all social forces)
    params["burst_rate"] = 3.3

    params["flee_range"] = 7  # Range at which prey detects pred/fishing-agent with prob=1
    # #################### P behavior
    # N_confu = N_sense where prob_confu is 0.5
    params["N_confu"] = 4
    # pred_kill
    # 0: no-kill, 1: radial-kill, 2: probabilistic (r<r_kill),
    # 3: prob.+confusion, 4: prob.+confu.+selection
    params["pred_kill"] = 1
    params["pred_move"] = 0 # 0:random move, 1: follow closest
    params["pred_time"] = pred_time
    params["pred_speed0"] = 15
    params["kill_range"] = 5 # range in prey get killed by pred- or fishing-agent
    params["kill_rate"] = 1

    # #################### F reaction on P
    params['burst_duration'] = 0.121 # burst-time
    params['env_strength'] = params['soc_strength']  # strength of envirnomental force
    params['prob_social'] = 0.8   # prob to do social update
    params['alphaTurn'] = 17.62    # beta-turn: dphi/dt = alphaTurn * F * dt / v
    params['BC'] = 6               # boundary: 0 periodic, 5 elastic circular tank, 6 inelastic circular tank
    params['size'] = 30            # radius of circle
    if mode in ['natPred', 'natPredNoConfu', 'sinFisher', 'mulFisher']:
        params['N'] = 30
        params['BC'] = 0               # periodic BC
        params['size'] = 100
        params['Npred'] = 1
        if mode == 'sinFisher':
            params["pred_move"] = 0
        elif mode == 'natPred':
            params["pred_move"] = 1
            params["pred_kill"] = 3
        elif mode == 'natPredNoConfu':
            params["pred_move"] = 1
            params["pred_kill"] = 2
        elif mode == 'mulFisher':
            params['Npred'] =  int(params['size'] / 2)
    return params


def dic2swarmdyn_command(dic):
    '''
    transforms dictionary of parameters to calls to swarmdyn.cpp
    INPUT:
        dic dictionary
            keys = parameter names
            values = parameter values
    '''
    #consistency check:
    for key in dic.keys():
        if 'time' in key:
            assert np.any(dic[key] >= 0), '{} has negative value {}'.format(key, dic[key])
    command = './swarmdyn'
    command += ' -L %g' % dic['size']
    command += ' -N %d' % dic['N']
    command += ' -n %d' % dic['Npred']
    command += ' -D %g' % dic['Dphi']
    command += ' -b %g' % dic['beta']
    command += ' -f %g' % dic['flee_range']
    command += ' -e %g' % dic['pred_time']
    command += ' -S %g' % dic['pred_speed0']
    command += ' -h %g' % dic['rep_range']
    command += ' -a %g' % dic['alg_range']
    command += ' -r %g' % dic['att_range']
    command += ' -H %g' % dic['soc_strength']
    command += ' -A %g' % dic['burst_rate']
    command += ' -R %g' % dic['env_strength']
    command += ' -o %g' % dic['output']
    command += ' -J %d' % dic['out_h5']
    command += ' -d %g' % dic['dt']
    command += ' -t %g' % dic['time']
    command += ' -B %d' % dic['BC']
    command += ' -I %d' % dic['IC']
    command += ' -l %s' % dic['path']
    command += ' -m %d' % dic['output_mode']
    command += ' -c %d' % dic['N_confu']
    command += ' -T %g' % dic['trans_time']
    command += ' -x %d' % dic['pred_kill']
    command += ' -X %d' % dic['pred_move']
    command += ' -O %g' % dic['kill_rate']
    command += ' -G %g' % dic['kill_range']
    command += ' -E %s' % dic['fileID']
    command += ' -Q %g' % dic['prob_social']
    command += ' -u %g' % dic['burst_duration']
    command += ' -Y %g' % dic['alphaTurn']
    return command

possible_modes = ['burst_coast', 'sinFisher', 'mulFisher', 'natPred', 'natPredNoConfu', 'exploration']
