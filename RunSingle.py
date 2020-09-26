'''
    RunSingle
    Start a simulation of the burst-coast model by launching cpp-code 'swarmdyn'
    with the here defined parameters.
    Also visualizes the simulation by calling the script "AnimateRun.py".
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
import os
import time as pytime
import matplotlib
if __name__ == '__main__':
    matplotlib.use('TkAgg')
import SwarmDynByPy as swarmPy
import AnimateRun

def main():
    # Input Parameters
    #########################################
    animate = True
    dockerName = None # alternatively 'gcc_docker' see README.md for usage (Docker-alternative)

    fine = 18
    pred_time = 0     # 120 (voro)
    record_time = 20   # 20j

    mode = 'mulFisher'
    mode = 'sinFisher'
    mode = 'natPred'
    mode = 'burst_coast'
    dic = dict()
    dic = swarmPy.get_base_params(pred_time, record_time, mode=mode)
    dic['output_mode'] = 1

    # burst-coast or fishing
    SelectionLine = 1   # 0:LH, 1:RH, 2:SH
    swarmPy.selectionLineParameter(dic, SelectionLine)
    dic['pred_speed0'] = 20
    dic['rep_range'] = 0.001
    dic['alg_range'] = 0.001
    dic['att_range'] = 0.001
    dic['prob_social'] = 1

    # Generate and Run Command
    #########################################
    command = swarmPy.dic2swarmdyn_command(dic)
    print(command)
    t0 = pytime.time()
    dockerCall = ''
    if dockerName is not None:
        dockerCall = 'docker run --rm -v "$PWD":/usr/src/myapp -w /usr/src/myapp {} '.format(dockerName)
    os.system(dockerCall + 'make cl;')
    os.system(dockerCall + command)
    t1 = pytime.time()
    print(t1-t0)

    if animate:
        print('Animating')
        AnimateRun.main(dic['size'])
    return

if __name__ == '__main__':
    main()
