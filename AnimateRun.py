'''
    AnimateRun
    How this program works:

    This program creates a visualization of a simulation. There are three different modes to visualize,
    one showing the global movement, one showing movement relative to center of mass, and one showing
    movement relative to a moving center of mass. You can adjust size, delaytime, and taillength of the
    visualization.

    In 'normal' mode, the visualization is simply displayed
    In 'pictures' mode, the visualization is displayed, and the images are saved in the current directory
    In 'movie' mode, the visualization is not displayed, and the images are made into an .avi file in the
    current directory

    v0.1, 13.5.2020

    (C) 2020 Pascal Klamser, Pawel Romanczuk, Ishan Levy

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
from time import sleep
import configparser
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import h5py
from pathlib import Path

def h5ReadDset(h5file, dname):
    data = h5file.get(dname)
    np_data = np.array(data)
    print(np_data)
    if gname in list(h5file.keys()):
        data2 = h5file.get(gname + '/datadouble')
        np_data2 = np.array(data2)
        print(np_data2)

def init_head_tails(data, size, colors, ax, cmap=None, scatter=None):
    '''
    INPUT:
        data.shape(Time, N, N_coord)
            data containing for "Time" time points the position for "N"
            agents.
        size float
            marker size of the plot
        colors.shape(N)
            array or list containing for each particle a colorcode (float) 
        ax matplotlib.axes.AxesSubplot object
            Subplot where the data is plotted
        cmap matplotlib.Colormap
            e.g.: cmap = plt.get_cmap('Reds')
        scatter boolean
            if true scatter is used instead of plot
    OUTPUT:
        heads
        tails
    '''
    if scatter is None:
        scatter = True
    if len(data.shape) == 2:
        x_init = data[0, 0]
        y_init = data[0, 1]
    if len(data.shape) == 3:
        x_init = data[0, :, 0]
        y_init = data[0, :, 1]
    if scatter:
        heads = ax.scatter(x_init, y_init, marker='o',
                           c=colors, cmap=cmap, s=size**2)
        tails = ax.scatter(x_init, y_init, marker='o',
                           c=colors, cmap=cmap, s=(size/2.)**2,
                           alpha=0.3)
    else:
        heads = ax.plot(x_init, y_init, marker='o',
                        color=colors, ms=size, linestyle='none')[0]
        tails = ax.plot(x_init, y_init, marker='o',
                        color=colors, ms=size/2., linestyle='none',
                        alpha=0.3)[0]
    return heads, tails


def update_head_tails(head, tail, data, s, tail_length, scatter=None):
    '''
    INPUT:
        head 1st output of init_head_tails
        tail 1st output of init_head_tails
        data.shape(Time, N, N_coord)
            data containing for "Time" time points the position for "N"
            agents.
        s int
            timestep s /in [0, Time-1]
        tail_length int
            length of tail
        scatter boolean
            if true scatter is used instead of plot
    '''
    if scatter is None:
        scatter = True
    endtail = s - tail_length
    if endtail < 0:
        endtail = 0
    if scatter:
        if len(data.shape) == 2:
            dat = data[endtail:s+1, :2]
        if len(data.shape) == 3:
            dat = data[endtail:s+1, :, :2]
        head.set_offsets(dat[-1])
        tail.set_offsets(dat[:-1].reshape(-1, 2))
    else:
        if len(data.shape) == 2:
            x_dat = data[endtail:s+1, 0]
            y_dat = data[endtail:s+1, 1]
        if len(data.shape) == 3:
            x_dat = data[endtail:s+1, :, 0]
            y_dat = data[endtail:s+1, :, 1]
        head.set_data(x_dat[-1], y_dat[-1])
        tail.set_data(x_dat[:-1], y_dat[:-1])

def main(simuMode=None):
    if simuMode is None:
        simuMode = 'natPred'
        simuMode = 'burst_coast'
    center = 0    # 0: TANK->center=(0, 0)+shows boundary 1: Moving Center 2: FISHING->fixed frame for periodic BC
    if simuMode == 'burst_coast':
        center = 0
    else:
        center = 2
    centersize = 100;    # defines
    tail_length = 5

    ShowBoundaryCircle = False   # shows circle of size = L around center (BC=5)
    if center == 0:
        ShowBoundaryCircle = True
    else:
        ShowBoundaryCircle = False 
    minframe = 0
    maxframe = -1 #-1 is all the frames
    mode = "pictures"
    mode = "movie"
    mode = "normal"
    moviename = "Animation"
    fps = 30 # 13
    timer = 1/(fps + fps)

    ####################################

    fileID = 'xx'
    folder = Path.cwd()

    ####### check which input shall be used ######
    try:
        pardat = np.genfromtxt(str(folder / ('part_'+fileID+'.dat')))
        inmode = "txt"  # if data is stored in hdf5(=h5) file or txt file
    except:
        inmode = "h5"  # if data is stored in hdf5(=h5) file or txt file

    ######## get parameters ##############
    config = configparser.RawConfigParser()
    config.read(str(folder /  'parameters.dat'))
    N = config.getint('Parameters', 'particles')
    L = config.getfloat('Parameters', 'size')
    BC = config.getint('Parameters', 'BC')
    Npred = config.getint('Parameters', 'Npred')
    sim_time = config.getfloat('Parameters', 'time')
    pred_time = config.getfloat('Parameters', 'pred_time')
    trans_time = config.getfloat('Parameters', 'trans_time')
    output = config.getfloat('Parameters', 'output')
    t0p = int((pred_time-trans_time)/output)
    if t0p < 0:
        t0p = 0
    totaltime = (sim_time-trans_time)/output
    print(('N={},L={},BC={},Npred={}'.format(N, L, BC, Npred)))
    if center == 2:
        centersize = L
    ######## load data ###################
    if inmode == "txt":
        pardat = np.genfromtxt(str(folder / ('part_'+fileID+'.dat')))
        frames_total = int(np.shape(pardat)[0]/N)
        if (maxframe == -1):
            maxframe = frames_total
        variables = np.shape(pardat)[1]
        print(('total frames: {}'.format(frames_total)))
        print('totaltime: ', totaltime)
        fdat = np.reshape(pardat, (frames_total, -1, variables))   # shape=time, N, variable
    elif inmode == "h5":
        f = h5py.File(str(folder / ('out_' + fileID + '.h5')), 'r')   # open existing file, must exist
        fdat = np.array(f.get('part'))
        frames_total = fdat.shape[0]
        if (maxframe == -1):
            maxframe = frames_total

    ######## load pred data ##############
    if(Npred != 0):
        if inmode == "txt":
            pdat = np.genfromtxt(str(folder / ('pred_'+fileID+'.dat')))
        elif inmode == "h5":
            pdat = np.array(f.get('pred'))
        if len(pdat.shape) == 2:
            pdat = np.expand_dims(pdat, axis=1)
            pDdat = np.expand_dims(pDdat, axis=1)

    ######## define colors for different agents ##############
    colors = list(range(N))
    colors = N * ['k']
    cmap = None
    ######## animate #####################
    fig, ax = plt.subplots(1, 1)
    ax.set_aspect('equal')
    # ax.axis('off')
    for spine in plt.gca().spines.values():
        spine.set_visible(False)
        plt.tick_params(top='on', bottom='on', left='on', right='on', labelleft='off', labelbottom='off')
    if (not center):
        ax.set_xlim(-L, L)
        ax.set_ylim(-L, L)
    else:
        ax.set_xlim(0, centersize)
        ax.set_ylim(0, centersize)
    if ShowBoundaryCircle:
        circle = mpatches.Circle((0, 0), L, color='b', fill=False)
        ax.add_artist(circle)
    plt.tight_layout()
    plt.show(False)

    plt.draw()

    preypos = np.array(fdat[:, :, 0:2])
    preypos *= 2*np.pi/L
    if(center == 1):
        preyposcos = np.cos(preypos)
        preypossin = np.sin(preypos)
        avposcos = np.average(preyposcos, axis=1)
        avpossin = np.average(preypossin, axis=1)
        avpos = np.arctan2(-avpossin, -avposcos) + np.pi
        avpos *= L/(2*np.pi)
        avpos -= centersize/2
        for i in range(len(fdat[0, :, 0])):
            fdat[:, i, 0:2] = ((fdat[:, i, 0:2]-avpos) + 0.5*L) % L - 0.5*L
        if(Npred != 0):
            for i in range(len(pdat[0, :, 0])):
                pdat[:, i, 0:2] = ((pdat[:, i, 0:2] - avpos[t0p:]) + 0.5*L) % L - 0.5*L

    elif(center == -1):
        coarsepos = np.zeros(np.shape(fdat[:, :, 0:2]))
        com = np.zeros((len(fdat[:]), 2))
        com[0] = np.mean(fdat[0, :, 0:2], axis=0)
        for t in range(1, len(preypos[:, 0, 0])):
                diffX = preypos[t, :, 0]-preypos[t-1, :, 0]
                diffY = preypos[t, :, 1]-preypos[t-1, :, 1]
                # check for jump in X and update coarse pos variable
                condXdecr = np.where(diffX > 0.5*L)
                condXincr = np.where(diffX < -0.5*L)
                coarsepos[t, condXdecr, 0] -= L
                coarsepos[t, condXincr, 0] += L
                # check for jump in Y and update coarse pos variable
                condYdecr = np.where(diffY > 0.5*L)
                condYincr = np.where(diffY < -0.5*L)
                coarsepos[t, condYdecr, 1] -= L
                coarsepos[t, condYincr, 1] += L
                com[t] = np.mean(fdat[t, :, 0:2], axis=0)

    radii = np.array([1/4., 10, 25])
    pixels = ax.transData.transform([(0, 1), (1, 0)])-ax.transData.transform((0, 0))
    xpix = pixels[1].max()
    ypix = pixels[0].max()
    size = radii*xpix/77.*100

    points, pointstail = init_head_tails(fdat, size[0],
                                         colors, ax, cmap=cmap)
    if (t0p < totaltime):
        Ppoints, Ppointstail = init_head_tails(pdat, size[0]*2,
                                               'r', ax, scatter=False )


    for s in range(minframe, maxframe):
        if (center == -1):
            ax.set_xlim(com[s, 0]-0.5*centersize, com[s, 0]+0.5*centersize)
            ax.set_ylim(com[s, 1]-0.5*centersize, com[s, 1]+0.5*centersize)
        if (s >= t0p and Npred > 0):
            update_head_tails(Ppoints, Ppointstail, pdat, s-t0p, tail_length, scatter=False)
        update_head_tails(points, pointstail, fdat, s, tail_length)
        ax.set_title('step=%d' % (s))
        if(mode != "movie"):
            fig.canvas.draw()
            sleep(timer)
        if(mode != "normal"):
            fig.savefig('f%06d.jpg' % s, dpi=300)

    if (mode == "movie"):
        os.system("mencoder mf://f*.jpg -mf fps={}:type=jpg  -vf scale=-10:-1 -ovc x264 -x264encopts  bitrate=2400 -o {}.avi".format(fps, moviename))
        os.system("rm f*.jpg")

if __name__ == '__main__':
    main()
