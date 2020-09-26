#!/usr/bin/env python
import numpy as np
import matplotlib
if __name__ == '__main__':
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import animation
import matplotlib.patches as mpatches
import h5py
from pathlib import Path
from functools import partial
import pdb


class datCollector:
    def __init__(self, dat):
        if dat is None:
            self.pos = None
            self.vel = None
        else:
            dat = dat.T
            self.pos = dat[:2].T
            self.vel = dat[2:4].T
            if len(dat) > 4:
                self.fitness = dat[4].T
                self.force = dat[5:].T
            dead = np.where( np.sum(self.pos, axis=-1) == 0 )
            self.pos[dead] = np.nan
            self.vel[dead] = np.nan


class headAndTail:
    def __init__(self, data, size, colors, ax, tail_length,
                 cmap=None, scatter=None, marker=None, delay=None):
        '''
        INPUT:
            data.shape(Time, N, N_coord)
                data containing for "Time" time points the position for "N"
                agents.
            size float
                radius of marker size
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
        if delay is None:
            delay = 0
        self.pos = data.pos
        self.tail_length = tail_length
        self.scatter = scatter
        self.delay = delay
        self.ax = ax
        self.sizeRaw = size
        if marker is None:
            marker = 'o'
        if self.scatter is None:
            self.scatter = True
        if self.scatter:
            self.createFunc = partial(ax.scatter, marker=marker, c=colors, cmap=cmap)
        else:
            self.createFunc = partial(ax.plot, marker=marker, c=colors, linestyle='none')


    def sizeAdjust(self):
        trans = self.ax.transData
        pixels = trans.transform([(0, 1), (1, 0)]) - trans.transform((0, 0))
        xpix = pixels[1].max()
        self.size = self.sizeRaw * xpix / 77. * 100
        if self.scatter:
            self.size = self.size ** 2


    def create(self):
        posInit = self.pos[0]   # use Time=1 for initialization
        if self.scatter:
            self.heads = self.createFunc(posInit.T[0], posInit.T[1], s=self.size)
            self.tails = self.createFunc(posInit.T[0], posInit.T[1],
                                         s=self.size/4, alpha=0.5)
        else:
            self.heads = self.createFunc(posInit.T[0], posInit.T[1], ms=self.size)[0]
            self.tails = self.createFunc(posInit.T[0], posInit.T[1],
                                         ms=self.size/2, alpha=0.5)[0]


    def update(self, s):
        self.sizeAdjust()
        s = s - self.delay
        if s == 0:
            self.create()
        elif s > 0:
            endtail = s - self.tail_length
            if endtail < 0:
                endtail = 0
            if self.scatter:
                tail = self.pos[endtail:s].reshape(-1, 2)
                self.heads.set_offsets(self.pos[s])
                self.tails.set_offsets(tail)
                self.heads.set_sizes(len(self.pos[s]) * [self.size])
                self.tails.set_sizes(len(tail) * [self.size / 4])
            else:
                self.heads.set_data(self.pos[s].T[0], self.pos[s].T[1])
                self.tails.set_data(self.pos[endtail:s].T[0],
                                    self.pos[endtail:s].T[1])
                self.heads.set_markersize(self.size)
                self.tails.set_markersize(self.size / 2)
        if s >= 0:
            return self.heads, self.tails


def pavas2colors(pavas):
    if np.std(pavas) > 1e-5:
        colors = np.squeeze(pavas)
        colors -= colors.min()
        colors /= colors.max()
    else:
        colors = 'k'
    return colors


class taskCollector:
    '''
    collects tasks (=objects with update(s) method returning artist objects)
    '''
    def __init__(self):
        self.taskList = []

    def append(self, task):
        self.taskList.append(task)

    def update(self, s):
        '''
        returns all output summarized
        '''
        artistObjs = ()
        for task in self.taskList:
            artist = task.update(s)
            if type(artist) != tuple:
                artist = tuple([artist])
            # print(type(artistObjs), type(artist), task.__class__)
            artistObjs += artist 
        return artistObjs


def UpdateViaAnimation(fig, tasks, tmin, tmax, fps=None, dpi=None,
                       mode=None, name=None, repeat=None):
    fps = setDefault(fps, 15)
    dpi = setDefault(dpi, 300)
    mode = setDefault(mode, 'normal')
    name = setDefault(name, 'Animation')
    name = setDefault(repeat, True)
    interval = 1000*(1/fps)
    anim = animation.FuncAnimation(fig, tasks.update, interval=interval,
                                   frames=range(tmin-1, tmax), repeat=repeat)
    if mode == 'movie':
        anim.save(name + '.mp4', writer='ffmpeg', dpi=dpi)
    elif mode == 'gif':
        anim.save(name + '.gif', writer='imagemagick', dpi=dpi)
    else:
        plt.show()


def main(size):
    if size is None:
        size = 25
    mode = 'normal' # 'normal', 'pictures', 'movie', 'gif'
    # # Definitions
    tail_length = 3
    fps = 30
    dpi = 200
    sizePrey = 1/4
    sizePred = sizePrey * 2
    name = 'Animation'
    cmap = plt.get_cmap('coolwarm') # alternatives 'bwr', 'Reds'
    folder = Path.cwd()

    # # Load data 
    f_h5 = folder / 'out_xx.h5'
    if f_h5.exists:
        with h5py.File(f_h5) as fh5:
            items = []
            fh5.visit(items.append)
            preys = datCollector( np.array(fh5['/part']) )
            if 'pred' in items:
                preds = datCollector( np.array(fh5['/pred']) )
    else:
        print('{} does not EXIST'.format(f_h5))
    pava_in_name = folder / 'pava_in_xx.in'
    colors = 'k'
    if pava_in_name.exists():
        pavas = np.loadtxt(str(pava_in_name))
        colors = pavas2colors(pavas)
    # get info from files
    time, N, _ = preys.pos.shape
    if 'pred' in items:
        timePred, Npred, _ = preds.pos.shape
        t0pred = time - timePred # predator can appear later

    # # Animation 
    f, ax = plt.subplots(1)
    # ax.axis('off')
    for spine in plt.gca().spines.values():
        spine.set_visible(False)
    plt.tick_params(top='on', bottom='on', left='on', right='on')   # , labelleft='off', labelbottom='off')
    ax.set_aspect('equal')
    ax.set_ylim(0, size)
    ax.set_xlim(0, size)
    if size < 100:  # burst_coast scenario -> draw Circle
        ax.set_ylim(-size, size)
        ax.set_xlim(-size, size)
        circle = mpatches.Circle((0, 0), size, color='b', fill=False)
        ax.add_artist(circle)
    # Collect update-tasks
    tasks = taskCollector()
    tasks.append( headAndTail(preys, sizePrey, colors, ax,
                              tail_length, cmap=cmap) )
    if 'pred' in items:
        tasks.append( headAndTail(preds, sizePred, 'r', ax, tail_length,
                                  scatter=False, delay=t0pred) )
    # animation
    interval = 1000*(1/fps)
    anim = animation.FuncAnimation(f, tasks.update, interval=interval,
                                   frames=range(0-1, time), repeat=False)
    if mode == 'movie':
        anim.save(name + '.mp4', writer='ffmpeg', dpi=dpi)
    elif mode == 'gif':
        anim.save(name + '.gif', writer='imagemagick', dpi=dpi)
    else:
        plt.show()


if __name__ == '__main__':
    main(None)
