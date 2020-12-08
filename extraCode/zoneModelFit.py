import numpy as np
import scipy
from scipy import optimize
from numba import njit, jit
import h5py
import csv
import os
import glob
import pickle
import scipy.stats as stats
from functools import partial
import pandas as pd
from scipy import stats
from scipy.spatial import Delaunay
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

from TrackScrap import general as gen

from animateSwarm import AnimateTools as at
from TsTools import generalPlots as genp
from TsTools import PPKfunctions as ppk
from TsTools import TSPosition as tsp
from TsTools import TSPositionPlus as tspp
from TsTools import Tracked as tr
from TsTools import networkProp as netp


def BurstConstellation(dat, datDic, tankCenter,
                       tankRadius, BlInPx):
    '''
    INPUT:
        dat
            DataPlus objects from TSPositionPlus
        BlInPx float
            the mean body length in pixel
            (needed to exclude burst positions where the individual
                is too close to the wall)
    '''
    
    pTime, pN, _ = dat.pos.shape
    pNbursts = len(dat.burstStart)
    posAtStart = dat.pos[dat.burstStart]
    velAtStart = dat.vel[dat.burstStart]
    velBursterEnd = dat.vel[dat.burstEnd, dat.bursterID]
    velBursterStart = dat.vel[dat.burstStart, dat.bursterID]
    posBursterStart = dat.pos[dat.burstStart, dat.bursterID]
    headingBursterStart = np.arctan2(velBursterStart[:, 1], velBursterStart[:, 0])
    headingBursterEnd = np.arctan2(velBursterEnd[:, 1], velBursterEnd[:, 0])
    
    # position relative to burster
    relPos2Burster = posAtStart - posBursterStart[:, None]
    relVel2Burster = velAtStart - velBursterStart[:, None]
    
    # # exclude bursterID from neighbour-lists
    # thereNNs = np.where(np.sum(relPos2Burster, axis=-1) != 0)
    # relPos2Burster = relPos2Burster[thereNNs].reshape(pNbursts, pN-1, 2)
    # relVel2Burster = relVel2Burster[thereNNs].reshape(pNbursts, pN-1, 2)
    
    # Create mask with True=Voronoi neighbor of bursting individual  
    mask = np.zeros((pNbursts, pN), dtype=bool)
    for i in range(pNbursts):
        NNlist = netp.VoroNeighborList(posAtStart[i])[dat.bursterID[i]]
        # now only take the NN and set the other data to None 
        mask[i, NNlist] = True
    relPos2Burster[~mask] = np.nan
    relVel2Burster[~mask] = np.nan
    
    
    # exclude too close to tankwall
    dist2WallBurster = posBursterStart - tankCenter[None, :]
    dist2WallBurster = tankRadius - np.sqrt(np.sum(dist2WallBurster**2, axis=1))
    okDist2Wall = dist2WallBurster > BlInPx*3
    # print('okDist2Wall: ', np.sum(okDist2Wall), 'relPos2Burster: ', len(relPos2Burster), len(relPos2Burster[okDist2Wall]))
    # exclude from the return values
    gen.append2key(datDic, 'relPos2Burster', relPos2Burster[okDist2Wall])
    gen.append2key(datDic, 'relVel2Burster', relVel2Burster[okDist2Wall])
    gen.append2key(datDic, 'posBursterStart', posBursterStart[okDist2Wall])
    gen.append2key(datDic, 'headingBursterStart', headingBursterStart[okDist2Wall])
    gen.append2key(datDic, 'headingBursterEnd', headingBursterEnd[okDist2Wall])


def threeZoneForceSlow(r_r, r_o, r_a, relPosNN, relVelNN):
    absolute = lambda vec: np.sqrt(np.sum(vec**2, axis=-1))
    Nsam, N, _ = relPosNN.shape
    dNN = np.sqrt(np.sum(relPosNN**2, axis=-1))
    forceDir = np.empty(Nsam)
    for i in range(Nsam):
        d = dNN[i]
        relPos = relPosNN[i]
        relVel = relVelNN[i]
        if np.any(d <= r_r):
            meanRepulse = -np.mean(relPos[d <= r_r], axis=0)
            force = meanRepulse
        elif np.any(d <= r_a):
            N_interact = np.sum(d<=r_a)
            N_o = np.sum(d<=r_o)
            N_a = N_interact - N_o
            # orientation
            meanOrient = np.empty(2)
            if N_o > 0:
                meanOrient = np.mean(relVel[(r_r < d) & (d <= r_o)], axis=0)
                meanOrient /= absolute(meanOrient)
            # attraction
            meanAttract = np.empty(2)
            if N_a > 0:
                meanAttract = np.mean(relPos[(r_o < d) & (d <= r_a)], axis=0)
                meanAttract /= absolute(meanAttract)
            force = np.average(np.vstack((meanOrient, meanAttract)),
                               weights=[N_o, N_a], axis=0)
        else:
            force = np.random.standard_normal(2)
        forceDir[i] = np.arctan2(force[1], force[0])
    return forceDir


@njit
def threeZoneForce(r_r, r_o, r_a, relPosNN, relVelNN):
    absolute = lambda vec: np.sqrt(np.sum(vec**2, axis=-1))
    Nsam, N, _ = relPosNN.shape
    dNN = np.sqrt(np.sum(relPosNN**2, axis=-1))
    forceDir = np.empty(Nsam)
    for i in range(Nsam):
        d = dNN[i]
        relPos = relPosNN[i]
        relVel = relVelNN[i]
        if np.any(d <= r_r):
            mask = d <= r_r
            meanRepulse = -np.sum(relPos[mask], axis=0) / np.sum(mask) 
            force = meanRepulse
        elif np.any(d <= r_a):
            N_interact = np.sum(d<=r_a)
            N_o = np.sum(d<=r_o)
            N_a = N_interact - N_o
            # orientation
            meanOrient = np.empty(2)
            if N_o > 0:
                mask = (r_r < d) & (d <= r_o)
                meanOrient = np.sum(relVel[mask], axis=0) / np.sum(mask) 
                meanOrient /= absolute(meanOrient)
            # attraction
            meanAttract = np.empty(2)
            if N_a > 0:
                mask = (r_o < d) & (d <= r_a)
                meanAttract = np.sum(relPos[mask], axis=0) / np.sum(mask) 
                meanAttract /= absolute(meanAttract)
            force = np.sum(np.vstack((N_o*meanOrient, N_a*meanAttract)),
                           axis=0) / (N_o + N_a)
        else:
            force = np.random.standard_normal(2)
        forceDir[i] = np.arctan2(force[1], force[0])
    return forceDir



def forceDifference(paras, extra_data=None):
    relPosNN, relVelNN, endDir = extra_data
    r_r = paras[0]
    r_o = r_r + paras[1]
    r_a = r_o + paras[2]
    predictDir = threeZoneForce(r_r, r_o, r_a, relPosNN, relVelNN)
    angleBetween = predictDir - endDir
    angleBetween = np.mod(angleBetween + np.pi, 2*np.pi) - np.pi 
    errors = np.abs(angleBetween)
    return np.mean(errors)

# voronoi-network code

def VoroNeighborList(pos):
    '''
    INPUT:
        pos shape=(N, 2)
            x, y positions of N agents
    '''
    N = len(pos)
    edges = DelaunayEdges(pos)
    NNlist = [[] for i in range(N)]
    for i, j in edges:
        NNlist[i] += [j]
        NNlist[j] += [i]
    return NNlist


def DelaunayEdges(pos):
    '''
    INPUT:
        pos shape=(N, 2)
            x, y positions of N agents
    '''
    tri = Delaunay(pos)
    edges = UniqueEdgesFromTriangulation(tri)
    return edges


def UniqueEdgesFromTriangulation(tri):
    '''
    Get unique edges from Delaunay triangulation for Voronoi neighbours 
    '''
    e1 = tri.simplices[:,0:2]
    e2 = tri.simplices[:,1:3]
    e3 = tri.simplices[:,::2]
    edges = np.vstack((e1,e2,e3))
    edges.sort(axis=1) # IMPORTANT for unique: [[3, 1], [1, 3]] -> [[1, 3], [1, 3]]
    # This was the original part by pawel.... but it did not work (it should work since haider uses it)
    # edges_c = np.ascontiguousarray(edges).view(np.dtype((np.void, edges.dtype.itemsize * edges.shape[1])))
    # _, idx = np.unique(edges_c, return_index=True)
    # edges_unique = edges[idx]
    edges_unique = np.unique(edges, axis=0)
    return edges_unique


def pcolor_bwr_zerowhite(f, axs, mat, xvals=None, yvals=None, cLabel=None,
                         cbar=None, Nticks=None, maxVal=None,
                         **kwgs):
    '''
    creates a pcolormesh plot with bwr-colormap 
    where white is exaclty at zero
    OR
    with Reds-colormap if only positive or negative data
    INPUT:
        maxVal double
            value at which the results are cut (better name: cutValue)
        kwgs dict
            keywords for colorbar creation
    '''
    # DEFAULTS
    cbar = gen.setDefault(cbar, True)
    Nticks = gen.setDefault(Nticks, 4)
    mini = np.nanmin(mat)
    maxi = np.nanmax(mat)
    maxse = np.max([np.abs(mini), maxi])
    maxVal = gen.setDefault(maxVal, maxse)
    if mini*maxi < 0: 
        c = axs.pcolormesh(mat, cmap=cm.bwr, vmin=-maxVal, vmax=maxVal)
    elif mini >= 0:
        cmap = truncate_colormap(cm.bwr, minval=0.5, maxval=1.0, n=500)
        c = axs.pcolormesh(mat, cmap=cmap, vmin=mini, vmax=maxVal)
        # c = axs.pcolormesh(mat, cmap=cm.Reds) # old version
    elif maxi <= 0:
        cmap = truncate_colormap(cm.bwr, minval=0, maxval=0.5, n=500)
        c = axs.pcolormesh(mat, cmap=cmap, vmin=-maxVal, vmax=maxi)
        # c = axs.pcolormesh(mat, cmap=cm.Blues_r) # old version
    else:
        c = axs.pcolormesh(mat, cmap=cm.Greys)
    if xvals is not None:
        plot_set_xticks(axs, Nticks, xvals)
    if yvals is not None:
        plot_set_yticks(axs, Nticks, yvals)
    if cbar:
        col = f.colorbar(c, ax=axs, **kwgs)
        if cLabel is not None:
            col.set_label(cLabel)
        ticks = col.get_ticks()
        if len(ticks) >=5:
            col.set_ticks(ticks[::2])
        return col


def truncate_colormap(cmap, minval=None, maxval=None, n=None):
    '''
    to truncate an existing colormap
    source: 
    https://stackoverflow.com/questions/18926031/how-to-extract-a-subset-of-a-colormap-as-a-new-colormap-in-matplotlib
    '''
    minval = gen.setDefault(minval, 0.0)
    maxval = gen.setDefault(maxval, 1.0)
    n = gen.setDefault(n, 100)
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def plot_set_xticks(axs, Nticks, values, noBox=None):
    # if box or pcolor plot: the tick should be in the middle of box
    if noBox is None:
        noBox = False
    # assert Nticks <= len(values), 'Nticks > len(values)'
    if Nticks > len(values):
        Nticks = len(values)
    if type(values[0]) != str:
        ticks = np.linspace(0, len(values)-1, Nticks, dtype=int)
        ticksID = ticks.copy()
        if not noBox:
            ticks = ticks.astype('float')
            ticks += 0.5 
        axs.set_xticks(ticks)
        axs.set_xticklabels(np.round(values[ticksID], 2))
    else:
        axs.set_xticks(np.arange(Nticks) + 0.5)
        axs.set_xticklabels(values, rotation='vertical')


def plot_set_yticks(axs, Nticks, values):
    # assert Nticks <= len(values), 'Nticks > len(values)'
    if Nticks > len(values):
        Nticks = len(values)
    if type(values[0]) != str:
        ticks = np.linspace(0, len(values)-1, Nticks, dtype=int)
        axs.set_yticks(ticks + 0.5)
        axs.set_yticklabels(np.round(values[ticks], 2))
    else:
        axs.set_yticks(np.arange(Nticks) + 0.5)
        axs.set_yticklabels(values)
