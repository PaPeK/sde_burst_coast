import os
import multiprocessing as mp
import numpy as np
from functools import partial
import glob
import sys
import time as pytime
import h5py
import pickle
import shutil
import matplotlib
from pathlib import Path
import SwarmDynByPy as sd
import general as gen

if __name__ == '__main__':
    matplotlib.use('Agg')

def mergeSimulations(params, outpath, paratuple):
    p_out = Path(outpath)
    dic = params.copy()
    dic.pop('para_name2')
    dirname = get_swarmdyn_out_groupname(dic, paratuple)
    f_h5 = p_out / ('out_' + dirname + '.h5')
    pattern = 'out_{}_{}*.h5'.format(dirname, params['para_name2'])
    f_splits = list(p_out.glob(pattern))
    create_or_check_swarmdyn_hdf5(outpath, dirname, dic, verb=False)

    with h5py.File(str(f_h5), 'r+') as h5f:
        for f_split in f_splits:
            with h5py.File(str(f_split), 'r+') as h5f_split:
                n_datasets = gen.h5AllDatasetNames(h5f_split, verbose=False)
                for n_dataset in n_datasets:
                    d_name = '{}/{}'.format(dirname, n_dataset.split('/')[-1])
                    data = h5f_split[n_dataset]
                    gen.h5ExtendDataSet2(h5f, d_name, data, g_name=None)
            f_split.unlink() # avoids redundant data


def increaseDictValues(dic, inc):
    for key in dic.keys():
        dic[key] += inc


def get_swarmdyn_out_groupname(dic, para_vals):
    '''
    checks in path if the value is close to an already existing value
    if yes: the string of this close value is used in groupname
        no: the string representation is used
    '''
    files = glob.glob( os.path.join(dic['path'], 'out*h5') )
    para_valsStr = []
    for i, para in enumerate(para_vals):
        para_nam = dic['para_name' + str(i)]
        files_filtered = [f for f in files if para_nam in f] # necessary
        there = []
        if len(files_filtered) > 0:
            paravals = np.array([gen.string_NrAfterWord(fil, para_nam) for fil in files_filtered])
            there = np.where(paravals == None)
            paravals = np.delete(paravals, there)
            paravals = np.unique(paravals)
            there = np.where( np.isclose(para, np.array(paravals, dtype=float)) )[0]
        if len(there) > 0:
            para_valsStr.append(paravals[there[0]])
        else:
            para_valsStr.append(str(para))

    name = ''
    for i, para in enumerate(para_vals):
        key_name = 'para_name' + str(i)
        if key_name in dic.keys():
            name += dic[key_name] +  '{}_'.format(para_valsStr[i])
        else:
            raise ValueError('para_name{} not in params.keys()'.format(i))
    if name == '':
        name = 'run_base_paras'
    else:
        name = name[:-1]  # exclude last "_" from file_name
    return name


def get_swarmdyn_out_dirname(dic, para_vals):
    groupname = get_swarmdyn_out_groupname(dic, para_vals)
    return 'out_' + groupname

def get_swarmdyn_out_name(dic, para_vals):
    name = get_swarmdyn_out_dirname(dic, para_vals)
    return os.path.join(dic['path'], name + '.h5')


def names2dict(names, pre = None):
    if pre is None:
        pre = ''
    keys = [pre + key for key in names]
    vals = list(range(len(keys)))
    dic = dict(zip(keys, vals))
    assert len(dic.keys()) == len(keys), 'len(dic.keys()) != len(keys) => probably redundant key'
    return dic


def get_out_swarm_dict(pre = None):
    keys0 = ['N',
             'pol_order',
             'L_norm',
             'avg_x[0]',
             'avg_x[1]',
             'avg_v[0]',
             'avg_v[1]',
             'avg_speed',
             'avg_vsquare',
             'elongation',
             'aspect_ratio',
             'a_com_maxIID',
             'Area_ConHull',
             'NND',
             'IID',
             'ND',
            ]
    dic = names2dict(keys0, pre=pre)
    return dic


def get_out_swarm_pred_dict(pre = None):
    keys0 = ['N_clu',
             'clu_avg_s',
             'clu_avg_vsquare',
             'clu_dpi',
             'pred.kills',
            ]
    dic = names2dict(keys0, pre=pre)
    return dic


def get_out_swarm_fishNet_dict(pre = None):
    keys0 = ['NfrontClose',
             'distNet',
             'Nfront',
             'Ndead',
            ]
    dic = names2dict(keys0, pre=pre)
    return dic


def join_enumerated_dicts(dic0, dic1):
    Nkeys0 = len(dic0.keys())
    assert Nkeys0 - 1 == max(dic0.values()), 'len(dic0.keys()) - 1 != dic0.values() => Wrong values'
    dic_new = dic0.copy()
    dic11 = dic1.copy()
    for key in dic11.keys(): 
        dic11[key] += Nkeys0
    dic_new.update(dic11)
    return dic_new


def get_out_part_dict(pre = None):
    keys0 = ['x0',
             'x1',
             'v0',
             'v1',
             'fitness',
             'force0',
             'force1',]
    dic = names2dict(keys0, pre=pre)
    return dic


def get_out_pred_dict(pre = None):
    keys0 = ['x0',
             'x1',
             'v0',
             'v1',
             'force0',
             'force1',]
    dic = names2dict(keys0, pre=pre)
    return dic


def h5CompareGrAttr2Dict(h5file, gname, dicce, verb=False):
    '''
    compares group values of group attributes of h5-file with
    values of keys of given dictionary 'dicce'
    '''
    if gname in list(h5file.keys()):
        grp = h5file.require_group(gname)
        # compare attributes:
        count = 0
        for i in list(grp.attrs.keys()):
            if i in dicce.keys():
                if type(dicce[i]) in [float, int] and i not in ["output_mode"]:
                    if dicce[i] != grp.attrs[i]:
                        print('MISMATCH OF ', i, ': ', dicce[i], '!=', grp.attrs.get(i))
                        count += 1
        if count == 0:
            if verb:
                print('all relevant attribtes MATCH')
        else:
            sys.exit()
    else:
        print('GROUP not in FILE')
        sys.exit()


def existingSamplesH5(path, params):
    dic = params.copy()  # to not modify original part
    dic.pop('para_name2', None)
    paratuple = [dic['para_values0'][0], dic['para_values1'][0]]
    dirname = get_swarmdyn_out_groupname(dic, paratuple)
    f_h5 = os.path.join(path, 'out_' + dirname + '.h5')
    existing_samples = 0
    if os.path.exists(f_h5):
        with h5py.File(f_h5, 'r+') as f:
            existing_samples = f['/' + dirname + '/swarm'].shape[0]
    return int(existing_samples)


def create_or_check_swarmdyn_hdf5(path, dirname, dic, verb=False):
    f_h5 = os.path.join(path, 'out_' + dirname + '.h5')
    Step_mean = (int((int(dic["time"]/dic["dt"])-1)/int(dic["output"]/dic["dt"])) -
                 int(int(dic["trans_time"]/dic["dt"])/int(dic["output"]/dic["dt"])) +
                 int(0 == dic["trans_time"]/dic["dt"]/int(dic["output"]/dic["dt"]) % 1))
    Step_pred = 1*Step_mean
    if dic["pred_time"] > dic["trans_time"]:
        Step_pred = (int((int(dic["time"]/dic["dt"])-1)/int(dic["output"]/dic["dt"])) -
                     int(int(dic["pred_time"]/dic["dt"])/int(dic["output"]/dic["dt"])) +
                     int(0 == dic["pred_time"]/dic["dt"]/int(dic["output"]/dic["dt"]) % 1))
    try:    # try creation of file and datasets
        with h5py.File(f_h5, 'w-') as f:
            grp = f.create_group(dirname)
            for i in list(dic.keys()):      # create attributes of group
                if type(dic[i]) != str:     # using str as attribute value is complicated
                    grp.attrs.create(i, dic[i])
            Nouts_swarm = len(get_out_swarm_dict())
            dims = np.array([0, Step_mean, Nouts_swarm])
            gen.h5createDataset(f, dirname, 'swarm', dims)
            if dic['Npred'] > 0:
                if dic["output_mode"] == 1:
                    Nouts = len(get_out_pred_dict())
                    dims = np.array([0, Step_pred, Nouts])
                    gen.h5createDataset(f, dirname, 'pred', dims)
                    gen.h5createDataset(f, dirname, 'predD', dims)
                if dic['Npred'] == 1 and dic['BC'] == -1:
                    Nouts_swarm_pred = len(get_out_swarm_pred_dict())
                    dims = np.array([0, Step_pred, Nouts_swarm_pred])
                    gen.h5createDataset(f, dirname, 'swarm_pred', dims)
                    gen.h5createDataset(f, dirname, 'swarm_predD', dims)
                else:
                    Nouts_swarm_pred = len(get_out_swarm_fishNet_dict())
                    dims = np.array([0, Step_pred, Nouts_swarm_pred])
                    gen.h5createDataset(f, dirname, 'swarm_fishNet', dims)
                    gen.h5createDataset(f, dirname, 'swarm_fishNetD', dims)
                dims = np.array([0, dic["N"], 2])
                gen.h5createDataset(f, dirname, 'endD', dims)
            dims = np.array([0, dic["N"], 2])
            gen.h5createDataset(f, dirname, 'end', dims)
            dims = np.array([0, dic["N"], 5])
            gen.h5createDataset(f, dirname, 'start', dims)
            if dic["output_mode"] == 1:
                Nouts = len(get_out_part_dict()) 
                dims = np.array([0, Step_mean, dic["N"], Nouts])
                gen.h5createDataset(f, dirname, 'part', dims)
                dims = np.array([0, Step_pred, dic["N"], Nouts])
                gen.h5createDataset(f, dirname, 'partD', dims)
        if verb:
            print('PPK: created ', f_h5)
        existing_samples = 0
    except:  # load file and check if same attributes where used
        if verb:
            print('PPK: try opening')
        with h5py.File(f_h5, 'r+') as f:
            h5CompareGrAttr2Dict(f, dirname, dic)
            existing_samples = f['/' + dirname + '/swarm'].shape[0]

    # if dic["output_mode"] == 1:
    #     try:
    #         with h5py.File(f_h5, 'w-') as f:
    #             grp = f.create_group(dirname)
    #             for i in list(dic.keys()):     # create attributes of group
    #                 if type(dic[i]) in [float, int]:
    #                     grp.attrs.create(i, dic[i])
    #             # CREATE PARTICLE DSET:
    #             Nouts = len(get_out_part_dict()) 
    #             dims = np.array([0, Step_mean, dic["N"], Nouts])
    #             gen.h5createDataset(f, dirname, 'part', dims)
    #             dims = np.array([0, Step_pred, dic["N"], Nouts])
    #             gen.h5createDataset(f, dirname, 'partD', dims)
    #             Nouts = len(get_out_pred_dict())
    #             dims = np.array([0, Step_pred, Nouts])
    #             gen.h5createDataset(f, dirname, 'pred', dims)
    #             gen.h5createDataset(f, dirname, 'predD', dims)
    #     except:
    #         with h5py.File(f_h5, 'r+') as f:
    #             h5CompareGrAttr2Dict(f, dirname, dic)
    return existing_samples


def MultipleRunSwarmdyn(params, runs, path, paratuple, verb=False):
    '''
    sequantially runs swarmdyn with same parameters
    ?params[out_h5] == 1?
        -yes: ?existing hdf5-datasets with same name?
            ->yes: ?same paras used?
               ->yes: use same dataset
               ->no: stop simulation
            -> no: create new hdf5-dataset
        -> no: save in new txt-file (with index corresponding to run)
    INPUT:
        params dict()
            contains all parameters needed for simulation
        runs int
           # of repeats of same parameters 
        path str
            output directory
        paratuple list
            para_values which shall be changed
            corresponding para_names are defined in params['para_name0']
            params['para_name1']
    OUTPUT:
        [time, paratuple]
            only needed to estimate computation time for specific
            parameters
    '''
    dic = params.copy()  # to not modify original part
    dirname = get_swarmdyn_out_groupname(dic, paratuple)
    if 'selectionLine' in dic.keys():
        sd.selectionLineParameter(dic, dic['selectionLine'])
    # needs repetition otherwise overwritten by line above
    dirname = get_swarmdyn_out_groupname(dic, paratuple)
    parameter_dependence(dic)
    dic['path'] = path
    dic['fileID'] = dirname
    if verb:
        print('Running {} runs for '.format(runs) + dirname)
        print(dic['path'])

    if dic['output_mode'] != 2:
        if dic['out_h5'] == 0:     # txt output, runs sorted in directories
            path = path + dirname + os.sep
            add = len(glob.glob(path + "*"))    # glob.glob(pattern) returns list of paths matching pattern
            if not os.path.exists(path):
                os.mkdir(path)
                add = 0
        else:     # hdf5 output, file extended for new run
            _ = create_or_check_swarmdyn_hdf5(path, dirname, dic)

    t0 = pytime.time()
    for i in range(runs):
        if dic['out_h5'] == 0 and dic['output_mode'] != 2: # txt output, runs sorted in directories
            dic['path'] = path + str(i+add) + os.sep
            if verb:
                print("run: ", i, dic['path'])
            if not os.path.exists(dic['path']):
                os.mkdir(dic['path'])
        command = sd.dic2swarmdyn_command(dic)
        os.system(command)
    t1 = pytime.time()
    return [t1 - t0, paratuple]


def DoubleSimulationScan(params, runs=20, outpath='./',
                         scantype='seq', no_processors=24, verb=None,
                         useAllProcesses=None):
    if verb is None:
        verb = False
    if useAllProcesses is None:
        useAllProcesses = False 
    dic = params.copy()  # to not modify original part
    # ensure no redundant computation: also causing Bug in multiprocessing
    dic['para_values0'] = np.unique(dic['para_values0'])
    dic['para_values1'] = np.unique(dic['para_values1'])
    #  Runs Simulation and Saves it in Folder
    checkRuns = True
    if not os.path.exists(outpath):
        os.mkdir(outpath)
        checkRuns = False 
    # check how many runs have already been done
    # ATTENTION: assumes that Scan was finished -> 1 hdf5-file represents all
    Nsplit = 1
    if checkRuns:
        existing_runs = existingSamplesH5(outpath, dic)
        runs -= existing_runs
    if(scantype == 'seq'):
        for p in dic['para_values0']:
            for q in dic['para_values1']:
                MultipleRunSwarmdyn(dic, runs, outpath, [p, q])
    else:
        total_no_p = len(dic['para_values0']) * len(dic['para_values1'])
        if useAllProcesses:
            Nsplit = int( no_processors / total_no_p )
            if(Nsplit > 1):
                no_processors = Nsplit * total_no_p 
                dic['para_name2'] = 'Nsplit'
                dic['para_values2'] = list(range(Nsplit)) 
                runs /= Nsplit
                if (runs % 1) > 0: # to ensure: runs * Nsplit >= runs_original
                    runs = runs + 1
                runs = int(runs)
        else:
            if(no_processors > total_no_p):
                no_processors = total_no_p
        parallel_run_func = partial(MultipleRunSwarmdyn, dic,
                                    runs, outpath)
        comp_pool = mp.Pool(no_processors)
        para_values = []
        # TODO: Stopped HERE!
        for i in dic['para_values0']:
            pa_tuples = list(zip([i]*len(dic['para_values1']), dic['para_values1']))
            if 'para_name2' in dic.keys():
                for j in dic['para_values2']:
                    para_values += [list(tup) + [j] for tup in pa_tuples]
            else:
                para_values += pa_tuples
        if verb:
            print('para_values: ', para_values)
        out = comp_pool.map(parallel_run_func, para_values, chunksize=1)   # para_values passed as tuple
        # JOIN SPLITTED RESULTS 
        if useAllProcesses and Nsplit > 1:
            para_values = [vals[:-1] for vals in para_values]  # remove para_value2
            para_values = np.unique(para_values, axis=0)
            parallel_merge = partial(mergeSimulations, dic,
                                     outpath)
            _ = comp_pool.map(parallel_merge, para_values, chunksize=1)
        comp_pool.terminate()
    return

def Scanner(outpath, params, scantype, no_processors, runs):
    # ensure reproducability
    f_current = Path(os.path.realpath(__file__))
    d_current = f_current.parent
    gen.copyIfNeeded(str(d_current / 'src'), outpath)
    gen.copyIfNeeded(str(d_current / 'swarmdyn'), outpath)
    gen.copyIfNeeded(str(d_current / 'Makefile'), outpath)
    gen.copyIfNeeded(str(d_current / 'SwarmDynByPy.py'), outpath)
    gen.copyIfNeeded(str(f_current), outpath)

    t_start = pytime.time()
    # params used:
    for key in params:
        print(key + " = " + str(params[key]))
    # Run Simulation and Move Files
    ############################################################
    DoubleSimulationScan(params, runs=runs, scantype=scantype,
                         no_processors=no_processors, outpath=outpath,
                         useAllProcesses=True)

    t_end = pytime.time()
    print(outpath + ' finished in ', (t_end - t_start)/60, ' minutes')



def parameter_dependence(dic):
    para_names = [dic['para_name0'], dic['para_name1']]
    if dic['alg_range'] > dic['att_range']:
        if 'att_range' in para_names:
            dic['alg_range'] = dic['att_range']
        else:
            dic['att_range'] = dic['alg_range']
    if dic['rep_range'] > dic['alg_range']:
        if 'alg_range' in para_names:
            dic['rep_range'] = dic['alg_range']
        else:
            dic['alg_range'] = dic['rep_range']
    if dic['env_strength'] != dic['soc_strength']:
        dic['env_strength'] = dic['soc_strength']


def run_func4list(func, outpaths, base_para, para_changes,
                  scantype, no_processors, runs):
    '''
    Helper functions which just repeatedly calls the
    same functions with modified parameters
    1. "base_para" contains the default parameters
    1.1. "base_para" is copied so the original is not modified
    2. "para_changes" is a dictionary with
        keys = parameters of "base_para" to be changed
        and
        values = list of the parameter values to substitute
    3. "outpaths" list of dictionary names where the repspective results will be
        to be saved
    '''
    for i in range(len(outpaths)):
        if not os.path.exists(outpaths[i]):
            os.mkdir(outpaths[i])
        newpara = base_para.copy()
        for k in para_changes.keys():
            assert len(para_changes[k]) == len(outpaths), (
                    'len(para_changes[{}]) != len(outpaths), {}'.format(k, para_changes[k]))
            if para_changes[k][i] is not None:
                newpara[k] = para_changes[k][i]
        f_paras = Path(outpaths[i]) / 'paras_independent.pkl'
        time_stamp = pytime.strftime('%y%m%d%H')
        f_para_stamped = f_paras.parent / ('paras_independent' + time_stamp + '.pkl')
        pickle.dump(newpara, f_paras.open('wb'))   # paras which will be used
        pickle.dump(newpara, f_para_stamped.open("wb"))
        print((("Outpath = {}, Processors = {}, Runs = {}, "
               ).format(outpaths[i], no_processors, runs)))
        func(outpaths[i], newpara, scantype,
             no_processors, runs)


possible_modes = ['burst_coast', 'sinFisher', 'mulFisher', 'natPred', 'natPredNoConfu', 'exploration']
# This function is executed if the script is called from terminal e.g. using > python RunSimulationScan.py
if __name__ == '__main__':
    ###########################################################
    #
    #How this program works:
    #
    #The ./swarmdyn is run while scanning two parameters of your choice. The output is saved to a local
    #folder and then copied to another (network) folder after the scan is completed. You can choose how
    #many runs you would like to do, and at what number you would like to start at (default is 0).
    #
    #Additionally you may choose to run the program in parallel, using as many cores as you would like.
    ###########################################################

    # to make sure new code is used!!!
    os.system("make")

    t0 = pytime.time()
    #Scan Values
    scantype = 'para'# 'para': parallel; 'seq':sequential (Better for debugging)
    no_processors = 38  # 66: for itb cluster, 40: for CHIPS
    runs = 2   # 40 

    f_current = Path(os.path.realpath(__file__))
    d_current = f_current.parent
    outpaths = [
                str(d_current / 'ScanOutput'),
                ]
    for i in range(len(outpaths)):
        if (not outpaths[i].endswith('/')):
            outpaths[i] += "/"

    rep = len(outpaths) # repeat: create the correct length of lists
    # Simulation BASE-Parameter Values
    pred_time = 20
    record_time = 100
    # AllOptiModes = ['WithAlphaSimplestI', 'WithAlphaSimplestII', 'WithAlphaAsPNAS']
    errOrder = False

    para_changes = dict()
    base_para = sd.get_base_params(pred_time, record_time)
    # non-explorativ: parameters based on fitting and optimization
    base_para['selectionLine'] = 0
    para_changes['para_name0']   = rep * ['prob_social']   # prob to do social update
    para_changes['para_values0'] = rep * [np.arange(0.05, 0.96, step=0.1)]
    para_changes['para_name1'] = rep * ['soc_strength']
    para_changes['para_values1'] = rep * [np.arange(70, 300, step=100)]

    for i in range(len(outpaths)):
        if (not outpaths[i].endswith('/')):
            outpaths[i] += "/"
    para_changes['path'] = outpaths


    # create Data:
    run_func4list(Scanner, outpaths, base_para, para_changes,
                  scantype, no_processors, runs)

    t1 = pytime.time()
    print('total time:', t1-t0)

    # Analyze Data, create Plots
    import EvaluateScan as es
    pic_noProcessors = int(1 * no_processors)
    if pic_noProcessors > len(outpaths):
        pic_noProcessors = len(outpaths)
    comp_pool = mp.Pool(pic_noProcessors)
    comp_pool.map(es.Evaluate, outpaths)
