import numpy as np
import re
import os
from pathlib import Path
import h5py

def string_NrAfterWord(string, word):
    nums = re.compile(word + r"[+-]?\d+(?:\.\d+)?")
    wordNr = nums.search(string)
    if wordNr is None:
        nr = None
    else:
        wordNr = wordNr.group(0)
        nr = wordNr.split(word)[1]
    return nr


def h5createDataset(h5file, gname, dname, dims, v=False):
    '''
    create a dataset with apropriate chunksize
    -1 chunk saved in continous array on disk
    -regulating chunksize via 2nd dimension
        -recommendation: 1KB - 10MB
        -note: if any element of chunk accessed -> whole chunk is loaded
    -compression: save disk space while sacrificing read speed!!!
        -gzip is most common one
        -compression level between [0,9]
    '''
    maxshap = list(dims)
    for i in range(len(maxshap)):
        if dims[i] == 0:
            maxshap[i] = None      # maxshape = None means unlimited extendable

    # chunk size may not be to large
    elems = 1
    chunks = np.array(dims, dtype=int)
    chunks[0] = 1   # otherwise zero chuck size
    # assert chunks[1] > 1 or dims == 2, 'chunks[1] <= 1, thus no chunkreduction along 2nd dim possible'
    if v:
        print("1.chunks: ", chunks)
    for i in range(1, len(dims)):
        elems *= dims[i]
    while elems > 10000:
        chunks[1] /= 2
        elems /= 2
        if chunks[1] == 1:
            elems = 10000
            print("ATTENTION chunks are still large")
    if v:
        print("2.chunks: ", chunks)

    h5file.create_dataset(gname+'/'+dname, tuple(dims),
                          maxshape=tuple(maxshap), chunks=tuple(chunks),
                          dtype='double', fillvalue=0,
                          compression="gzip", compression_opts=9)


def copyIfNeeded(src, dst):
    '''
    copies file or directory to 'dst' if it does not exist
    otherwise
    do not copy
    INPUT:
        src str
            file tobe copied
        dst str
            directory in which the file/dir to be copied
    '''
    p_src = Path(src)
    p_dst = Path(dst) / p_src.parts[-1]
    if not p_dst.exists():
        if p_src.is_dir():
            shutil.copytree( src, str(p_dst) )
        elif p_src.is_file():
            shutil.copy( src, dst )


def NanAverage(data, weights, verb=None):
    '''
    as np.average you can use weights but this also neglects nans
    INPUT:
        data.shape(samples, vari) or (samples)
        weitghts.shape
    '''
    if verb is None:
        verb=False
    assert 0 < len(data.shape) and len(data.shape) < 3 , 'len(data.shape) != 1 or 2'
    if len(data.shape) == 1:
        data = data.reshape(len(data), 1)
    samples, varis = data.shape
    aver = np.empty(varis)
    for i, dat in enumerate(data.T):
        indi = ~np.isnan(dat)
        if indi.sum() > 0:
            aver[i] = np.average(dat[indi], weights=weights[indi])
        else:
            aver[i] = np.nan
            if verb:
                print('nan at vari:', i)
    return aver


def silentRemove(f_name):
    try:
        os.remove(f_name)
    except OSError as e:
        if e.errno != errno.ENOENT:
            raise


def setDefault(x, val):
    if x is None:
        x = val
    return x
