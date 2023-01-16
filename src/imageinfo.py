"""Evaluate images in the directories. 
Specifically, look for sets of images taken with a particular objective
that are grouped together in one subdirectory,

"""

import argparse
from pathlib import Path
from typing import Union

import ephys.ephysanalysis as EP
import numpy as np
from dateutil.tz import tzlocal
from ephys.tools import tifffile
from nwbinspector import inspect_all, inspect_nwb

AR = EP.acq4read.Acq4Read()

class EvalImages():
    def __init__(self, topdir:Path):
        if not topdir.is_dir():
            raise ValueError("Input path is not a directory")
        self.recurse_dirs(topdir, evalfunc=self.examine_image)

    def recurse_dirs(self, directory, evalfunc:object):
        pdir = Path(directory)
        if pdir.is_dir():
            # print("dir: ", pdir)
            files = list(pdir.glob('*.tif'))
            ntiffs = len(files)
            if len(files) > 0:
                print("\nExamining: ", str(pdir))
                dirindex = AR.readDirIndex(pdir)
                # print(dirindex.keys())
            tifinfo = {}
            for f in files:
                fn, s = self.examine_image(filename=f)
                if fn in dirindex.keys():
                    tifinfo[fn] = {'shape': s, 
                        'objective': dirindex[fn]['objective'], 
                        'position': dirindex[fn]['transform']['pos'],
                        'scale': dirindex[fn]['transform']['scale']}

            if len(tifinfo) > 0:
                # find out what image sizes are present
                # and what objectives are present
                sizes = []
                mags = []
                for img in tifinfo.keys():

                    print('   ', img, tifinfo[img]['shape'])
                    print('            ', tifinfo[img]['objective'])
                    print('            ', tifinfo[img]['position'], tifinfo[img]['scale'])
                    if tifinfo[img]['shape'] not in sizes:
                        sizes.append(tifinfo[img]['shape'])
                    if tifinfo[img]['objective'] not in mags:
                        mags.append(tifinfo[img]['objective'])
                print("Images sizes: ", sizes)
                print("Image mags: ", mags)
                print("\n")
            for d in pdir.glob('*/'):
                if d.is_dir():
                    self.recurse_dirs(d, evalfunc)

    def examine_image(self, filename):
        with open(filename, 'rb') as fh:
            tif = tifffile.imread(fh)
        return  filename.name, tif.shape


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate Image files in a folder")
    parser.add_argument(dest="inputfile", type=str, help="input filename")
    args = parser.parse_args()
    
    EvalImages(Path(args.inputfile))
