"""Display data from an NWB file generated by acq4_to_NWB
Includes analysis with IVSummary in ephys.ephysanalysis for current-clamp data


"""
import argparse
import datetime
from dataclasses import dataclass
from pathlib import Path
from typing import Union

import ephys.ephysanalysis as EP
import matplotlib.pyplot as mpl
import numpy as np
import pynwb as NWB
from dateutil.tz import tzlocal
from pynwb import NWBHDF5IO

# f = "test_data/2017.05.01~S0C0~CCIV_short_000.nwb"
# f = "test_data/2022.10.10~S1C0~CCIV_long_000.nwb"
class ReadNWB():
    def __init__(self, filename:Union[str, Path, None]=None):
        if filename is None:
            return
        else:
            self.readfile(f=filename)

    def readfile(self, f: Union[Path, str, None]=None):
        """_summary_

        Params
        ------
            f (Union[Path, str, None], optional): File path. Defaults to None.
        
        Returns
        -------
        Nothing
        """        

        io = NWBHDF5IO(Path(f), mode="r")
        nwbfile = io.read()
        di = nwbfile.get_intracellular_recordings()

        recmode = nwbfile.acquisition['Ics1'].comments

        if recmode == 'CC':  # get command and response data .. 
            self.cmd = nwbfile.acquisition['Ics1'].data[:]
            self.resp = nwbfile.acquisition['Vcs1'].data[:]
        elif recmode == 'VC':
            self.resp = nwbfile.acquisition['Ics1'].data[:]
            self.cmd = nwbfile.acquisition['Vcs1'].data[:]
            
        self.sample_rate = nwbfile.acquisition['Vcs1'].rate  # data sample rate in Hz
        self.protocol = nwbfile.protocol
        self.time = np.linspace(0., self.resp.shape[0]/self.sample_rate, self.resp.shape[0])
        # make an ephys Clamps structure
        self.EPC = EP.MakeClamps.MakeClamps()
        self.EPC.set_clamps(self.time, self.resp, 
            cmddata = self.cmd, 
            tstart_tdur=[0.1, 0.1], 
            protocol=self.protocol,
            dmode=recmode)
        self.EPC.setNWB(nwbmode=True, nwbfile=nwbfile)

        self.EPC.getClampData()

        if self.EPC.mode == "CC":
            IVS = EP.IVSummary.IVSummary(datapath = None, altstruct=self.EPC, file=f)
            IVS.plot_mode(mode="normal")
            IVS.compute_iv()
        elif self.EPC.mode == "VC":
            VCS = EP.VCTraceplot.VCTraceplot(datapath=None, altstruct=self.EPC, file=f)
            VCS.plot_vc()
            
        # IVS.plot_iv()


    def plot_traces(self):
        fig, ax = mpl.subplots(2,1)
        for i in range(self.EPC.traces.shape[0]):
            ax[0].plot(self.EPC.time, self.EPC.traces.view(np.ndarray)[i,:])
            ax[1].plot(self.EPC.time, self.EPC.cmd_wave.view(np.ndarray)[i,:])
        mpl.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Display an NWB file")
    parser.add_argument(dest="inputfile", type=str, help="input filenbame")
    args = parser.parse_args()
    
    R = ReadNWB(args.inputfile)
