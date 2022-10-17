import datetime
from dataclasses import dataclass
from pathlib import Path
from typing import Union

import ephys.ephysanalysis as EP
import numpy as np
from dateutil.tz import tzlocal
import pynwb as NWB
from pynwb import NWBHDF5IO
import matplotlib.pyplot as mpl


f = "test_data/2017.05.01~S0C0~CCIV_short_000.nwb"

class ReadNWB():
    def __init__(self):
        pass

    def readfile(self, f: Union[Path, str, None]=None):

        io = NWBHDF5IO(Path(f), mode="r")
        nwbfile = io.read()
        di = nwbfile.get_intracellular_recordings()

        self.istim = nwbfile.acquisition['Ics1'].data[:]
        self.vresp = nwbfile.acquisition['Vcs1'].data[:]
        self.sample_rate = nwbfile.acquisition['Vcs1'].rate  # data sample rate in Hz
        self.protocol = nwbfile.protocol
        self.time = np.linspace(0., self.vresp.shape[0]/self.sample_rate, self.vresp.shape[0])
        # make an ephys Clamps structure
        self.EPC = EP.MakeClamps.MakeClamps()
        self.EPC.set_clamps(self.time, self.vresp, cmddata = self.istim, tstart_tdur=[0.1, 0.1], protocol=self.protocol)
        self.EPC.setNWB(nwbmode=True)

        self.EPC.getClampData()
        IVS = EP.IVSummary.IVSummary(datapath = None, altstruct=self.EPC, file=f)
        IVS.plot_mode(mode="normal")
        IVS.compute_iv()

        # IVS.plot_iv()


    def plot_traces(self):
        fig, ax = mpl.subplots(2,1)
        for i in range(self.EPC.traces.shape[0]):
            ax[0].plot(self.EPC.time, self.EPC.traces.view(np.ndarray)[i,:])
            ax[1].plot(self.EPC.time, self.EPC.cmd_wave.view(np.ndarray)[i,:])
        mpl.show()

if __name__ == "__main__":
    R = ReadNWB()
    R.readfile(f)
    # R.plot_traces()