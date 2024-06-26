"""Display data from an NWB file generated by acq4_to_NWB
    Handles opto mapping as stored by acq4_to_NWB.py
"""

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Union
import re
import matplotlib.pyplot as mpl
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import numpy as np
from pynwb import NWBHDF5IO

teststim = "Vc_LED_stim_single_000:ostimled_0"
# These regular expressions are used to parse the acquisition and stimulus names
# Note that we define these in the same way in acq4_to_NWB.py, but it is 
# possible that you might want to add a "field" to the name.
ostim_re = re.compile(
    r"(?P<protocol>[\w_]+):(?P<name>[\w\s]+),\s+sweep=(?P<sweep>\d+)"
)
estim_re = re.compile(r"(?P<protocol>[\w_]+):(?P<name>[\w\s]+)_(?P<sweep>\d+)")

testacq = "Vc_LED_stim_single_000:Ics1_0"
acq_re = re.compile(r"(?P<protocol>[\w_]+):(?P<name>[\w\s]+)_(?P<sweep>\d+)")
test2stim = "Vc_LED_stim_single_000:Vcs1_0"

# stim = stim_re.match(test2stim)
# print(stim)
# if stim is None:
#     stim=stim_re2.match(test2stim)
# print(stim.groupdict())

# exit()


@dataclass()
class acq4_data:  # dataclass to hold the data from a single acq4 protocol
    nwb_name: str = ""  # name of the nwb file
    protocol: str = ""  # name of the acq4 protocol in this data set
    mode: str = ""  # recording mode
    nsweeps: int = 0  # number of sweeps
    response: object = None  # response data  # (becomes np.ndarray later)
    command: object = None  # electrical command data
    ostim: object = None  # optical stimulus data
    ostim_timebase: object = None  # timebase for optical stimulus
    optospots: object = None  # optical spot data
    timebase: object = None
    laser_command: object = None
    laser_timebase: object = None


class ReadNWB:
    def __init__(self, filename: Union[str, Path, None] = None):
        if filename is None:
            return
        else:
            self.readfile(f=filename)

    def split_acqname(self, acqname: str):
        """_summary_
        acquisition names:
            'Vc_LED_stim_single_000:Ics1_0'
            split to protocol, name, sweep, and get the recording mode from the comments field.

        Params
        ------
            acqname (str): acquisition name

        Returns
        -------
        Tuple[str, str, str, str]: protocol, sweep, protonum, recmode
        """
        acqdict = acq_re.match(acqname).groupdict()
        recmode = self.nwbfile.acquisition[acqname].comments
        return acqdict["protocol"], acqdict["sweep"], acqdict["name"], recmode

    def split_stimname(self, stimname: str):
        """_summary_
        stimulus names:
            'Vc_LED_stim_single_000:ChR2 Stimulation, LED 470 nm, sweep=0
            split to protocol, name, sweep, and get the recording mode from the comments field.

        Params
        ------
            acqname (str): acquisition name

        Returns
        -------
        Tuple[str, str, str, str]: protocol, sweep, protonum, recmode
        """
        recmode = self.nwbfile.stimulus[stimname].comments
        stim = estim_re.match(stimname)

        stimdict = stim.groupdict()
        if stimdict["name"].startswith("opto"):
            stimdict["stimtype"]="optical"
        elif stimdict["name"].startswith("OptoSeries"):
            stimdict["stimtype"]="OptoSeries"
        else:
            stimdict["stimtype"]="electrical"
        # print("stimdict: ", stimdict)
        return (
            stimdict,
            recmode,
        )

    def readfile(self, f: Union[Path, str, None] = None):
        """_summary_

        Params
        ------
            f (Union[Path, str, None], optional): File path. Defaults to None.

        Returns
        -------
        Nothing
        """
        self.filename = f
        if not Path(f).is_file():
            print(f"File {f:s} not found.")
            return
        else:
            print("File found.")
        io = NWBHDF5IO(Path(f), mode="r")
        nwbfile = io.read()
        self.nwbfile = nwbfile
        di = nwbfile.get_intracellular_recordings()
        root_table = nwbfile.get_icephys_meta_parent_table()
        # print("Root table : ", root_table.neurodata_type)
        if root_table.has_foreign_columns():
            print("Has Foreign Columns:", root_table.has_foreign_columns())
            print("Foreign Columns:", root_table.get_foreign_columns())

        linked_tables = root_table.get_linked_tables()
        # Print the links
        for i, link in enumerate(linked_tables):
            print(
                "%s (%s, %s) ----> %s"
                % (
                    "    " * i,
                    link.source_table.name,
                    link.source_column.name,
                    link.target_table.name,
                )
            )

        # print("acquisition: ", nwbfile.acquisition)
        acquisition_names = nwbfile.acquisition.keys()
        stimulus_names = nwbfile.stimulus.keys()

        doneprots = []
        # print("\nacquistion_names: ", acquisition_names)
        # print("\nstimulus names: ", stimulus_names)
        if len(acquisition_names) == 0:
            print("No acquisition data found.")
            print(nwbfile)
            return

        # get all of the stimulation data in the nwb file
        e_stimdict = {}
        o_stimdict = {}
        o_spotdict = {}
        for stim in stimulus_names:
            stimdict, recmode = self.split_stimname(stim)
            if stimdict["stimtype"] == "electrical":
                e_stimdict[(stimdict["protocol"], stimdict["sweep"])] = stim
            elif stimdict["stimtype"] == "optical":
                o_stimdict[(stimdict["protocol"], stimdict["sweep"])] = stim
            elif stimdict["stimtype"] == "OptoSeries":
                o_spotdict[(stimdict["protocol"], stimdict["sweep"])] = stim
        # electrical stimulus (currents in current-clamp, voltages in voltage clamp)
        # and optical stimuli (if they exist) are paired within a protocol with the
        # acquistion of voltage or current.
        # The pairing is uniquely matched up by protocol name (unique) and sweep number.

        self.protocols = (
            {}
        )  # dict pointing to each recording/stimulating dataset for a given protocol

        for acq in acquisition_names:
            protocol, sweep, acqname, recmode = self.split_acqname(acq)
            # print("acq: ", acq, "acqname: ", acqname, "protocol: ", protocol, "recmode: ", recmode)
            if protocol not in self.protocols:
                self.protocols[protocol] = acq4_data(
                    nwb_name=acq, protocol=protocol, mode="CC"
                )  # initialize storage

        # print("\nPROTOCOLS:\n", self.protocols)
        # we add the protocols to a list - note that the initial list have empty
        # data (None) until we populate
        for k, protocol in self.protocols.items():
            # print("protocol: ", protocol)
            self.add_protocol(nwbfile, protocol, e_stimdict, o_stimdict, o_spotdict)

    def add_protocol(self, nwbfile, protocol, e_stimdict: dict, o_stimdict: dict, o_spotdict: dict):
        """add the data from this protocol to the dataclass object"""
        # print("\nNWB_FILE: ", dir(nwbfile))
        recmode = nwbfile.acquisition[protocol.nwb_name].comments
        acquisition_names = nwbfile.acquisition.keys()
        stimulus_names = nwbfile.stimulus.keys()
        # print("\nstimulus names: ", stimulus_names)
        # print("\ne stim dict: ", e_stimdict)
        # print("\no stim dict: ", o_stimdict)
        # look up the matching e and o stimuli

        for i, acq in enumerate(acquisition_names):
            if not acq.startswith(
                protocol.protocol
            ):  # filter to get data from our current protocol
                # print("\nacq; ", acq, "protocol: ", protocol)
                continue
            # print("\nacq name: ", acq, protocol.nwb_name)
            protocol_name, sweep, acqname, recmode = self.split_acqname(acq)
            estim = None
            if (protocol_name, sweep) in e_stimdict:
                estim = e_stimdict[(protocol_name, sweep)]
            optostim = None
            if (protocol_name, sweep) in o_stimdict:
                # print("getting optostim")
                optostim = o_stimdict[(protocol_name, sweep)]
            optospot = None
            if (protocol_name, sweep) in o_spotdict:
                optospot = o_spotdict[(protocol_name, sweep)]
            # print("estim: ", estim, "\noptostim: ", optostim)
            # print("add_protocol: ", optospot )
            
            if recmode in ["CC", "IC"]:  # get command and response data
                response = np.array(nwbfile.acquisition[acq].data)
                command = np.array(nwbfile.stimulus[estim].data)
                try:
                    optostimdata = np.array(nwbfile.stimulus[optostim].data)
                    optospotdata = np.array(nwbfile.stimulus[optospot].data)
                except:
                    optostimdata = None
                    optospotdata = None
            
            elif recmode == "VC":
                response = np.array(nwbfile.acquisition[acq].data)
                command = np.array(nwbfile.stimulus[estim].data)
                # print("vc ostim: ", protocol.ostim)
                # print("vc optospots: ", protocol.optospots)
                # print("Stimulus:\n", nwbfile.stimulus.keys())

                if optostim is not None:
                    optostimdata = np.array(nwbfile.stimulus[optostim].data)
                    location = eval(nwbfile.stimulus[optospot].site.location)  # location of the spot
                    optospotdata = np.array(nwbfile.stimulus[optospot].data)
                else:
                    optostimdata = None
                    optospotdata = None
                    location = None
            if protocol.response is not None:
                # print("protocol.response : ", i, protocol.response.shape)
                # print("response: ", i, response.shape)
                protocol.response = np.vstack((protocol.response, response))
            else:
                protocol.response = response
            if protocol.command is not None:
                protocol.command = np.vstack((protocol.command, command))
            else:
                protocol.command = command

            if optostimdata is not None:
                ostim_timebase = np.linspace(
                    0.0,
                    optostimdata.shape[0] / nwbfile.stimulus[optostim].rate,
                    optostimdata.shape[0],
                )
                if protocol.ostim is not None:
                    protocol.ostim = np.vstack((protocol.ostim, optostimdata))
                    protocol.optospots = np.vstack((protocol.optospots, location))
                    protocol.ostim_timebase = np.vstack((protocol.ostim_timebase, ostim_timebase))
                else:
                    protocol.ostim = optostimdata
                    protocol.ostim_timebase = ostim_timebase
                    protocol.optospots = location

        self.sample_rate = self.nwbfile.acquisition[
            protocol.nwb_name
        ].rate  # data sample rate in Hz
        self.protocol = nwbfile.protocol
        # print("Response shape: ", protocol.response.shape)
        if protocol.response.ndim == 1:  # handle special case where there is just one trace.
            protocol.response = protocol.response.reshape(1, protocol.response.shape[0])
            protocol.command = protocol.command.reshape(1, protocol.command.shape[0])
            if protocol.ostim is not None:
                protocol.ostim = protocol.ostim.reshape(1, protocol.ostim.shape[0])
                protocol.ostim_timebase = protocol.ostim_timebase.reshape(
                    1, protocol.ostim_timebase.shape[0]
                )

        protocol.timebase = np.linspace(
            0.0, protocol.response.shape[1] / self.sample_rate, protocol.response.shape[1]
        )
        self.nwbfile = nwbfile
        self.response = protocol.response
        self.plot_traces(recmode, protocol, optostim)

    def onpick(self, event):
        # print("onpick event: ", event.ind)
        # print(self.response.shape)
        for i in event.ind:
            if i > 0:
                continue
            print(i, self.response[event.ind[i], :].shape)
            self.trace_plot[0].set_ydata(self.response[event.ind[i]])
            self.ELabel.set_text(f"Trace {event.ind[i]:d}")
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def plot_traces(self, recmode: str, protocol: acq4_data, optostim: object = None):
        """Minimal plotter to show the data in an nwb file generated by acq4_to_NWB.py"""

        fig, axs = mpl.subplot_mosaic(
            [["A", "D"], ["A", "D"], ["A", "D"], ["A", "D"], ["B", "E"], ["C", "F"]],
            constrained_layout=True,
            figsize=(8, 10),
        )
        self.fig = fig
        self.axs = axs
        for ax in ["A", "B", "C", "E", "F"]:
            axs[ax].spines["right"].set_visible(False)
            axs[ax].spines["top"].set_visible(False)
        fig.suptitle(
            f"{str(self.filename):s}  {protocol.nwb_name:s}", fontsize=8, fontweight="demibold"
        )
        if protocol.protocol.startswith(("CCIV", "Ic", "IC")):
            delta_I = 0.075
        elif protocol.protocol.startswith(("Vc", "VC", "Map")):
            delta_I = 100e-12
        else:
            raise ValueError(f"Protocol starting with : {protocol.protocol:s} not handled.")
        # print("protocol reponse shape: ", protocol.response.shape)
        for i in range(protocol.response.shape[0]):
            axs["A"].plot(protocol.timebase, protocol.response[i, :] + i * delta_I, linewidth=0.33)
            axs["B"].plot(protocol.timebase, protocol.command[i, :], linewidth=0.33)

            if protocol.ostim is not None:
                axs["C"].plot(protocol.ostim_timebase[i, :], protocol.ostim[i, :])
        axs["F"].plot(protocol.timebase, np.mean(protocol.response, axis=0), linewidth=0.5)
        self.trace_plot = axs["E"].plot(
            protocol.timebase, protocol.response[0, :], "g-", linewidth=0.5
        )
        self.ELabel = self.axs["E"].text(
            0.05, 0.95, f"Trace: {0:d}", ha="left", va="top", transform=self.axs["E"].transAxes
        )
        axs["A"].set_xlabel("T (s)")
        axs["F"].set_xlabel("T (s)")
        if recmode.lower() in ["cc", "ic"]:
            axs["A"].set_ylabel("V (V)")
            axs["B"].set_ylabel("I (A)")
            axs["F"].set_ylabel("V (V)")
        elif recmode.lower() in ["vc"]:
            axs["A"].set_ylabel("I (A)")
            axs["B"].set_ylabel("V (V)")
            axs["F"].set_ylabel("mean I (A)")
        if protocol.ostim is not None:
            axs["C"].set_xlabel("T (s)")
            axs["C"].set_ylabel("Intensity (AU)")

        # draw a map of the points by current
        if protocol.ostim is not None:
 
            opos = protocol.optospots
            # .reshape(len(self.opto.control), 2)
            spotsize = self.nwbfile.stimulus[optostim].comments.split(":")[-1]
            osize = float(spotsize) * 1e3  # bring to match units
            patches = []
            # create a circle to show the position of each stimulus
            for isp in range(opos.shape[0]):
                circle = mpatches.Circle(
                    (opos[isp, 0] / 1e3, opos[isp, 1] / 1e3),
                    radius=osize / 2.0,
                    color="blue",
                    ec="none",
                )
                patches.append(circle)
            collection = PatchCollection(patches, cmap=mpl.cm.cool, alpha=0.6)
            istim = int(0.5 * self.sample_rate)
            onems = int(0.001 * self.sample_rate)  # self.EPC.sample_interval)
            twentyms = 20 * onems
            colors = np.sum(protocol.response[:, istim : istim + twentyms], axis=0) - np.mean(
                protocol.response[:, istim - twentyms : istim], axis=0
            )  # , np.linspace(0, 1, len(patches))
            colors = colors / np.min(colors)
            collection.set_array(colors)
            collection.set_picker(True)
            axs["D"].add_collection(collection)
            axs["D"].set_facecolor("white")
            axs["D"].set_aspect("equal")
            axs["D"].scatter(
                np.mean(opos[:, 0]) / 1e3, np.mean(opos[:, 1] / 1e3), s=0.1, marker="o"
            )
            fig.canvas.mpl_connect("pick_event", self.onpick)

        mpl.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Display an NWB file")
    parser.add_argument(dest="inputfile", type=str, help="input filenbame")
    args = parser.parse_args()
    if not Path(args.inputfile).is_file():
        print(f"File {args.inputfile} not found.")
        exit(1)
    print("trying")
    R = ReadNWB(args.inputfile)
    # R.plot_traces()
