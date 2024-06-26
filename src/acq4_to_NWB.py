"""
This program provides functions to tead and convert acq4 data to NWB format
that is acceptable for the Dandi repository.
It was specifically designed for the auditory cortex experiments with Kato lab,
but should be generally useful for converting other kinds of acq4 experiments.
A single file can be converted by calling ConvertFile(filename)
The conversion is checked after generating the outputfile with nwbinspector.

V0.2.0, 24 June 2024

NWBFile expects:
class pynwb.file.NWBFile(
    Required:
        session_description (str) :: a description of the session where this data was generated
        identifier (str) :: a unique text identifier for the file
        session_start_time (datetime) :: the start date and time of the recording session
        
    Optional (we try to include these):
        experimenter (tuple or list or str) :: name of person who performed experiment
        experiment_description (str) :: general description of the experiment
        institution (str) :: institution(s) where experiment is performed
        lab (str) :: lab where experiment was performed
        source_script (str) :: Script file used to create this NWB file.
        source_script_file_name (str) :: Name of the source_script file

        devices (list or tuple) :: Device objects belonging to this NWBFile
        subject (Subject) :: subject metadata
        protocol (str) :: Experimental protocol, if applicable. E.g., include IACUC protocol
        acquisition (list or tuple) :: Raw TimeSeries objects belonging to this NWBFile
        stimulus (list or tuple) :: Stimulus TimeSeries objects belonging to this NWBFile : Electrical
        opto: OptogeneticStimulusSites and waveforms that belong to this NWBFile.
        lab_meta_data (list or tuple) :: an extension that contains lab-specific meta-data
        intracellular_recordings (IntracellularRecordingsTable) :: the IntracellularRecordingsTable table that belongs to this NWBFile
        icephys_simultaneous_recordings (SimultaneousRecordingsTable) :: the SimultaneousRecordingsTable table that belongs to this NWBFile
        icephys_sequential_recordings (SequentialRecordingsTable) :: the SequentialRecordingsTable table that belongs to this NWBFile
        icephys_repetitions (RepetitionsTable) :: the RepetitionsTable table that belongs to this NWBFile
        icephys_experimental_conditions (ExperimentalConditionsTable) :: the ExperimentalConditionsTable table that belongs to this NWBFile

         units (Units) :: A table containing unit metadata

        pharmacology (str) :: Description of drugs used, including how and when they were administered. Anesthesia(s), painkiller(s), etc., plus dosage, concentration, etc.
        virus (str) :: Information about virus(es) used in experiments, including virus ID, source, date made, injection location, volume, etc.

    Optional and not included:
        file_create_date (ndarray or list or tuple or Dataset or StrDataset or HDMFDataset or AbstractDataChunkIterator or datetime) :: the date and time the file was created and subsequent modifications made
        timestamps_reference_time (datetime) ::  date and time corresponding to time zero of all timestamps; defaults to value of session_start_time
        session_id (str) :: lab-specific ID for the session
        keywords (ndarray or list or tuple or Dataset or StrDataset or HDMFDataset or AbstractDataChunkIterator) :: Terms to search over
        related_publications (tuple or list or str) :: Publication information.PMID, DOI, URL, etc. If multiple, concatenate together and describe which is which. such as PMID, DOI, URL, etc
        slices (str) :: Description of slices, including information about preparation thickness, orientation, temperature and bath solution
        data_collection (str) :: Notes about data collection and analysis.
        surgery (str) :: Narrative description about surgery/surgeries, including date(s) and who performed surgery.
        stimulus_notes (str) :: Notes about stimuli, such as how and where presented.
        analysis (list or tuple) :: result of analysis
        stimulus (list or tuple) :: Stimulus TimeSeries objects belonging to this NWBFile
        stimulus_template (list or tuple) :: Stimulus template TimeSeries objects belonging to this NWBFile
        epochs (TimeIntervals) :: Epoch objects belonging to this NWBFile
        epoch_tags (tuple or list or set) :: A sorted list of tags used across all epochs
        trials (TimeIntervals) :: A table containing trial data
        invalid_times (TimeIntervals) :: A table containing times to be omitted from analysis
        intervals (list or tuple) :: any TimeIntervals tables storing time intervals
        processing (list or tuple) :: ProcessingModule objects belonging to this NWBFile
        electrodes (DynamicTable) :: the ElectrodeTable that belongs to this NWBFile
        electrode_groups (Iterable) :: the ElectrodeGroups that belong to this NWBFile
        sweep_table (SweepTable) :: the SweepTable that belong to this NWBFile
        imaging_planes (list or tuple) :: ImagingPlanes that belong to this NWBFile
        ogen_sites (list or tuple) :: OptogeneticStimulusSites that belong to this NWBFile

        scratch (list or tuple) :: scratch data
       icephys_filtering (str) :: [DEPRECATED] Use IntracellularElectrode.filtering instead. Description of filtering used.
        ic_electrodes (list or tuple) :: DEPRECATED use icephys_electrodes parameter instead. IntracellularElectrodes that belong to this NWBFile

        
    Mapping ACQ4 to NWB:
    ====================
    Acq4 stores data in a hierarchical directory structure.
        The top level is the DAY of the experiment.
        The second level is the SLICE number.
        The third level is the CELL number.
        The fourth level is the PROTOCOL name, with an appended number in the format "_000".
            Repeats of a protocol will increment this number, although if a protocol is ended
            early, the sequence may not have fixed steps.
        Within a protocol, there is a series of subdirectories, each containing
            the data for a single sweep.
        The protocol "sweep" subdirectories are numbered according to the nesting of repititions. If there is only 
        one repitition, the subdirectory is named "000". If there are multiple repitions, the
        Subdirectories are named "000_000", "001_000", (for example showing two reps of the same stimulus).
            etc. If there is another parameter that is varied, it is added to the name, e.g. "000_000_000" (this
            is rarely if ever done).
    
    Mapping this to the NWB structure:
        1. Experimental conditions table - not used here.
        2. Repetitions Table: Should represent the repetitions of a protocol.
        3. Sequential Recordings Table: Should represent the sweeps within a repetition (usually parametric
        variation of current, laser spot position, or other stimulus parameter). This corresponds to a "protocol"
        in the acq4 structure.
        4. Simultaneous Recordings Table: Should represent the simultaneous recordings from multiple electrodes. We
        are not doing this, and it can be skipped according to the pynwb documentation.
        5. Intracellular Recordings Table: mapping recordings from multiple electrodes to the individual intracellular
        recordings. This is not used here (we are not recording from multiple electrodes).
        6. PatchClampSeries: This is the lowest level, that corresponds to a single trial or sweep in acq4 - typically
        the current or voltage trace, and the stimulus command. This is the level at which we will store the traces.

        Also included are optical stimluation data, stored in ogen structures (data and series). 

We use the NAME of each entry to parse and reassemble the data for display. See display_nwb.py 
for an example of how we use this. The name is formatted as a string that can be parsed (in regex) to extract the
protocol, data type, and sweep number. The protocol and sweep numbers are used to bind the stimulus
(optical, electrical) and recordings (electrical) together. The optogenetic data in the ogen structures similarly
has names that match and hold other parameters such as the spot location (if using the laser with a scanning
mirror), power, etc. )

        
Support::

    NIH grants:
    DC RF1 NS128873 (Kato, Manis, MPI, 2022-) Cortical circuits for the integration of parallel short-latency auditory pathways.
    DC R01 DC019053 (Manis, 2020-2025) Cellular mechanisms of auditory information processing. 

Copyright 2022-2024 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 

"""

import datetime
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Literal, Tuple, Union

import numpy as np
import pynwb as NWB
from pynwb.core import DynamicTable, VectorData
from dateutil.tz import tzlocal
from ephys import datareaders as DR
from nwbinspector import inspect_all, inspect_nwbfile
from pynwb import NWBHDF5IO


def def_lister():
    return []


@dataclass
class ExperimentInfo:
    """Data class to hold the metadata for the experiment."""

    description: str = ""
    protocol: str = ""
    time: str = ""
    experimenters: list = field(default_factory=def_lister)
    lab: str = ""
    institution: str = "UNC Chapel Hill"
    experiment_description: str = ""
    sessionid: str = "0"
    notes: str = ""
    subject: str = ""


print("acq4tonwb")


class ACQ4toNWB:
    def __init__(self, out_file_path: Union[str, Path, None] = None):
        """Convert a file from ACQ4 (www.acq4.org) format to NWB format (www.nwb.org)
        See the documentation above for the mapping between these formats.

        Args:
            out_file_path (Union[str, Path, None], optional): The path to the output files. Defaults to None.
        """
        self.AR = DR.acq4_reader.acq4_reader()  # get acq4 reader from the ephys package.

        self.out_file_path = out_file_path
        # All of the data here is single electrode, and "MultiClamp1.ma" is the
        # name of the data file holding the ephys data. The ".ma" acq4 files are held in "metaarray" format, but
        # these are basically hdf5 files. They are accompanied by .index files with additional metadata.
        # if a different electrode/amplifier is used, change [1] to [2] etc. in self.ampchannels
        self.set_amplifier_name("MultiClamp1.ma")
        self.laser_device = None
        self.LED_device = None
        self.ampchannels = [1]
        # self.manager = NWB.get_manager()
        self.ID = 0

    def _get_slice_cell(self, f: Union[str, Path]):
        """pull the slice and cell directories from the full file name
        Acq4 uses a fairly strict hierichal structure of:
        Day
            slice number 0
                cell number 0
                cell number 1
            slice number 1
                cell number 0
            ... etc.


        Args:
            f (Union[str, Path], optional): Full file path.

        Returns:
            slicen, cell (string represtations of partial path values)
        """
        f = Path(f)
        protocol = f.stem
        cellp = f.parent
        cell = cellp.stem
        slicenp = cellp.parent
        slicen = slicenp.stem
        day = slicenp.parent
        return slicen, cell

    def _get_short_name(
        self,
        f: Union[str, Path],
        expt_id: Union[str, None] = None,
        rig_id: Union[str, None] = None,
    ):
        """
        Convert a data path of this format (used by acq4):
        f = Path('/Volumes/Pegasus/all/the/leading/pathnames/2017.05.01_000/slice_000/cell_000')
        To:
        2017.05.01~S0C0 ("short form")

         Args:
            f (Union[str, Path], optional): Full file path.
            expt_id: str or None: an experiment label to put in the filename.
                ignored if None
            rig_id: str or None: a rig label to put in the filename
                ignored if None

        Returns:
            short path name string
        """
        if f is None:
            raise ValueError(f"Input file to get_shortname is NONE")
        fp = Path(f).parts
        cell = fp[-1]
        slicenp = fp[-2]
        day = fp[-3]
        if expt_id == None:
            expt_id = ""
        if rig_id == None:
            rig_id = ""
        foname = str(Path(expt_id, rig_id, day[:10], "S" + slicenp[-2:] + "C" + cell[-2:])).replace(
            "/", "~"
        )
        return foname

    def set_amplifier_name(self, amplifier_name: str = "MultiClamp1.ma"):
        """Convenience function to change the amplifier name
        This is the name of the acq4 file that holds the recorded data, and
        may also hold the stimulus files. It also has some associated metadata,
        in the .index file that is associated with the .ma file.


        Args:
            datatype (str): name of the amplifier
        """
        self.AR.setDataName(amplifier_name)

    def ISO8601_age(self, agestr, strict: bool = False):
        """Convert free-form age designators to ISO standard, e.g.:
            postnatal day 30 mouse = P30D  (or P30W, or P3Y)
            Ranges are P1D/P3D if bounded, or P12D/ if not known but have lower bound.

        Params:
            agestr (str): age string from the file

        Returns:
            str: sanitized age string
        """

        agestr = agestr.replace("p", "P")
        agestr = agestr.replace("d", "D")
        if strict:
            if agestr.find("ish") or agestr.find("?"):
                raise ValueError(
                    "Age with 'ish or '?' is not a valid age descriptor; please fix in the acq4 data."
                )
        if "P" not in agestr:
            agestr = "P" + agestr
        if "D" not in agestr:
            agestr = agestr + "D"
        if agestr == "PD":
            agestr = "P9999D"  # no age specified
        return agestr

    def find_name_in_path(self, sourcepath: Union[str, Path], name: Union[str, None] = None):
        """find_name_in_path Search the path string for a name, and return the name if found.
        This is used to find "Rig", "Experiment", "Slice", "Cell" etc. in the path string.

        Parameters
        ----------
        sourcepath : Union[str, Path]
            The full path
        name : Union[str, None], optional
            The string to search for in the path, by default None

        Returns
        -------
        _type_
            str: the name found in the path

        Raises
        ------
        ValueError
            when name is ambiguous in the path (non-unique)
        """
        if name is None:
            return ""
        pathparts = Path(sourcepath).parts
        name_id = [r for r in pathparts if r.startswith(name)]
        if name_id == []:
            name_id = ""
        elif len(name_id) > 1:
            raise ValueError("Name is ambiguous in the path: try unique pattern")
        else:
            name_id = name_id[0]
        return name_id

    def acq4tonwb(
        self,
        experiment_name: str,
        path_to_cell: Union[str, Path],
        protocols: list,
        outfilename: Union[Path, str, None] = None,
        appendmode: bool = False,
        recordingmode: Literal["IC", "CC", "VC", "I=0"] = "CC",
    ):
        """Convert one cell directory to an NWB file.
        The NWB file will contain all of the recordings from ONE cell.

        Args:
            path_to_cell (string or Path): Full path and name of the cell directory
            outfilename (Union[Path, str, None], optional): NWB output filename. Defaults to None.
            recordingmode (Literal[&quot;IC&quot;, &quot;CC&quot;, &quot;VC&quot;, &quot;I, optional):
              _description_. Defaults to 0"]="CC".

        Returns:
            None or outputfilename
                None indicates failure to read the data.
                an outputfile indicates conversion success.
        """

        rig_id = self.find_name_in_path(path_to_cell, "Rig")
        expt_id = self.find_name_in_path(path_to_cell, experiment_name)
        if outfilename is None:
            outfilename = Path(
                "test_data",
                self._get_short_name(path_to_cell, expt_id, rig_id) + ".nwb",
            )

        print("NWB filename: ", outfilename)
        # Read the acq4 index metadata for the day info, the current slice, and the current cell.
        # assumes that the currdir is the cell directory.
        info = self.AR.readDirIndex(currdir=path_to_cell.parent.parent)["."]
        self.day_index = info
        slice_index = self.AR.readDirIndex(currdir=path_to_cell.parent)["."]
        cell_index = self.AR.readDirIndex(currdir=path_to_cell)["."]

        # We should check for ".mosaic" files in the slice index. If these are present,
        # they are a pyqtgraph configuration file (renders as a dictionary)
        # that drives the acq4 mosaic module to "process" the image data.
        # The mosaic module takes the dictionary, and from the raw image
        # data, fiduciary points, and some other settings, generates a
        # an image that can be a reconstruction of the experiment.
        #
        # If there is afile with the same name as the mosaic file, then
        # is an image of the processed images etc., and we will add it to
        # the repsitory as a single RGB image.

        data_date = datetime.date.fromtimestamp(info["__timestamp__"])

        # get the cell file metadata and clean up the representation.
        age = self.ISO8601_age(info["age"])
        if "sex" not in info.keys() or info["sex"] == "":
            info["sex"] = "U"
        else:
            info["sex"] = info["sex"].upper()
        if "weight" not in info.keys():
            info["weight"] = None
        if "species" not in info.keys() or info["species"] == "":
            info["species"] = "Mus musculus"
        if age is not None:
            dob = datetime.datetime.combine(
                data_date - datetime.timedelta(days=int(age[1:-2])),
                datetime.time(),
                tzlocal(),
            )
        # dobstr = dob.strftime("%d/%m/%Y%Z")
        else:
            dob = None
        if "mouse" in info["species"] or "Mouse" in info["species"]:
            info["species"] = "Mus musculus"
        subject_id = info["animal identifier"]
        dset = Path(path_to_cell).parts
        if subject_id.strip() in [None, "?", "", "NA"]:
            subject_id = str(Path(dset[-2]))
        session_id = self._get_short_name(path_to_cell, rig_id=rig_id, expt_id=expt_id)

        if "type 1" not in list(cell_index.keys()):
            ctypes = "Not specified"
        else:
            ctypes = f"{cell_index['type 1']:s} and {cell_index['type 2']:s}"

        if "notes" not in info.keys():
            info["notes"] = "  No Notes"

        # populate the NWB subject object from the acq4 metadata

        subject = NWB.file.Subject(
            age=self.ISO8601_age(info["age"]),
            description=info["strain"],
            genotype=info["genotype"],
            sex=info["sex"],
            species=info["species"],
            subject_id=subject_id,
            weight=info["weight"],
            date_of_birth=dob,
        )

        # Populate the NWB amplifier (device) object
        device = NWB.device.Device(
            "MC700B",
            "Current and Voltage clamp amplifier",
            "Axon Instruments (Molecular Devices)",
        )

        # check the acq4 day description field - if it is empty, populate with a default
        if "description" not in info.keys():
            info["description"] = "No description provided"

        # populate the NWB metadata
        self.NWBFile = NWB.NWBFile(
            identifier=str(uuid.uuid4()),  # random uuid for this dataset.
            session_start_time=datetime.datetime.fromtimestamp(info["__timestamp__"], tz=tzlocal()),
            session_id=f"{session_id:s}",
            session_description=info["description"],
            keywords=[
                "mouse",
                "intracellular recording",
                "channelrhodopsins",
                "pathway tracing",
                "auditory cortex",
                "auditory thalamus",
                "medial geniculate",
                "inferior colliculus",
                "brachium of the inferior colliculus",
                "cochlear nucleus",
                "dorsal cochlear nucleus",
                "AAV",
                "synaptic potentials",
                "optogenetics",
            ],
            notes=f"Cell Type: {ctypes:s}\n" + info["notes"],
            protocol="21.123",  # str(path_to_cell.name),
            timestamps_reference_time=datetime.datetime.now(tzlocal()),
            experimenter=[
                "Kasten, Michael R.",
                "Garcia, Michelee",
                "Kline, Amber",
                "Tsukano, Hiroaki",
                "Kato, Hirouki",
                "Manis, Paul B.",
            ],
            lab="Manis Lab",
            institution="UNC Chapel Hill",
            experiment_description="1 R01 RF1NS128873",  # Kato/Manis grant.
            subject=subject,
        )

        self.NWBFile.add_device(device)
        self.sweep_counter = 0  # cumulative counter of traces/sweeps that are stored in this file.

        # In acq4, data is held in subdirectories, one for each "protocol" that was run.
        # Protocols are named according to their 'function' (the type of manipulation/recording
        # that is done), and have an appended number indicating when they are repeated (e.g., _002).
        # Note the the appended numbers may not be sequential - when protocols are stopped prematurely,
        # the protocol is excluded from further analysis, and so is not included here.

        # Now build data structures according to recording mode for each protocol
        print("protocols: ", protocols)

        for protonum, protocol in enumerate(protocols):
            print("Protocol: ", protonum, protocol)

            protocol = protocol.strip()
            path_to_protocol = Path(path_to_cell, protocol)
            print(f"    Datafile: {path_to_protocol!s}\n{" "*14:s}Exists: {path_to_protocol.is_dir()!s}")

            # get the protocol data set
            self.AR.setProtocol(Path(path_to_cell, protocol))
            proto_index = self.AR.readDirIndex(currdir=path_to_protocol)["."]

            recordingmode = proto_index["devices"]["MultiClamp1"]["mode"]
            assert recordingmode in ["IC", "I=0", "CC", "VC"]
            match recordingmode:
                case "IC" | "I=0" | "CC":
                    self.get_CC_data(
                        path_to_cell=Path(path_to_cell),
                        protocol=protocol,
                        info=info,
                        device=device,
                        protonum=protonum,
                    )
                case "VC":
                    self.get_VC_data(
                        path_to_cell=Path(path_to_cell),
                        protocol=protocol,
                        info=info,
                        device=device,
                        protonum=protonum,
                    )
                case _:
                    print(f"Recording mode {recordingmode:s} is not implemented")

        self.NWBFile.generate_new_id()
        with NWB.NWBHDF5IO(str(outfilename), "w") as io:
            io.write(self.NWBFile)
        return outfilename

    def get_one_recording(
        self,
        recording_mode: str,
        series_description: str,
        path_to_cell: Path,
        protocol: str,
        info: dict,
        device: NWB.device,
        acq4_source_name: str = "MultiClamp1.ma",
        electrode: int = 1,
        protonum: int = 0,
    ):
        """Make NWB intracellular electrophysiology tables from a
        single acq4 current clamp or voltage clamp protocol

        Args:
            recording_mode: (str): Either CC or VC for current or voltage clamp
            series_description (str): Name for this  series
            path_to_cell (Path): full path of the data, from which the protocol is pulled
            info (dict): acq4 "info" block provided with this protocol
            device (NWB.device): NWB device structure
            acq4_source_name (str, optional): name of device in acq4. Defaults to "MultiClamp1.ma".
            electrode (int, optional): # of recording electrode. Defaults to 1.
        """

        match recording_mode:
            case "CC":
                stim_name = f"{protocol:s}:Ics{electrode:d}"
                rec_name = f"{protocol:s}:Vcs{electrode:d}"
                opto_name = f"{protocol:s}:OptoSeries"
                stim_units = "amperes"
                rec_units = "volts"

            case "VC":
                stim_name = f"{protocol:s}:Vcs{electrode:d}"
                rec_name = f"{protocol:s}:Ics{electrode:d}"
                opto_name = f"{protocol:s}:OptoSeries"
                stim_units = "volts"
                rec_units = "amperes"

        self.AR.setDataName(acq4_source_name)
        # print(f"\nLooking for a recording from {acq4_source_name:s}")
        dataok = self.AR.getData(check=True)  # just check, don't read raw arrays
        if not dataok:
            return None  # probably wrong channel?
        dataok = self.AR.getData()
        if not dataok:
            print("Error reading data: ", path_to_cell)
            return None

        slicen, cell = self._get_slice_cell(f=path_to_cell)
        for d in self.AR.clampInfo["dirs"]:
            datainfo = self.AR.getDataInfo(Path(d, self.AR.dataname))

        # make the electrode
        elec1 = NWB.icephys.IntracellularElectrode(
            name=rec_name,
            description=f"Sutter 1.5 mm patch as intracellular electrode: {self.day_index['internal']:s}",
            cell_id=f"{slicen:s}_{cell:s}",
            filtering=str(datainfo[1]["ClampState"]["ClampParams"]["PrimarySignalLPF"]),
            device=device,
        )
        # self.NWBFile.add_icephys_electrode(elec1)  # not well documented!

        step_times = [
            np.uint64(self.AR.tstart * self.AR.samp_rate),
            np.uint64(self.AR.tend * self.AR.samp_rate),
        ]

        recordings = []  # keep a list of the recordings

        print("    # sweeps: ", self.AR.cmd_wave.shape[0])
        for sweepno in range(self.AR.cmd_wave.shape[0]):
            
            match recording_mode:
                case "CC":
                    istim = NWB.icephys.CurrentClampStimulusSeries(
                        name=f"{stim_name:s}_{sweepno:d}",  # :sweep{np.uint32(sweepno):d}:protonum{np.uint32(protonum):d}",
                        data=np.array(self.AR.cmd_wave[sweepno, :]),
                        electrode=elec1,
                        gain=datainfo[1]["ClampState"]["extCmdScale"],  # from acq4
                        stimulus_description=str(self.AR.protocol.name),
                        description="Control values are the current injection step times in seconds",
                        control=step_times,
                        control_description=["tstart", "tend"],
                        comments=recording_mode,
                        unit=stim_units,
                        starting_time=info["__timestamp__"],
                        rate=self.AR.sample_rate[0],
                        sweep_number=np.uint32(sweepno),
                    )

                    vdata = NWB.icephys.CurrentClampSeries(
                        name=f"{rec_name:s}_{sweepno:d}",  # sweep{np.uint32(sweepno):d}:protonum{np.uint32(protonum):d}",
                        description=series_description,
                        data=self.AR.data_array[sweepno, :],
                        unit=rec_units,
                        electrode=elec1,
                        gain=datainfo[1]["ClampState"]["primaryGain"],
                        bias_current=datainfo[1]["ClampState"]["ClampParams"]["Holding"],
                        bridge_balance=datainfo[1]["ClampState"]["ClampParams"]["BridgeBalResist"],
                        capacitance_compensation=datainfo[1]["ClampState"]["ClampParams"][
                            "NeutralizationCap"
                        ],
                        stimulus_description="Current Steps",
                        conversion=1.0,
                        timestamps=None,
                        starting_time=info["__timestamp__"],
                        rate=self.AR.sample_rate[0],
                        comments=recording_mode,
                        control=None,
                        control_description=None,
                        sweep_number=np.uint32(sweepno),
                    )

                    IR_sweep_index = self.NWBFile.add_intracellular_recording(
                        electrode=elec1,
                        response=vdata,
                        stimulus=istim,
                        id=sweepno + self.sweep_counter,
                    )

                    odata, osite = self.capture_optical_stimulation(
                        path_to_cell=path_to_cell,
                        protocol=protocol,
                        recording_mode=recording_mode,
                        sweepno=sweepno,
                        electrode=elec1,
                    )
                    if odata is not None:
                        self.NWBFile.add_intracellular_recording(
                            electrode=elec1,
                            response=vdata,
                            stimulus=odata,
                            id=sweepno + self.sweep_counter + 10000,
                        )
                        self.NWBFile.add_ogen_site(osite)
                        ogenseries = NWB.ogen.OptogeneticSeries(
                            name=f"{opto_name:s}_{sweepno:d}",
                            description="Optogenetic stimulation",
                            data=[70.0],
                            site=osite,
                            rate=1.0,
                        )
                        self.NWBFile.add_stimulus(ogenseries)

                case "VC":
                    vstim = NWB.icephys.VoltageClampStimulusSeries(
                        name=f"{stim_name:s}_{sweepno:d}",
                        stimulus_description=str(self.AR.protocol.name),
                        description=series_description,
                        control=step_times,
                        control_description=["tstart", "tend"],
                        comments=recording_mode,
                        data=np.array(self.AR.cmd_wave[sweepno, :]),  # stim_units,
                        starting_time=info["__timestamp__"],
                        rate=self.AR.sample_rate[0],
                        electrode=elec1,
                        gain=datainfo[1]["ClampState"]["extCmdScale"],
                        sweep_number=np.uint32(sweepno),
                    )
                    cparams = datainfo[1]["ClampState"]["ClampParams"]

                    idata = NWB.icephys.VoltageClampSeries(
                        name=f"{rec_name:s}_{sweepno:d}",
                        description=series_description,
                        data=self.AR.data_array[sweepno, :],
                        unit=rec_units,
                        electrode=elec1,
                        gain=datainfo[1]["ClampState"]["primaryGain"],
                        stimulus_description=str(path_to_cell.name),
                        capacitance_fast=cparams["FastCompCap"],
                        capacitance_slow=cparams["SlowCompCap"],
                        resistance_comp_bandwidth=cparams["RsCompBandwidth"],
                        resistance_comp_correction=cparams["RsCompCorrection"],
                        resistance_comp_prediction=0.0,  # not recorded in acq4, we rarely use this.
                        whole_cell_capacitance_comp=cparams["WholeCellCompCap"],
                        whole_cell_series_resistance_comp=cparams["WholeCellCompResist"],
                        resolution=np.NaN,
                        conversion=1.0,
                        timestamps=None,
                        starting_time=info["__timestamp__"],
                        rate=self.AR.sample_rate[0],
                        comments=recording_mode,
                        control=None,
                        control_description=None,
                        sweep_number=np.uint32(sweepno),
                        offset=datainfo[1]["ClampState"]["holding"],
                    )

                    IR_sweep_index = self.NWBFile.add_intracellular_recording(
                        electrode=elec1,
                        stimulus=vstim,
                        response=idata,
                        id=sweepno + self.sweep_counter,
                    )
                    odata, osite = self.capture_optical_stimulation(
                        path_to_cell=path_to_cell,
                        protocol=protocol,
                        recording_mode=recording_mode,
                        sweepno=sweepno,
                        electrode=elec1,
                    )
                    if odata is not None:
                        self.NWBFile.add_intracellular_recording(
                            electrode=elec1,
                            response=idata,
                            stimulus=odata,
                            id=sweepno + self.sweep_counter + 10000,
                        )
                        self.NWBFile.add_ogen_site(osite)
                        ogenseries = NWB.ogen.OptogeneticSeries(
                            name=f"{opto_name:s}_{sweepno:d}",
                            description="Optogenetic stimulation",
                            data=[70.0],
                            site=osite,
                            rate=1.0,
                        )
                        self.NWBFile.add_stimulus(ogenseries)

                case "I=0":
                    print(f"Recording mode {recording_mode:s} not implemented")
                    return None
            recordings.append(IR_sweep_index)

        self.sweep_counter += sweepno + 1
        IR_simultaneous_index = self.NWBFile.add_icephys_simultaneous_recording(
            recordings=recordings,
        )

        # make the sequential group (we only have one "simultaneous" recording)
        IR_sequence_index = self.NWBFile.add_icephys_sequential_recording(
            simultaneous_recordings=[IR_simultaneous_index],
            stimulus_type="square",
        )
        IR_run_index = self.NWBFile.add_icephys_repetition(
            sequential_recordings=[IR_sequence_index]
        )
        self.NWBFile.add_icephys_experimental_condition(repetitions=[IR_run_index])

    def get_CC_data(
        self, path_to_cell: Path, protocol: str, info: dict, device: NWB.device, protonum: int = 0
    ):
        """Check if a protocol is a one of our known
        current-clamp protocols, and if so, put the data into the
        NWB format

        This adds the cc protocol information directly to the self.NWBFile object.

        Args:
            path_to_cell: Path, : full path to the protocol
            info: dict, : acq4 metadata dictionary
            device: NWB.device
        Returns:
            Nothing
        """
        protoname = protocol
        series_description_description = None
        if protoname.startswith("CCIV"):
            series_description_description = "Current Clamp IV series"
        elif protoname.startswith("Ic_LED"):
            series_description_description = (
                "Current Clamp with LED Illumination for optical stimulation"
            )
        elif protoname.startswith("Map_NewBlueLaser_IC"):
            series_description_description = "Current Clamp with Laser scanning photostimulation"
        if series_description_description is None:
            return None
        for ampchannel in self.ampchannels:
            self.get_one_recording(
                recording_mode="CC",
                path_to_cell=path_to_cell,
                series_description=series_description_description,
                protocol=protocol,
                info=info,
                device=device,
                acq4_source_name=f"MultiClamp{ampchannel:d}.ma",
                electrode=ampchannel,
                protonum=protonum,
            )

    def get_VC_data(
        self, path_to_cell: Path, protocol: str, info: dict, device: NWB.device, protonum: int
    ):
        """Check if a protocol is a one of our known
        voltage-clamp protocols, and if so, put the data into the
        NWB format

        This adds the vc protocol information directly to the self.NWBFile object.

         Args:
            path_to_cell: Path, : full path to the protocol
            info: dict, : acq4 metadata dictionary
            device: NWB.device
        Returns:
            Nothing
        """

        vcseries_description = None
        if protocol.startswith("VCIV"):
            vcseries_description = "Voltage Clamp series"
        elif protocol.startswith(("Vc_LED", "VC_LED")):
            vcseries_description = "Voltage Clamp with LED Illumination for optical stimulation"
        elif protocol.startswith("Map_NewBlueLaser_VC"):
            vcseries_description = "Voltage Clamp with Laser scanning photostimulation"
        if vcseries_description is None:
            return None

        for ampchannel in self.ampchannels:
            self.get_one_recording(
                recording_mode="VC",
                series_description=vcseries_description,
                path_to_cell=path_to_cell,
                protocol=protocol,
                info=info,
                device=device,
                acq4_source_name=f"MultiClamp{ampchannel:d}.ma",
                electrode=ampchannel,
                protonum=protonum,
            )

    def capture_optical_stimulation(
        self,
        path_to_cell: Path,
        protocol: str,
        recording_mode: str = "VC",
        sweepno: int = 0,
        electrode: str = "None",
    ):
        # capture optical stimulation information as well

        # The optical stimulation data is held as a VC stimulus series
        # This should be changed in the future
        # The "control" element holds a 2-d array of sites for laser scanning
        # (sorry, the whole array is repeated in every trial).
        # nwb wants the control data to be in integer format, so we scale to microns

        odata = None
        osite = None
        #  laser scanning:
        if self.AR.getLaserBlueCommand():
            if self.laser_device is None:
                self.laser_device = self.NWBFile.create_device(name="Oxxius Laser", description="450 nm 70 mW single-photon Laser")
            scinfo = self.AR.getScannerPositions()
            sites = [self.AR.scanner_info[spot]["pos"] for spot in self.AR.scanner_info.keys()]
            for i, site in enumerate(sites):
                sites[i] = (
                    np.int32(sites[i][0] * 1e6),
                    np.int32(sites[i][1] * 1e6),
                )
                # sites = map(sites, lambda x: (np.int32(x[0]*1e6), np.int32(x[1]*1e6)))
            # sites = str(sites)  # must convert to string to pass nwb tests.
            location_column = VectorData(
                name="sites",
                data=sites,
                description="x, y positions in microns for laser scanning photostimulation",
            )
            control_description = "Laser scanning (spot) photostimulation positions (in microns)"
            light_source = "450 nm 70 mW Laser"

            if recording_mode == "VC":
                odata = NWB.icephys.VoltageClampStimulusSeries(
                    name=f"{protocol:s}:opto_{sweepno:d}",
                    description=f"Optical stimulation control waveform and xy positions for {light_source:s}",
                    data=self.AR.LaserBlue_pCell,
                    rate=self.AR.LaserBlue_sample_rate[0],
                    electrode=electrode,
                    gain=1.0,
                    stimulus_description=str(path_to_cell.name),
                    conversion=1.0,
                    # control=sites,
                    control_description=control_description,
                    comments=f"{self.AR.scanner_spotsize:e}",
                    unit="volts",
                )

            if recording_mode == "CC":
                odata = NWB.icephys.CurrentClampStimulusSeries(
                    name=f"{protocol:s}:opto_{sweepno:d}",
                    description=f"Optical stimulation control waveform and xy positions for {light_source:s}",
                    data=self.AR.LaserBlue_pCell,
                    rate=self.AR.LaserBlue_sample_rate[0],
                    electrode=electrode,
                    gain=1.0,
                    stimulus_description=str(path_to_cell.name),
                    conversion=1.0,
                    # control=sites,
                    control_description=control_description,
                    comments=f"{self.AR.scanner_spotsize:e}",
                    unit="amperes",
                )
            osite = NWB.ogen.OptogeneticStimulusSite(
                name=f"{protocol:s}:sweep={sweepno:d}",
                device=self.laser_device,
                description="Laser Scanning Photostimulation site",
                excitation_lambda=450.0,
                location=str(sites[sweepno]),
            )
        # widefiled LED:
        if self.AR.getLEDCommand():
            light_source = "470 nm Thorlabs LED"
            if self.LED_device is None:
                self.LED_device = self.NWBFile.create_device(name=f"{light_source:s}", description="Widefield LED")
            sites = [[0,0]]*np.array(self.AR.LED_Raw).shape[0]  # also possible to get fov from acq4... 
            spotsize = 1e-4  # 100 microns.
            control_description = "Widefield LED through objective"

            if recording_mode == "VC":
                odata = NWB.icephys.VoltageClampStimulusSeries(
                    name=f"{protocol:s}:opto_{sweepno:d}",
                    description=f"Optical stimulation waveform for {light_source:s}",
                    data=np.array(self.AR.LED_Raw)[:, sweepno],  # assuming is constant...
                    rate=self.AR.LED_sample_rate[0],
                    electrode=electrode,
                    gain=1.0,
                    stimulus_description=str(path_to_cell.name),
                    conversion=1.0,
                    # control=sites,
                    control_description=control_description,
                    comments=f"Spotsize: {spotsize:e}",
                    unit="volts",
                )
            if recording_mode == "CC":
                odata = NWB.icephys.CurrentClampStimulusSeries(
                    name=f"{protocol:s}:opto_{sweepno:d}",
                    description=f"Optical stimulation waveform",
                    data=np.array(self.AR.LED_Raw)[:, sweepno],  # assuming is constant...
                    rate=self.AR.LED_sample_rate[0],
                    electrode=electrode,
                    gain=1.0,
                    stimulus_description=str(path_to_cell.name),
                    conversion=1.0,
                    # control=sites,
                    control_description=control_description,
                    comments=f"Spotsize: {spotsize:e}",
                    unit="amperes",
                )

            osite = NWB.ogen.OptogeneticStimulusSite(
                name=f"{protocol:s}:optosite_{sweepno:d}",
                device=self.LED_device,
                description="Widefield optical stimulation",
                excitation_lambda=470.0,
                location=str(sites[sweepno]),
            )
        return odata, osite


# Read the data back in
def validate(testpath, nwbfile):
    with NWBHDF5IO(testpath, "r") as io:
        infile = io.read()

        # assert intracellular_recordings
        assert np.all(
            infile.intracellular_recordings.id[:] == nwbfile.intracellular_recordings.id[:]
        )

        # Assert that the ids and the VectorData, VectorIndex, and table target of the
        # recordings column of the Sweeps table are correct
        assert np.all(
            infile.icephys_simultaneous_recordings.id[:]
            == nwbfile.icephys_simultaneous_recordings.id[:]
        )
        assert np.all(
            infile.icephys_simultaneous_recordings["recordings"].target.data[:]
            == nwbfile.icephys_simultaneous_recordings["recordings"].target.data[:]
        )
        assert np.all(
            infile.icephys_simultaneous_recordings["recordings"].data[:]
            == nwbfile.icephys_simultaneous_recordings["recordings"].data[:]
        )
        assert (
            infile.icephys_simultaneous_recordings["recordings"].target.table.name
            == nwbfile.icephys_simultaneous_recordings["recordings"].target.table.name
        )

        # Assert that the ids and the VectorData, VectorIndex, and table target of the simultaneous
        #  recordings column of the SweepSequences table are correct
        assert np.all(
            infile.icephys_sequential_recordings.id[:]
            == nwbfile.icephys_sequential_recordings.id[:]
        )
        assert np.all(
            infile.icephys_sequential_recordings["simultaneous_recordings"].target.data[:]
            == nwbfile.icephys_sequential_recordings["simultaneous_recordings"].target.data[:]
        )
        assert np.all(
            infile.icephys_sequential_recordings["simultaneous_recordings"].data[:]
            == nwbfile.icephys_sequential_recordings["simultaneous_recordings"].data[:]
        )
        assert (
            infile.icephys_sequential_recordings["simultaneous_recordings"].target.table.name
            == nwbfile.icephys_sequential_recordings["simultaneous_recordings"].target.table.name
        )

        # Assert that the ids and the VectorData, VectorIndex, and table target of the
        # sequential_recordings column of the Repetitions table are correct
        assert np.all(infile.icephys_repetitions.id[:] == nwbfile.icephys_repetitions.id[:])
        assert np.all(
            infile.icephys_repetitions["sequential_recordings"].target.data[:]
            == nwbfile.icephys_repetitions["sequential_recordings"].target.data[:]
        )
        assert np.all(
            infile.icephys_repetitions["sequential_recordings"].data[:]
            == nwbfile.icephys_repetitions["sequential_recordings"].data[:]
        )
        assert (
            infile.icephys_repetitions["sequential_recordings"].target.table.name
            == nwbfile.icephys_repetitions["sequential_recordings"].target.table.name
        )

        # Assert that the ids and the VectorData, VectorIndex, and table target of the
        # repetitions column of the Conditions table are correct
        assert np.all(
            infile.icephys_experimental_conditions.id[:]
            == nwbfile.icephys_experimental_conditions.id[:]
        )
        assert np.all(
            infile.icephys_experimental_conditions["repetitions"].target.data[:]
            == nwbfile.icephys_experimental_conditions["repetitions"].target.data[:]
        )
        assert np.all(
            infile.icephys_experimental_conditions["repetitions"].data[:]
            == nwbfile.icephys_experimental_conditions["repetitions"].data[:]
        )
        assert (
            infile.icephys_experimental_conditions["repetitions"].target.table.name
            == nwbfile.icephys_experimental_conditions["repetitions"].target.table.name
        )
        assert np.all(
            infile.icephys_experimental_conditions["tag"][:]
            == nwbfile.icephys_experimental_conditions["tag"][:]
        )


def ConvertFile(
    experiment_name: str,
    filename: Union[str, Path],
    protocols: list,
    appendmode: bool = False,
    outputpath: str = "test_data",
    device="MultiClamp1.ma",
    mode="CC",
):
    # experiment_name = None
    A2N = ACQ4toNWB(outputpath)
    A2N.set_amplifier_name(device)
    NWBFile = A2N.acq4tonwb(
        experiment_name=experiment_name,
        path_to_cell=filename,
        protocols=protocols,
        appendmode=appendmode,
        recordingmode=mode,
    )
    print("NWB file written:  ", NWBFile)
    return NWBFile


def main():
    """Read the pandas database for this set of experiments.
    The pandas database is generated by dataSummary.py (from github.com/pbmanis/ephys/util).
    This database is a complete record of the experiments, not filtered by any criteria.
    Each row of the database corresponds to one day's recordings - usually one subject
    The column "data_complete" holds the list of protocols that were completed and are
    potentially suitable for further analysis.
    """
    import pandas as pd

    HK_db = pd.read_pickle(
        # "/Users/pbmanis/Desktop/Python/HK_Collab/datasets/HK_fastcortex_cortex.pkl"
        # "/Users/pbmanis/Desktop/Python/HK_Collab/datasets/HK_fastcortex_DCN.pkl"
        "/Users/pbmanis/Desktop/Python/HK_Collab/datasets/Thalamocortical/HK_TC_datasummary-2024.02.10.pkl"
    )
    experiment_name = "Thalamocortical"
    print("Converting acq4 to NWB from database: ", HK_db.cell_id)

    def select(protocols: list):
        ok_prots = []
        for protocol in protocols:
            if (
                protocol.startswith("CCIV")
                or protocol.startswith("Map")
                or protocol.startswith("Vc_LED")
                or protocol.startswith("Ic_LED")
                or protocol.startswith("VC_LED")
            ):
                ok_prots.append(protocol)
        return ok_prots

    NWBFiles = []
    files = []

    def cvt(row, experiment_name, rigs: Union[Tuple, List] = None, row_num: int = 0):
        """Do the file conversion for each complete data set in a given database row

        Args:
            row (Pandas row (series)): A row containing information about the day
        """

        path = row.data_directory
        day = row.date
        if rigs is not None:  # restrict by rig
            rig_ok = False
            for rig in rigs:
                if day.startswith(rig):  # only do data from a given rig.
                    rig_ok = True
            if not rig_ok:
                return  # did not match rigs
        slice = row.slice_slice
        cell = row.cell_cell
        protocols = row.data_complete.split(",")
        protocols = [p.strip() for p in protocols if p.strip() not in [" ", "", None]]
        protocols = select(protocols)
        if len(protocols) == 0:
            return
        print("\nCell id: ", row.cell_id,  " index: ", row_num)
        print("     Protocols: ", protocols)
        file = Path(path, day, slice, cell)
        files.append(file)
        NWBFile = ConvertFile(experiment_name=experiment_name, filename=file, protocols=protocols)
        NWBFiles.append(NWBFile)

    for index in HK_db.index:
        if index <= 64:  # just for a test...
            continue
        print(HK_db.iloc[index].cell_id)
        cvt(HK_db.iloc[index], experiment_name, rigs=["Rig2", "Rig4"], row_num=index)

    # HK_db = HK_db.apply(cvt, axis=1, experiment_name=experiment_name, rig="Rig2")
    # check the files:
    print("\n", "=" * 80)
    # ConvertFile.validate(NWBFile[0])

def check_conversion(nwbf):
    try:
        results = list(inspect_nwbfile(nwbfile_path=nwbf))
        if len(results) == 0:
            print("    Conversion OK")
        else:
            print(f"    Error in conversion detected: for file {i:d}: {nwbf!s} ")
            print(results)
            exit()
    except:
        print(f"    Error in conversion detected: for file  {nwbf!s} ")
        print(results)
        exit()

def check_conversions():
    print("Checking conversions")
    NWBFiles = []
    filepath = "/Users/pbmanis/Desktop/Python/Dandi/test_data"
    for f in Path(filepath).rglob("*.nwb"):
        NWBFiles.append(f)
    # print("NWBFiles: ", NWBFiles)
    for i, nwbf in enumerate(NWBFiles):
        if nwbf is None:
            continue
        check_conversion(nwbf)


if __name__ == "__main__":
    main()  # convert a bunch of files
    check_conversions()
    # check_conversion("/Users/pbmanis/Desktop/Python/Dandi/test_data/Thalamocortical~Rig2~2023.03.31~S01C00.nwb")  # check the conversion
