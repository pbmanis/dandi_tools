"""
Read and convert acq4 data to NWB format that is acceptable for the Dandi repository
Specifically designed for the auditory cortex experiments with Kato lab.
A single file can be converted by calling ConvertFile(filename)
The conversion is checked after generating the outputfile with nwbinspector.

NWFFile expects:
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
        lab_meta_data (list or tuple) :: an extension that contains lab-specific meta-data
        icephys_electrodes (list or tuple) :: IntracellularElectrodes that belong to this NWBFile.
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
        intracellular_recordings (IntracellularRecordingsTable) :: the IntracellularRecordingsTable table that belongs to this NWBFile
        icephys_simultaneous_recordings (SimultaneousRecordingsTable) :: the SimultaneousRecordingsTable table that belongs to this NWBFile
        icephys_sequential_recordings (SequentialRecordingsTable) :: the SequentialRecordingsTable table that belongs to this NWBFile
        icephys_repetitions (RepetitionsTable) :: the RepetitionsTable table that belongs to this NWBFile
        icephys_experimental_conditions (ExperimentalConditionsTable) :: the ExperimentalConditionsTable table that belongs to this NWBFile

        icephys_filtering (str) :: [DEPRECATED] Use IntracellularElectrode.filtering instead. Description of filtering used.
        ic_electrodes (list or tuple) :: DEPRECATED use icephys_electrodes parameter instead. IntracellularElectrodes that belong to this NWBFile

Support::

    NIH grants:
    DC RF1 NS128873 (Kato, Manis, MPI, 2022-) Cortical circuits for the integration of parallel short-latency auditory pathways.
    DC R01 DC019053 (Manis, 2020-2025) Cellular mechanisms of auditory information processing. 

Copyright 2022-2023 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 

"""
import datetime
from dataclasses import dataclass, field
from pathlib import Path
from typing import Union, Literal
import uuid
import ephys.datareaders as DR
import numpy as np
import pynwb as NWB
from dateutil.tz import tzlocal
from nwbinspector import inspect_all
from pynwb import NWBHDF5IO


def def_lister():
    return []


@dataclass
class ExperimentInfo:
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


class ACQ4toNWB:
    def __init__(self, out_file_path: Union[str, Path, None] = None):
        """Convert a file from ACQ4 (www.acq4.org) format to NWB format (www.nwb.org)

        Args:
            out_file_path (Union[str, Path, None], optional): The path to the output files. Defaults to None.
        """
        self.AR = DR.acq4_reader.acq4_reader()  # get acq4 reader

        self.out_file_path = out_file_path
        # All of the data here is single electrode, and "MultiClamp1.ma" is the
        # name of the data file holding the ephys data. The acq4 files are held in "metaarray" format, but
        # can be read as hdf5 files.
        self.set_data_name("MultiClamp1.ma")
        self.manager = NWB.get_manager()

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
        Convert a name of this format:
        f = Path('/Volumes/Pegasus/all/the/leading/pathnames/2017.05.01_000/slice_000/cell_000')
        To:
        2017.05.01~S0C0

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
        # protocol = f.stem
        # cellp = f.stem
        cell = fp[-1]
        slicenp = fp[-2]
        day = fp[-3]
        if expt_id == None:
            expt_id = ""
        if rig_id == None:
            rig_id = ""
        foname = str(
            Path(expt_id, rig_id, day[:10], "S" + slicenp[-2:] + "C" + cell[-2:])
        ).replace("/", "~")
        return foname

    def set_data_name(self, datatype: str):
        """Convenience function to change the amplifier data name

        Args:
            datatype (str): name of the amplifier
        """
        self.AR.setDataName("MultiClamp1.ma")

    def ISO8601_age(self, agestr):
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
        if "P" not in agestr:
            agestr = "P" + agestr
        if "D" not in agestr:
            agestr = agestr + "D"
        if agestr == "PD":
            agestr = "P9999D"  # no age specified
        return agestr

    def find_name_in_path(
        self, sourcepath: Union[str, Path], name: Union[str, None] = None
    ):
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
        """Convert one protocol directory to an NWB file

        Args:
            path_to_cell (string or Path): Full path and name of the protocol directory
            outfilename (Union[Path, str, None], optional): NWB output filename. Defaults to None.
            recordingmode (Literal[&quot;IC&quot;, &quot;CC&quot;, &quot;VC&quot;, &quot;I, optional): _description_. Defaults to 0"]="CC".

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

        info = self.AR.readDirIndex(currdir=path_to_cell.parent.parent)["."]
        slice_index = self.AR.readDirIndex(currdir=path_to_cell.parent)["."]
        cell_index = self.AR.readDirIndex(currdir=path_to_cell)["."]

        # data in self.AR.data_array, self.AR.time_base, stimuli in self.AR.cmd_wave
        data_date = datetime.date.fromtimestamp(info["__timestamp__"])

        # Clean up metadata
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

        # make NWB subject

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

        device = NWB.device.Device(
            "MC700B",
            "Current and Voltage clamp amplifier",
            "Axon Instruments (Molecular Devices)",
        )
        self.NWBFile = NWB.NWBFile(
            identifier=str(uuid.uuid4()),  # random uuid for this dataset.
            session_start_time=datetime.datetime.fromtimestamp(
                info["__timestamp__"], tz=tzlocal()
            ),
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
            experimenter=["Kasten, Michael R.", "Manis, Paul B."],
            lab="Manis Lab",
            institution="UNC Chapel Hill",
            experiment_description="1 R01 RF1NS128873",
            subject=subject,
        )

        self.NWBFile.add_device(device)

        # build data structures according to recording mode for each protocol
        for protocol in protocols:
            protocol = protocol.strip()
            path_to_protocol = Path(path_to_cell, protocol)
            print("datafile: ", path_to_protocol, path_to_protocol.is_dir())
            self.AR.setProtocol(Path(path_to_cell, protocol))
            dataok = self.AR.getData()

            if not dataok:
                print("Error reading data: ", path_to_protocol)
                return None

            proto_index = self.AR.readDirIndex(currdir=path_to_protocol)["."]

            recordingmode = proto_index["devices"]["MultiClamp1"]["mode"]
            assert recordingmode in ["IC", "I=0", "CC", "VC"]
            if recordingmode in ["IC", "I=0", "CC"]:
                self.get_CC_data(
                    path_to_cell=path_to_cell,
                    protocol=protocol,
                    info=info,
                    device=device,
                )
            elif recordingmode == "VC":
                self.get_VC_data(
                    path_to_cell=path_to_cell,
                    protocol=protocol,
                    info=info,
                    device=device,
                )
            else:
                print(f"Recording mode {recordingmode:s} is not implemented yet")
                return None

        with NWB.NWBHDF5IO(str(outfilename), "w", manager=self.manager) as io:
            io.write(self.NWBFile)
        return outfilename

    def get_one_CC_recording(
        self,
        ccseries: str,
        path_to_cell: Path,
        protocol: str,
        info: dict,
        device: NWB.device,
        acq4_source_name: str = "MultiClamp1.ma",
        electrode: int = 1,
    ):
        """Make a current clamp NWB series from a protocol

        Args:
            ccseries (str): Name for this CC series
            path_to_cell (Path): full path of the data, from which the protocol is pulled
            info (dict): acq4 "info" block provided with this protocol
            device (NWB.device): NWB device structure
            acq4_source_name (str, optional): name of device in acq4. Defaults to "MultiClamp1.ma".
            electrode (int, optional): # of recording electrode. Defaults to 1.
        """
        self.AR.setDataName(acq4_source_name)
        print(f"\nLooking for a recordings from {acq4_source_name:s}")
        dataok = self.AR.getData(check=True)  # just check, don't read raw arrays
        if not dataok:
            print(f"   Source data not found, exiting")
            exit()
        self.AR.getData()
        slicen, cell = self._get_slice_cell(f=path_to_cell)
        s = []
        for d in self.AR.clampInfo["dirs"]:
            datainfo = self.AR.getDataInfo(Path(d, self.AR.dataname))

        elec1 = NWB.icephys.IntracellularElectrode(
            name=f"{protocol:s}:elec{electrode:d}",
            description="Sutter 1.5 mm patch as intracellular electrode",
            cell_id=f"{slicen:s}_{cell:s}",
            filtering=str(datainfo[1]["ClampState"]["ClampParams"]["PrimarySignalLPF"]),
            device=device,
        )
        self.NWBFile.add_icephys_electrode(elec1)  # not well documented!

        step_times = [
            np.uint64(self.AR.tstart * self.AR.samp_rate),
            np.uint64(self.AR.tend * self.AR.samp_rate),
        ]
        istim = NWB.icephys.CurrentClampStimulusSeries(
            name=f"{protocol:s}:Ics{electrode:d}",
            stimulus_description=str(self.AR.protocol.name),
            description="Control values are the current injection step times in seconds",
            control=step_times,
            control_description=["tstart", "tend"],
            comments="CC",
            data=np.array(self.AR.cmd_wave).T,
            unit="amperes",
            starting_time=info["__timestamp__"],
            rate=self.AR.sample_rate[0],
            electrode=elec1,
            gain=datainfo[1]["ClampState"]["extCmdScale"],
            sweep_number=np.uint32(1),
        )

        vdata = NWB.icephys.CurrentClampSeries(
            name=f"{protocol:s}:Vcs1",
            description=ccseries,
            data=self.AR.data_array.T,
            unit="volts",
            electrode=elec1,
            gain=datainfo[1]["ClampState"]["primaryGain"],
            bias_current=datainfo[1]["ClampState"]["ClampParams"]["Holding"],
            bridge_balance=datainfo[1]["ClampState"]["ClampParams"]["BridgeBalResist"],
            capacitance_compensation=datainfo[1]["ClampState"]["ClampParams"][
                "NeutralizationCap"
            ],
            stimulus_description="Current Steps",
            resolution=np.NaN,
            conversion=1.0,
            timestamps=None,
            starting_time=info["__timestamp__"],
            rate=self.AR.sample_rate[0],
            comments="CC",
            control=None,
            control_description=None,
            sweep_number=np.uint(1),
        )
        # capture optical stimulation information as well

        scinfo = self.AR.getScannerPositions(dataname="Laser-Blue-raw.ma")
        lblue = self.AR.getLaserBlueCommand()
        odata = None
        # print("scinfo, lblue: ", scinfo, lblue)
        if scinfo and lblue:  # we have optical stimulation in the mix
            if "led" in protocol.lower():
                sites = [(0, 0)]
                control_description = "Widefield LED through objective"
                light_source = "470 nm LED"
            elif "laser" in protocol.lower():
                sites = [
                    self.AR.scannerinfo[spot]["pos"]
                    for spot in self.AR.scannerinfo.keys()
                ]
                for i, site in enumerate(sites):
                    sites[i] = (
                        np.int32(sites[i][0] * 1e6),
                        np.int32(sites[i][1] * 1e6),
                    )
                    # sites = map(sites, lambda x: (np.int32(x[0]*1e6), np.int32(x[1]*1e6)))
                control_description = (
                    "Laser scanning (spot) photostimulation positions (in microns)"
                )
                light_source = "450 nm 70 mW Laser"
            else:
                raise ValueError(
                    f"Protocol must include led or laser. got {protocol:s}"
                )

            # The optical stimulation data is held as a VC stimulus series
            # This should be changed in the future to properly link the
            # vcseries above with the optical stimulation using a table.
            # For now, assume they are parallel.
            # The "control" element holds a 2-d array of sites.
            # nwb wants the control data to be in integer format, so we scale to microns

            odata = NWB.icephys.CurrentClampStimulusSeries(
                name=f"{protocol:s}:ChR2 Stimulation {electrode:d}",
                description=f"Optical stimulation control waveform and xy positions for {light_source:s}",
                data=self.AR.LaserBlue_pCell.T,
                rate=self.AR.LBR_sample_rate[0],
                electrode=elec1,
                gain=1.0,
                stimulus_description=str(path_to_cell.name),
                conversion=1.0,
                control=sites,
                control_description=control_description,
                comments=f"{self.AR.spotsize:e}",
                unit="volts",
            )
        self.NWBFile.add_acquisition(istim)
        self.NWBFile.add_acquisition(vdata)
        if odata is not None:
            self.NWBFile.add_acquisition(odata)

    def get_CC_data(
        self, path_to_cell: Path, protocol: str, info: dict, device: NWB.device
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
        ccseries_description = None
        if protoname.startswith("CCIV"):
            ccseries_description = "Current Clamp IV series"
        elif protoname.startswith("Ic_LED"):
            ccseries_description = (
                "Current Clamp with LED Illumination for optical stimulation"
            )
        elif protoname.startswith("Map_NewBlueLaser_IC"):
            ccseries_description = "Current Clamp with Laser scanning photostimulation"
        if ccseries_description is None:
            return None
        for ampchannel in [1]:
            self.get_one_CC_recording(
                path_to_cell=path_to_cell,
                ccseries=ccseries_description,
                protocol=protocol,
                info=info,
                device=device,
                acq4_source_name=f"MultiClamp{ampchannel:d}.ma",
                electrode=ampchannel,
            )

    def get_one_VC_recording(
        self,
        vcseries: str,
        path_to_cell: Path,
        protocol: str,
        info: dict,
        device: NWB.device,
        acq4_source_name: str = "MultiClamp1.ma",
        electrode: int = 1,
    ):
        """
        Check if a protocol is a one of our known
        voltage-clamp protocols, and if so, put the data into the
        NWB format

        This adds the vc protocol information directly to the self.NWBFile object.

        Args:
            vcseries (str): Name for this VC series
            path_to_cell (Path): full path of the data, from which the protocol is pulled
            info (dict): acq4 "info" block provided with this protocol
            device (NWB.device): NWB device structure
            acq4_source_name (str, optional): name of device in acq4. Defaults to "MultiClamp1.ma".
            electrode (int, optional): # of recording electrode. Defaults to 1.
        """

        slicen, cell = self._get_slice_cell(f=path_to_cell)
        path_to_protocol = Path(path_to_cell, protocol)
        s = []
        for d in self.AR.clampInfo["dirs"]:
            datainfo = self.AR.getDataInfo(Path(d, self.AR.dataname))

        elec = NWB.icephys.IntracellularElectrode(
            name=f"{protocol:s}:elec{electrode:d}",
            description="Sutter 1.5 mm patch as intracellular electrode",
            cell_id=f"{slicen:s}_{cell:s}",
            filtering=str(datainfo[1]["ClampState"]["ClampParams"]["PrimarySignalLPF"]),
            device=device,
        )
        self.NWBFile.add_icephys_electrode(elec)  # not well documented!

        step_times = [
            np.uint64(self.AR.tstart * self.AR.samp_rate),
            np.uint64(self.AR.tend * self.AR.samp_rate),
        ]
        vstim = NWB.icephys.VoltageClampStimulusSeries(
            name=f"{protocol:s}:Vcs{electrode:d}",
            stimulus_description=str(self.AR.protocol.name),
            description=vcseries,
            control=step_times,
            control_description=["tstart", "tend"],
            comments="VC",
            data=np.array(self.AR.cmd_wave).T,
            unit="volts",
            starting_time=info["__timestamp__"],
            rate=self.AR.sample_rate[0],
            electrode=elec,
            gain=datainfo[1]["ClampState"]["extCmdScale"],
            sweep_number=np.uint32(1),
        )
        cparams = datainfo[1]["ClampState"]["ClampParams"]

        idata = NWB.icephys.VoltageClampSeries(
            name=f"{protocol:s}:Ics{electrode:d}",
            description=vcseries,
            data=self.AR.data_array.T,
            unit="amperes",
            electrode=elec,
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
            comments="VC",
            control=None,
            control_description=None,
            sweep_number=np.uint64(1),
            offset=datainfo[1]["ClampState"]["holding"],
        )

        # capture optical stimulation information as well
        scinfo = self.AR.getScannerPositions(dataname="Laser-Blue-raw.ma")
        lblue = self.AR.getLaserBlueCommand()
        odata = None
        # print("scinfo, lblue: ", scinfo, lblue)
        if scinfo and lblue:  # we have optical stimulation in the mix
            if "led" in protocol.lower():
                sites = [(0, 0)]
                control_description = "Widefield LED through objective"
                light_source = "470 nm LED"
            elif "laser" in protocol.lower():
                sites = [
                    self.AR.scannerinfo[spot]["pos"]
                    for spot in self.AR.scannerinfo.keys()
                ]
                for i, site in enumerate(sites):
                    sites[i] = (
                        np.int32(sites[i][0] * 1e6),
                        np.int32(sites[i][1] * 1e6),
                    )
                    # sites = map(sites, lambda x: (np.int32(x[0]*1e6), np.int32(x[1]*1e6)))
                control_description = (
                    "Laser scanning (spot) photostimulation positions (in microns)"
                )
                light_source = "450 nm 70 mW Laser"
            else:
                raise ValueError("Protocol must include led or laser")

            # The optical stimulation data is held as a VC stimulus series
            # This should be changed in the future to properly link the
            # vcseries above with the optical stimulation using a table.
            # For now, assume they are parallel.
            # The "control" element holds a 2-d array of sites.
            # nwb wants the control data to be in integer format, so we scale to microns

            odata = NWB.icephys.VoltageClampStimulusSeries(
                name=f"{protocol:s}:ChR2 Stimulation {electrode:d}",
                description=f"Optical stimulation control waveform and xy positions for {light_source:s}",
                data=self.AR.LaserBlue_pCell.T,
                rate=self.AR.LBR_sample_rate[0],
                electrode=elec,
                gain=1.0,
                stimulus_description=str(path_to_cell.name),
                conversion=1.0,
                control=sites,
                control_description=control_description,
                comments=f"{self.AR.spotsize:e}",
                unit="volts",
            )

        self.NWBFile.add_acquisition(vstim)
        self.NWBFile.add_acquisition(idata)
        if odata is not None:
            self.NWBFile.add_acquisition(odata)

    def get_VC_data(
        self, path_to_cell: Path, protocol: str, info: dict, device: NWB.device
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
        elif protocol.startswith("Vc_LED"):
            vcseries_description = (
                "Voltage Clamp with LED Illumination for optical stimulation"
            )
        elif protocol.startswith("Map_NewBlueLaser_VC"):
            vcseries_description = "Voltage Clamp with Laser scanning photostimulation"
        if vcseries_description is None:
            return None

        for ampchannel in [1, 2]:
            self.get_one_VC_recording(
                vcseries=vcseries_description,
                path_to_cell=path_to_cell,
                protocol=protocol,
                info=info,
                device=device,
                acq4_source_name=f"MultiClamp{ampchannel:d}.ma",
                electrode=ampchannel,
            )


def ConvertFile(
    experiment_name: str,
    filename: Union[str, Path],
    protocols: str,
    appendmode: bool = False,
    outputpath: str = "test_data",
    device="MultiClamp1.ma",
    mode="CC",
):
    # experiment_name = None
    A2N = ACQ4toNWB(outputpath)
    A2N.set_data_name(device)
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
       "/Users/pbmanis/Desktop/Python/HK_Collab/datasets/HK_fastcortex_cortex.pkl"
        # "/Users/pbmanis/Desktop/Python/HK_Collab/datasets/HK_fastcortex_DCN.pkl"
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
            ):
                ok_prots.append(protocol)
        return ok_prots

    NWBFiles = []

    def cvt(row, experiment_name):
        """Do the file conversion for each complete data set in a given database row

        Args:
            row (Pandas row (series)): A row containing information about the day
        """
        path = row.data_directory
        day = row.date
        slice = row.slice_slice
        cell = row.cell_cell
        protocols = row.data_complete.split(",")
        protocols = [p.strip() for p in protocols if p.strip() not in [" ", "", None]]
        protocols = select(protocols)
        if len(protocols) == 0:
            return
        print("\n   Cell id: ", row.cell_id)
        print("   protocols", protocols)
        file = Path(path, day, slice, cell)
        NWBFile = ConvertFile(
            experiment_name=experiment_name, filename=file, protocols=protocols
        )
        NWBFiles.append(NWBFile)

    HK_db = HK_db.apply(cvt, axis=1, experiment_name=experiment_name)
    # check the files:
    print("\n", "=" * 80)
    for nwbf in NWBFiles:
        if nwbf is None:
            continue
        results = list(inspect_all(path=nwbf))
        if len(results) == 0:
            print("    Conversion OK")
        else:
            print("    Error in conversion detected: ")
            print(results)


if __name__ == "__main__":
    main()  # convert a bunch of files
    # f1 = Path('/Volumes/Pegasus_002/ManisLab_Data3/Kasten_Michael/HK_Collab/Thalamocortical/Rig4/2022.11.29_000/slice_001/cell_002/Vc_LED_stim_onebig_000')
    # f1 = Path(
    #     "/Volumes/Pegasus_002/ManisLab_Data3/Kasten_Michael/HK_Collab/Thalamocortical/Rig4/2021.05.04_000/slice_000/cell_000/Map_NewBlueLaser_VC_increase_1ms_HK_001"
    # )
    # # f2 = Path('/Volumes/Pegasus_002/ManisLab_Data3/Kasten_Michael/HK_Collab/Thalamocortical/Rig4/2022.10.10_000/slice_001/cell_000/Vc_LED_stim_wut_000')
    # f2 = Path('/Volumes/Pegasus_002/ManisLab_Data3/Kasten_Michael/HK_Collab/Thalamocortical/Rig4/2022.10.10_000/slice_001/cell_000/Ic_LED_stim_000')
    # # f4 = Path('/Volumes/Pegasus_002/ManisLab_Data3/Kasten_Michael/HK_Collab/Thalamocortical/Rig4/2022.10.10_000/slice_001/cell_000/CCIV_1nA_max_1s_pulse_posonly_000')
    # f3 = Path('/Volumes/Pegasus/ManisLab_Data3/Kasten_Michael/2017.02.23_000/slice_000/cell_000/Map_NewBlueLaser_VC_single_MAX_005')
    # f4 = Path('/Volumes/Pegasus_002/ManisLab_Data3/Kasten_Michael/NF107Ai32_Het/2019.04.16_000/slice_002/cell_000/VCIV_simple_001')
    # convert files from a list
    # for f in [f1]:
    #     ConvertFile(filename=f, mode="VC")
