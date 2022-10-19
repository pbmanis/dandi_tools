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

Copyright 2022 Paul B. Manis
Distributed under MIT/X11 license. See license.txt for more infomation. 

"""
import datetime
from dataclasses import dataclass
from pathlib import Path, PurePath
from typing import Union, Literal

import ephys.ephysanalysis as EP
import numpy as np
import pynwb as NWB
from dateutil.tz import tzlocal
from ephys.ephysanalysis import MatdatacRead as MDR
from nwbinspector import inspect_all, inspect_nwb


@dataclass
class ExperimentInfo:
    description: str = ""
    protocol: str = ""
    time: str = ""
    experimenters: str = ""
    lab: str = ""
    institution: str = "UNC Chapel Hill"
    experiment_description: str = ""
    sessionid: str = "0"
    notes: str = ""
    subject: str = ""

class ACQ4toNWB:
    def __init__(self, out_file_path: Union[str, Path, None] = None):
        self.AR = EP.acq4read.Acq4Read() # get acq4 reader

        self.out_file_path = out_file_path
        self.set_data_name("MultiClamp1.ma")


    def _get_slice_cell(self, f: Union[str, Path, None] = None):
        f = Path(f)
        protocol = f.stem
        cellp = f.parent
        cell = cellp.stem
        slicenp = cellp.parent
        slicen = slicenp.stem
        day = slicenp.parent
        return slicen, cell


    def _get_short_name(self, f: Union[str, Path, None] = None):
        """
        Convert a name of this format:
        f = Path('/Volumes/Pegasus/ManisLab_Data3/Kasten_Michael/NF107Ai32Het/2017.05.01_000/slice_000/cell_000/CCIV_short_000')
        To:
        2017.05.01~S0C0~CCIV_short_000
        """
        if f is None:
            raise ValueError(f"Input file to get_shortname is NONE")
        f = Path(f)
        protocol = f.stem
        cellp = f.parent
        cell = cellp.stem
        slicenp = cellp.parent
        slicen = slicenp.stem
        day = slicenp.parent
        foname = str(
            Path(day.name[:10], "S" + slicen[-1] + "C" + cell[-1], protocol)
        ).replace("/", "~")
        return foname

    def set_data_name(self, datatype:str):
        self.AR.setDataName("MultiClamp1.ma")

    def ISO8601_age(self, agestr):
        """Convert somewhat random age designators to ISO standard, e.g.:
            postnatal day 30 mouse = P30D  (or P30W, or P3Y)
            Ranges are P1D/P3D if bounded, or P12D/ if not known but have lower bound.

        Params:
            agestr (str): age string from the file

        Returns:
            str: sanitized age string
        """

        agestr = agestr.replace('p', 'P')
        agestr = agestr.replace('d', 'D')
        if 'P' not in agestr:
            agestr = 'P' + agestr
        if 'D' not in agestr:
            agestr = agestr + "D"
        if agestr == "PD":
            agestr = "P9999D"  # no age specified
        return agestr


    def acq4tonwb(self, protocolname, outfilename: Union[Path, str, None] = None,
            recordingmode:Literal["CC", "VC", "I=0"]="CC"):
        print("datafile: ", protocolname, protocolname.is_dir())
        if outfilename is None:
            outfilename = self._get_short_name(protocolname)
        print("Out file: ", outfilename)

        self.AR.setProtocol(protocolname)
        dataok = self.AR.getData()

        if not dataok:
            return
        info = self.AR.readDirIndex(currdir=protocolname.parent.parent.parent)["."]
        slice_index = self.AR.readDirIndex(currdir=protocolname.parent.parent)["."]
        cell_index = self.AR.readDirIndex(currdir=protocolname.parent)["."]
        proto_index = self.AR.readDirIndex(currdir=protocolname)["."]
        #
        # data in self.AR.data_array, self.AR.time_base, stimuli in self.AR.cmd_wave
        data_date = datetime.date.fromtimestamp(info["__timestamp__"])
        try:
            age = int(
            info["age"]
            .replace("p", "")
            .replace("P", "")
            .replace("~", "")
            .replace("+", "")
            )
        except:
            age = None
        if "sex" not in info.keys() or info["sex"] == "":
            info["sex"] = "U"
        else:
            info["sex"] = info["sex"].upper()
        if "weight" not in info.keys():
            info["weight"] = None
        if "species" not in info.keys() or info["species"] == "":
            info["species"] = 'Mus musculus'
        if age is not None:
                dob = datetime.datetime.combine(
                data_date - datetime.timedelta(days=age), datetime.time()
            )
            # dobstr = dob.strftime("%d/%m/%Y%Z")
        else:
            dob = None
        if "mouse" in info["species"] or "Mouse" in info["species"]:
            info["species"] = 'Mus musculus'
        subject_id = info['animal identifier']
        if subject_id is None or subject_id == "":
            dset = Path(protocolname).parts
            subject_id = str(Path(*dset[-6:-3]))

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

        #     print(info)
        if "type 1" not in list(cell_index.keys()):
            ctypes = "Not specified"
        else:
            ctypes = f"{cell_index['type 1']:s} and {cell_index['type 2']:s}"
        
        if "notes" not in info.keys():
            info['notes'] = "  No Notes"

        self.NWBFile = NWB.NWBFile(
            session_start_time=datetime.datetime.fromtimestamp(info['__timestamp__'], tz=tzlocal()),
            session_id=f"{0:d}",
            session_description=info['description'],
            keywords=["mouse", "intracellular recording", "channelrhodopsins", "pathway tracing",
                "auditory cortex", "auditory thalamus", "medial geniculate", "inferior colliculus",
                "brachium of the inferior colliculus", "cochlear nucleus", "dorsal cochlear nucleus",
                "AAV", "synaptic potentials", "optogenetics"],
            notes = f"Cell Type: {ctypes:s}\n  Full Protocol Path: {str(protocolname):s}\n" + info['notes'],
            identifier=info['animal identifier'],
            protocol = str(protocolname.name),
            timestamps_reference_time=datetime.datetime.now(tzlocal()),
            experimenter="Kasten, Michael R.",
            lab="Manis Lab",
            institution="UNC Chapel Hill",
            experiment_description="RF1NS128873",
            subject=subject,
        )
        device = NWB.device.Device("MC700B", "Current and Voltage clamp amplifier", "Axon Instruments (Molecular Devices)")
        self.NWBFile.add_device(device)

        if recordingmode == "CC":
            self.get_CC_data(protocolname=protocolname, info=info, device=device)
        elif recordingmode == "VC":
            self.get_VC_data(protocolname=protocolname, info=info, device=device)
        else:
            print(f"Recording mode {recordingmode:s} is not implemented yet")
            return None
        outfile = Path("test_data", str(outfilename) + ".nwb")
        print(f"writing to: {str(outfile):s}")
        with NWB.NWBHDF5IO(str(outfile), "w") as io:
            io.write(self.NWBFile)
        return outfile

    def get_one_CC_recording(self, ccseries:str, 
            protocolname:Path,
            info:dict, device:NWB.device,
            acq4_source_name:str="MultiClamp1.ma",
            electrode:int=1):
        
        self.AR.setDataName(acq4_source_name)
        print(f"\nLooking for a recordings from {acq4_source_name:s}")
        dataok = self.AR.getData(check=True)
        if not dataok:
            print("   source not found, skipping")
            return
        self.AR.getData()
        slicen, cell = self._get_slice_cell(f=protocolname)
        s = []
        for d in self.AR.clampInfo['dirs']:
            datainfo = self.AR.getDataInfo(Path(d, self.AR.dataname))
            # print('datainfo: ', datainfo, "\n")
            # s.append(AR.getStim(Path(d, self.AR.dataname)))

        elec1 = NWB.icephys.IntracellularElectrode(
            name = f"elec{electrode:d}",
            description="Sutter 1.5 mm patch as intracellular electrode",
            cell_id = f"{slicen:s}_{cell:s}",
            filtering = str(datainfo[1]['ClampState']['ClampParams']['PrimarySignalLPF']),
            device=device,
        )
        self.NWBFile.add_icephys_electrode(elec1)  # not well documented!
        
        step_times = [np.uint64(self.AR.tstart*self.AR.samp_rate), np.uint64(self.AR.tend*self.AR.samp_rate)]
        istim = NWB.icephys.CurrentClampStimulusSeries(
            name = f"Ics{electrode:d}",
            stimulus_description = str(self.AR.protocol.name),
            description = "control values are the current injection step times in seconds",
            control = step_times,
            control_description = ["tstart", "tend"],
            comments = "CC",
            data = np.array(self.AR.cmd_wave).T,
            unit = "amperes",
            starting_time = info["__timestamp__"],
            rate = self.AR.sample_rate[0],
            electrode = elec1,
            gain = datainfo[1]['ClampState']['extCmdScale'],
            sweep_number = np.uint32(1),
        )

        vdata = NWB.icephys.CurrentClampSeries(
            name = "Vcs1",
            description = ccseries,
            data = self.AR.data_array.T,
            unit = "volts",
            electrode = elec1,
            gain = datainfo[1]['ClampState']['primaryGain'],
            bias_current = datainfo[1]['ClampState']['ClampParams']['Holding'],
            bridge_balance = datainfo[1]['ClampState']['ClampParams']['BridgeBalResist'],
            capacitance_compensation = datainfo[1]['ClampState']['ClampParams']['NeutralizationCap'],
            stimulus_description = "Current Steps",
            resolution = np.NaN,
            conversion = 1.0,
            timestamps = None,
            starting_time = info["__timestamp__"],
            rate = self.AR.sample_rate[0],
            comments = "CC",
            control = None,
            control_description = None,
            sweep_number = np.uint(1),
        )
        self.NWBFile.add_acquisition(istim)
        self.NWBFile.add_acquisition(vdata)

    def get_CC_data(self, protocolname:Path, info:dict, device:NWB.device):
        """Check if a protocol is a one of our known
        current-clamp protocols, and if so, put the data into the
        NWB format

        This adds the cc protocol information directly to the self.NWBFile object.

        Returns:
            Nothing
        """
        protoname = protocolname.name
        ccseries_description = None
        if protoname.startswith("CCIV"):
            ccseries_description = "Current Clamp IV series"
        elif protoname.startswith("Ic_LED"):
            ccseries_description = "Current Clamp with LED Illumination for optical stimulation"
        elif protoname.startswith("Map_NewBlueLaser_IC"):
            ccseries_description = "Current Clamp with Laser scanning photostimulation"
        if ccseries_description is None:
            return None
        for ampchannel in [1,2]:
            self.get_one_CC_recording(protocolname=protocolname,
                ccseries=ccseries_description,
                info=info, 
                device = device, 
                acq4_source_name=f"MultiClamp{ampchannel:d}.ma",
                electrode = ampchannel)

 
    def get_one_VC_recording(self, vcseries:str, 
            protocolname:Path,
            info:dict, device:NWB.device,
            acq4_source_name:str="MultiClamp1.ma",
            electrode:int=1):


        slicen, cell = self._get_slice_cell(f=protocolname)
        s = []
        for d in self.AR.clampInfo['dirs']:
            datainfo = self.AR.getDataInfo(Path(d, self.AR.dataname))
            # print('datainfo: ', datainfo, "\n")
            # s.append(AR.getStim(Path(d, self.AR.dataname)))

        elec = NWB.icephys.IntracellularElectrode(
            name=f"elec{electrode:d}",
            description="Sutter 1.5 mm patch as intracellular electrode",
            cell_id = f"{slicen:s}_{cell:s}",
            filtering = str(datainfo[1]['ClampState']['ClampParams']['PrimarySignalLPF']),
            device=device,
        )
        self.NWBFile.add_icephys_electrode(elec)  # not well documented!
        
        step_times = [np.int64(self.AR.tstart*self.AR.samp_rate), np.int64(self.AR.tend*self.AR.samp_rate)]
        vstim = NWB.icephys.VoltageClampStimulusSeries(
            name = f"Vcs{electrode:d}",
            stimulus_description = str(self.AR.protocol.name),
            description = vcseries,
            control = step_times,
            control_description = ["tstart", "tend"],
            comments = "VC",
            data = np.array(self.AR.cmd_wave).T,
            unit = "volts",
            starting_time = info["__timestamp__"],
            rate=self.AR.sample_rate[0],
            electrode = elec,
            gain = datainfo[1]['ClampState']['extCmdScale'],
            sweep_number = np.uint32(1),
        )
        cparams = datainfo[1]['ClampState']['ClampParams']

        idata = NWB.icephys.VoltageClampSeries(
            name=f"Ics{electrode:d}",
            description = vcseries,
            data = self.AR.data_array.T,
            unit = "amperes",
            electrode = elec,
            gain = datainfo[1]['ClampState']['primaryGain'],
            stimulus_description=str(protocolname.name),
            capacitance_fast = cparams['FastCompCap'],
            capacitance_slow = cparams['SlowCompCap'],
            resistance_comp_bandwidth = cparams['RsCompBandwidth'],
            resistance_comp_correction = cparams['RsCompCorrection'],
            resistance_comp_prediction = 0., # not recorded in acq4, we rarely use this.
            whole_cell_capacitance_comp = cparams['WholeCellCompCap'],
            whole_cell_series_resistance_comp = cparams['WholeCellCompResist'],
            resolution=np.NaN,
            conversion=1.0,
            timestamps=None,
            starting_time=info["__timestamp__"],
            rate=self.AR.sample_rate[0],
            comments="VC",
            control=None,
            control_description=None,
            sweep_number=np.uint(1),
            offset = datainfo[1]['ClampState']['holding']
        )
        self.NWBFile.add_acquisition(vstim)
        self.NWBFile.add_acquisition(idata)
    
    def get_VC_data(self, protocolname:Path, info:dict, device:NWB.device):
        """Check if a protocol is a one of our known
        voltage-clamp protocols, and if so, put the data into the
        NWB format

        This adds the vc protocol information directly to the self.NWBFile object.

        Returns:
            Nothing
        """
        protoname = protocolname.name
        print('protoname')
        vcseries_description = None
        if protoname.startswith("VCIV"):
            vcseries_description = "Voltage Clamp series"
        elif protoname.startswith("Vc_LED"):
            vcseries_description = "Voltage Clamp with LED Illumination for optical stimulation"
        elif protoname.startswith("Map_NewBlueLaser_VC"):
            vcseries_description = "Voltage Clamp with Laser scanning photostimulation"
        if vcseries_description is None:
            return None
        print("vcseries_description: ", vcseries_description)

        for ampchannel in [1,2]:
            self.get_one_VC_recording(vcseries=vcseries_description, 
                protocolname=protocolname,
                info=info, device=device,
                acq4_source_name=f"MultiClamp{ampchannel:d}.ma",
                electrode=ampchannel)


def ConvertFile(filename:Union[str, Path], outputpath:str="test_data", device='MultiClamp1.ma', mode="CC"):
    A2N = ACQ4toNWB(outputpath)
    A2N.set_data_name(device)
    NWBFile = A2N.acq4tonwb(f, recordingmode=mode)
    print("Nwb file generated: ", NWBFile)
    results = list(inspect_all(path=NWBFile))
    if len(results) == 0:
        print("Conversion OK")
    else:
        print("Error in conversion detected: ")
        print(results)

if __name__ == "__main__":
#    f = Path('/Volumes/Pegasus/ManisLab_Data3/Kasten_Michael/NF107Ai32_Het/2017.05.01_000/slice_000/cell_000/CCIV_short_000')
    # f = Path('/Volumes/Pegasus_002/ManisLab_Data3/Kasten_Michael/HK_collab_ICinj/Thalamocortical/2022.10.10_000/slice_001/cell_000/CCIV_long_000')
    # f = Path('/Volumes/Pegasus_002/ManisLab_Data3/Kasten_Michael/2017.02.23_000/slice_000/cell_000/Map_NewBlueLaser_VC_single_MAX_005')
    f = Path('/Volumes/Pegasus_002/ManisLab_Data3/Kasten_Michael/NF107Ai32_Het/2019.04.16_000/slice_002/cell_000/VCIV_simple_001')
    ConvertFile(filename=f, mode="VC")
