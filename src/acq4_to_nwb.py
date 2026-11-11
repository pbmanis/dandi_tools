"""
This program provides functions to tead and convert acq4 data to NWB format that
is acceptable for the Dandi repository. It was specifically designed for the
auditory cortex experiments with Kato lab, but should be generally useful for
converting other kinds of acq4 experiments. A single file can be converted by
calling ConvertFile(filename) The conversion is checked after generating the
outputfile with nwbinspector.

V0.2.0, 24 June 2024

NWBFile expects: class pynwb.file.NWBFile(
    Required:
        session_description (str) :: a description of the session where this
        data was generated identifier (str) :: a unique text identifier for the
        file session_start_time (datetime) :: the start date and time of the
        recording session

    Optional (we try to include these):
        experimenter (tuple or list or str) :: name of person who performed
        experiment experiment_description (str) :: general description of the
        experiment institution (str) :: institution(s) where experiment is
        performed lab (str) :: lab where experiment was performed source_script
        (str) :: Script file used to create this NWB file.
        source_script_file_name (str) :: Name of the source_script file

        devices (list or tuple) :: Device objects belonging to this NWBFile
        subject (Subject) :: subject metadata protocol (str) :: Experimental
        protocol, if applicable. E.g., include IACUC protocol acquisition (list
        or tuple) :: Raw TimeSeries objects belonging to this NWBFile stimulus
        (list or tuple) :: Stimulus TimeSeries objects belonging to this NWBFile
        : Electrical opto: OptogeneticStimulusSites and waveforms that belong to
        this NWBFile. lab_meta_data (list or tuple) :: an extension that
        contains lab-specific meta-data intracellular_recordings
        (IntracellularRecordingsTable) ::
            the IntracellularRecordingsTable table that belongs to this NWBFile
        icephys_simultaneous_recordings (SimultaneousRecordingsTable) ::
            the SimultaneousRecordingsTable table that belongs to this NWBFile
        icephys_sequential_recordings (SequentialRecordingsTable) ::
            the SequentialRecordingsTable table that belongs to this NWBFile
        icephys_repetitions (RepetitionsTable) ::
            the RepetitionsTable table that belongs to this NWBFile
        icephys_experimental_conditions (ExperimentalConditionsTable) ::
            the ExperimentalConditionsTable table that belongs to this NWBFile

        units (Units) :: A table containing unit metadata

        pharmacology (str) :: Description of drugs used, including how and when
             they were administered. Anesthesia(s), painkiller(s), etc., plus
             dosage, concentration, etc.
        virus (str) :: Information about virus(es) used in experiments,
            including virus ID, source, date made, injection location, volume,
            etc.

    Optional and not included:
        file_create_date (ndarray or list or tuple or Dataset or StrDataset
             or HDMFDataset or AbstractDataChunkIterator or datetime) :: the
             date and time the file was created and subsequent modifications
             made
        timestamps_reference_time (datetime) ::
            date and time corresponding to time zero of all timestamps; defaults
            to value of session_start_time
        session_id (str) :: lab-specific ID for the session keywords (ndarray or
        list or tuple or Dataset or StrDataset or HDMFDataset
            or AbstractDataChunkIterator) :: Terms to search over
        related_publications (tuple or list or str) :: Publication
        information.PMID, DOI, URL, etc. If multiple, concatenate together and
        describe which is which. such as PMID, DOI, URL, etc slices (str) ::
        Description of slices, including information about preparation
        thickness, orientation, temperature and bath solution data_collection
        (str) :: Notes about data collection and analysis. surgery (str) ::
        Narrative description about surgery/surgeries, including date(s) and who
        performed surgery. stimulus_notes (str) :: Notes about stimuli, such as
        how and where presented. analysis (list or tuple) :: result of analysis
        stimulus (list or tuple) :: Stimulus TimeSeries objects belonging to
        this NWBFile stimulus_template (list or tuple) :: Stimulus template
        TimeSeries objects belonging to this NWBFile epochs (TimeIntervals) ::
        Epoch objects belonging to this NWBFile epoch_tags (tuple or list or
        set) :: A sorted list of tags used across all epochs trials
        (TimeIntervals) :: A table containing trial data invalid_times
        (TimeIntervals) :: A table containing times to be omitted from analysis
        intervals (list or tuple) :: any TimeIntervals tables storing time
        intervals processing (list or tuple) :: ProcessingModule objects
        belonging to this NWBFile electrodes (DynamicTable) :: the
        ElectrodeTable that belongs to this NWBFile electrode_groups (Iterable)
        :: the ElectrodeGroups that belong to this NWBFile sweep_table
        (SweepTable) :: the SweepTable that belong to this NWBFile
        imaging_planes (list or tuple) :: ImagingPlanes that belong to this
        NWBFile ogen_sites (list or tuple) :: OptogeneticStimulusSites that
        belong to this NWBFile

        scratch (list or tuple) :: scratch data
       icephys_filtering (str) :: [DEPRECATED] Use
       IntracellularElectrode.filtering instead. Description of filtering used.
        ic_electrodes (list or tuple) :: DEPRECATED use icephys_electrodes
        parameter instead. IntracellularElectrodes that belong to this NWBFile


    Mapping ACQ4 to NWB:
    ====================
    Acq4 stores data in a hierarchical directory structure.
        The top level is the DAY of the experiment. The second level is the
        SLICE number. The third level is the CELL number. The fourth level is
        the PROTOCOL name, with an appended number in the format "_000".
            Repeats of a protocol will increment this number, although if a
            protocol is ended early, the sequence may not have fixed steps.
        Within a protocol, there is a series of subdirectories, each containing
            the data for a single sweep.
        The protocol "sweep" subdirectories are numbered according to the
        nesting of repititions. If there is only one repitition, the
        subdirectory is named "000". If there are multiple repitions, the
        Subdirectories are named "000_000", "001_000", (for example showing two
        reps of the same stimulus).
            etc. If there is another parameter that is varied, it is added to
            the name, e.g. "000_000_000" (this is rarely if ever done).

    Mapping this to the NWB structure:
        1. Experimental conditions table - not used here.
        2. Repetitions Table: one row for each protocol type, holding the
           repeats of that protocol (e.g., CCIV_long_000 and CCIV_long_001).
        3. Sequential Recordings Table: one row for each acq4 "protocol",
           holding its sweeps (usually parametric variation of current, laser
           spot position, or other stimulus parameter).
        4. Simultaneous Recordings Table: one row per sweep. We record from a
           single electrode, so each row holds one intracellular recording.
        5. Intracellular Recordings Table: one row per sweep, pairing the
           stimulus and response series with the electrode (one electrode per
           amplifier channel for the cell).
        6. PatchClampSeries: This is the lowest level, that corresponds to a
           single trial or sweep in acq4 - typically the current or voltage
           trace, and the stimulus command. This is the level at which we
           store the traces.

        Also included are optical stimulation data, stored in ogen structures
        (OptogeneticStimulusSite and OptogeneticSeries).

We use the NAME of each entry to parse and reassemble the data for display. See
display_nwb.py for an example of how we use this. The name is formatted as a
string that can be parsed (in regex) to extract the protocol, data type, and
sweep number. The protocol and sweep numbers are used to bind the stimulus
(optical, electrical) and recordings (electrical) together. The optogenetic data
in the ogen structures similarly has names that match and hold other parameters
such as the spot location (if using the laser with a scanning mirror), power,
etc. )


UPLOAD TO DANDI
----------------
1. make sure you have the latest DANDI (pip install dandi==0.67.3 or whatever is
   latest)
2. maker sure you have the latest dandischema (pip install dandischema==0.11.0
   or whatever is latest)
 Then at the terminal:
```
1. dandi download https://dandiarchive.org/dandiset/<dataset_id>/draft
2. cd <dataset_id> dandi organize <source_folder> -f dry   # dry run to see if it works
3. dandi organize <source_folder>
4. dandi validate dandi upload

``` Support::
    (this is also information that needs to be manually filled out at DANDI)

    NIH grants:
    DC RF1 NS128873 (Kato, Manis/Kasten, MPI, 2022-) Cortical circuits for the integration of parallel short-latency auditory pathways.
    DC R01 DC019053 (Manis, 2020-2025) Cellular mechanisms of auditory information processing.

    Ethics:
    iacuc@med.unc.edu
    IACUC protocol 24.083 (Manis)
    https://research.unc.edu/iacuc/

Copyright 2022-2025 Paul B. Manis Distributed under MIT/X11 license. See
license.txt for more information.
ORCID:0000-0003-0131-8961
ROR.org id: https://ror.org/0130frc33

"""

import argparse
import collections
import concurrent.futures
import datetime
import uuid
from dataclasses import dataclass, field
from pathlib import Path
import re
from typing import List, Literal, Tuple, Union

import numpy as np
import pynwb as NWB
from dateutil.tz import tzlocal

# from ephys import datareaders as DR
from nwbinspector import Importance, inspect_all, inspect_nwbfile, load_config
from pynwb import NWBHDF5IO
from pynwb.core import DynamicTable, VectorData
import pyqtgraph.configfile
from src.acq4_reader import acq4_reader
from pylibrary.tools import cprint
CP = cprint.cprint

# the project directory (holds src/ and config/); relative paths in the configuration
# are relative to this directory.
project_root = Path(__file__).resolve().parent.parent
default_config_file = Path(project_root, "config", "dandi_configuration.cfg")


def read_config(config_file: Union[str, Path, None] = None) -> dict:
    """Read the dataset configuration file (pyqtgraph configfile format)

    Parameters
    ----------
    config_file : str or pathlib.Path, optional
        the configuration file. Defaults to config/dandi_configuration.cfg
        in the project directory.

    Returns
    -------
    dict
        Dataset name -> its settings (see config/dandi_example_configuration.cfg).

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    """
    config_file = Path(config_file) if config_file is not None else default_config_file
    if not config_file.is_file():
        raise FileNotFoundError(
            f"Configuration file {config_file!s} not found; see config/dandi_example_configuration.cfg"
        )
    return pyqtgraph.configfile.readConfigFile(str(config_file))

# re_datasets = re.compile(r"^([a-z_]){1,3}", re.IGNORECASE) # (?P<id>)[0-9]{6}$", re.IGNORECASE)  # match the acq4 data files

# m = re_datasets.match("HK_Cre_ugh_001422")  # should be True
# print(m)
# print(m.group(2))
# print(m.groups())
# print(m.group("id"))  # should be 001422

re_pnd = re.compile(
    r"^P?(?P<age>[0-9]{1,9})\s*(?P<greaterthan>[\+\>]?)\s*(?P<ish>ish|\?)?\s*(?P<unit>[DYMW]?)$",
    re.IGNORECASE,
)  # match the acq4 age strings
age_unit_days = {"D": 1.0, "W": 7.0, "M": 30.44, "Y": 365.25}


def parse_age(agestr):
    """Parse free-form age designators, e.g.:
        postnatal day 30 mouse = P30D  (or P30W, or P3Y)
        An age with a lower bound only (P150+) becomes the open interval P150D/.
        Approximate ages (P150ish, P23?) are not valid ISO 8601, so the qualifier
        is returned separately and the ISO age is given without it.

    Parameters
    ----------
    agestr : str
        age string from the file

    Returns
    -------
    iso_age : str or None
        ISO 8601 age string, or None if no age was given
    age_days : float or None
        age in days, or None if the age is not exact
    qualifiers : list of str
        list of qualifiers: "lower bound", "approximate"

    Raises
    ------
    ValueError
        If the age string is not in a recognized format.
    """
    if agestr is None or str(agestr).strip() == "":
        return None, None, []
    age_match = re_pnd.match(str(agestr).strip())
    if age_match is None:
        raise ValueError(f"Age string {agestr!s} is not in a recognized format; please fix in the acq4 data.")
    unit = (age_match["unit"] or "D").upper()
    iso_age = f"P{age_match['age']:s}{unit:s}"
    age_days = int(age_match["age"]) * age_unit_days[unit]
    qualifiers = []
    if age_match["greaterthan"]:
        iso_age += "/"  # ISO 8601 interval with an open upper bound
        qualifiers.append("lower bound")
        age_days = None
    if age_match["ish"]:
        qualifiers.append("approximate")
        age_days = None
    return iso_age, age_days, qualifiers


def ISO8601_age(agestr):
    """Convert free-form age designators to an ISO 8601 duration (see parse_age)

    Parameters
    ----------
    agestr : str
        age string from the file

    Returns
    -------
    str or None
        The ISO 8601 age.
    """
    return parse_age(agestr)[0]


def sex_code(sex):
    """Map free-form sex designators to the NWB/DANDI codes M, F, O or U

    Parameters
    ----------
    sex : str or None
        The sex as written in the data (only the first letter is used).

    Returns
    -------
    str
        "M", "F", "O" or "U" (unknown, also for anything unrecognized).
    """
    sex = str(sex).strip().upper() if sex is not None else ""
    if sex[:1] in ["M", "F", "O", "U"]:
        return sex[:1]
    return "U"


re_weight = re.compile(r"^\s*(?P<value>\d+\.?\d*)\s*(?P<unit>kg|grams|g)?", re.IGNORECASE)


def weight_with_units(weight):
    """Convert free-form weight designators to a standard format, e.g.:
        20 g, 20g, 20 grams, 0.02 kg, 20 -> "20.0 g"

    Returns None if the weight is missing, zero (not measured), or cannot be interpreted.

    Parameters
    ----------
    weight : str, int or float
        The weight as written in the data.

    Returns
    -------
    str or None
        The weight in grams, e.g. "20.0 g".
    """
    if isinstance(weight, (int, float)):
        value = float(weight)
    elif isinstance(weight, str):
        wgroups = re_weight.match(weight)
        if wgroups is None:
            return None
        value = float(wgroups["value"])
        if (wgroups["unit"] or "g").lower() == "kg":  # reexpress in grams
            value *= 1000.0
    else:
        return None
    if np.isnan(value) or value <= 0.0:
        return None
    return f"{value} g"


# The acq4 protocols that can be converted, by name prefix:
#     (recording mode, series description, stimulus type)
known_protocols = {
    "CCIV": ("CC", "Current Clamp IV series", "current steps"),
    "Ic_LED": ("CC", "Current Clamp with LED Illumination for optical stimulation", "LED"),
    "Map_NewBlueLaser_IC": ("CC", "Current Clamp with Laser scanning photostimulation", "laser scanning"),
    "VCIV": ("VC", "Voltage Clamp series", "voltage steps"),
    "Vc_LED": ("VC", "Voltage Clamp with LED Illumination for optical stimulation", "LED"),
    "VC_LED": ("VC", "Voltage Clamp with LED Illumination for optical stimulation", "LED"),
    "Map_NewBlueLaser_VC": ("VC", "Voltage Clamp with Laser scanning photostimulation", "laser scanning"),
}


def protocol_type(protocol: str):
    """Return (recording mode, series description, stimulus type) for a known protocol, or None

    Parameters
    ----------
    protocol : str
        The acq4 protocol name; it is matched by prefix against ``known_protocols``.

    Returns
    -------
    tuple of str or None
        (recording mode "CC" or "VC", series description, stimulus type).
    """
    for prefix, ptype in known_protocols.items():
        if protocol.startswith(prefix):
            return ptype
    return None

def def_lister():
    """Default factory for list fields of :class:`ExperimentInfo`.

    Returns
    -------
    list
        A new empty list.
    """
    return []


@dataclass
class ExperimentInfo:
    """Data class to hold the metadata for the experiment.

    Attributes
    ----------
    description, protocol, time, lab, institution, experiment_description, sessionid, notes, subject : str
        Session metadata.
    experimenters : list of str
        The experimenters.
    """

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
    """Convert the recordings of acq4 cell directories to NWB files (one file per cell).

    See the module docstring for how the acq4 directory structure maps to the NWB icephys tables.

    Parameters
    ----------
    output_path : str or pathlib.Path, optional
        Directory for the NWB files.
    overwrite : bool, optional
        Replace existing NWB files (default False: skip cells that have one).
    testmode : bool, optional
        Read the data and metadata, but do not write files (default False).
    """
    def __init__(
        self,
        output_path: Union[str, Path, None] = None,
        overwrite: bool = False,
        testmode: bool = False,
    ):
        """Convert a file from ACQ4 (www.acq4.org) format to NWB format (www.nwb.org)
        See the documentation above for the mapping between these formats.

        Parameters
        ----------
        output_path : str or pathlib.Path, optional
            The path to the output files. Defaults to None.
        overwrite : bool, optional
            Replace existing NWB files (default False).
        testmode : bool, optional
            Do not write files (default False).
        """
        self.AR = acq4_reader()  # get acq4 reader from the ephys package.

        self.out_file_path = output_path
        self.overwrite = overwrite
        self.testmode = testmode
        # All of the data here is single electrode, and "MultiClamp1.ma" is the
        # name of the data file holding the ephys data. The ".ma" acq4 files are held in "metaarray" format, but
        # these are basically hdf5 files. They are accompanied by .index files with additional metadata.
        # if a different electrode/amplifier is used, change [1] to [2] etc. in self.ampchannels
        self.set_amplifier_name("MultiClamp1.ma")  # also sets self.ampchannels
        self.laser_device = None
        self.LED_device = None
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

        Parameters
        ----------
        f : str or pathlib.Path
            Full file path.

        Returns
        -------
        tuple of str
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
        """Convert a data path of this format (used by acq4):
        f = Path('/Volumes/Pegasus/all/the/leading/pathnames/2017.05.01_000/slice_000/cell_000')
        To:
        2017.05.01~S0C0 ("short form")

        Parameters
        ----------
        f : str or pathlib.Path
            Full file path.
        expt_id : str or None
            an experiment label to put in the filename.
            ignored if None
        rig_id : str or None
            a rig label to put in the filename
            ignored if None

        Returns
        -------
        str
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

        The amplifier device in acq4 (e.g., "MultiClamp1") is the file name without
        the extension, and the trailing number is the amplifier channel (electrode).

        Parameters
        ----------
        amplifier_name : str
            name of the amplifier data file, e.g. "MultiClamp1.ma"
        """
        self.AR.setDataName(amplifier_name)
        amp_match = re.match(r"^(?P<prefix>.*?)(?P<channel>\d+)$", Path(amplifier_name).stem)
        if amp_match is None:
            raise ValueError(f"Amplifier name {amplifier_name:s} must end in a channel number, e.g. MultiClamp1.ma")
        self.amplifier_prefix = amp_match["prefix"]  # e.g. "MultiClamp"
        self.ampchannels = [int(amp_match["channel"])]




    def find_name_in_path(self, sourcepath: Union[str, Path], name: Union[str, None] = None):
        """find_name_in_path Search the path string for a name, and return the name if found.
        This is used to find "Rig", "Experiment", "Slice", "Cell" etc. in the path string.

        Parameters
        ----------
        sourcepath : str or pathlib.Path
            The full path
        name : str or None, optional
            The string to search for in the path, by default None

        Returns
        -------
        str
            the name found in the path

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
        output_path: Union[Path, str] = None,
        outfilename: Union[Path, str, None] = None,
        appendmode: bool = False,
        dataset: dict = None,
    ):
        """Convert one cell directory to an NWB file.
        The NWB file will contain all of the recordings from ONE cell.
        The recording mode (CC or VC) of each protocol is read from the acq4 data.

        Parameters
        ----------
        experiment_name : str
            The dataset name (a top-level entry in the configuration file).
        path_to_cell : str or pathlib.Path
            Full path and name of the cell directory
        protocols : list of str
            The protocol directories to convert.
        output_path : str or pathlib.Path, optional
            Directory for the NWB file.
        outfilename : str or pathlib.Path, optional
            NWB output filename. Defaults to None (made from the path).
        appendmode : bool, optional
            Not used (default False).
        dataset : dict, optional
            The dataset's settings from the configuration file (keywords, experimenters, grants,
            device, light-source powers, ...).

        Returns
        -------
        pathlib.Path or None
            None indicates failure to read the data.
            an outputfile indicates conversion success.
        """
        keywords = dataset.get("keywords", None)
        experimenters = dataset.get("experimenters", None)
        experiment_description = dataset.get("experiment_description", None)
        iacuc_protocols = dataset.get("iacuc_protocols", None)
        assert keywords is not None, "List of keywords must be specified"
        assert experimenters is not None, "List of experimenters/authors must be specified"
        assert (
            experiment_description is not None
        ), "Experiment Description (grant number) must be specified"
        assert iacuc_protocols is not None, "IACUC protocol(s) must be specified"
        assert output_path is not None, "Output path must be specified"
        # nominal light source powers (W), used to scale the optogenetic command waveforms
        self.laser_power_W = float(dataset.get("laser_power_W", 0.070))
        self.led_power_W = float(dataset.get("led_power_W", np.nan))  # NaN if not set
        path_to_cell = Path(path_to_cell)
        output_path = Path(output_path)
        rig_id = self.find_name_in_path(path_to_cell, "Rig")
        expt_id = self.find_name_in_path(path_to_cell, experiment_name)
        if outfilename is None:
            if not output_path.exists():
                print(f"Output path {output_path} is not an existing directory, creating it.")
                output_path.mkdir(parents=True, exist_ok=True)
            outfilename = Path(
                output_path,
                self._get_short_name(path_to_cell, expt_id, rig_id) + ".nwb",
            )
        outfilename = Path(outfilename)

        print("NWB filename: ", outfilename)
        print("overwrite: ", self.overwrite)
        print("testmode: ", self.testmode)

         # Read the acq4 index metadata for the day info, the current slice, and the current cell.
        # assumes that the currdir is the cell directory.
        info = self.AR.readDirIndex(currdir=path_to_cell.parent.parent)["."]
        self.day_index = info
        slice_index = self.AR.readDirIndex(currdir=path_to_cell.parent)["."]
        cell_index = self.AR.readDirIndex(currdir=path_to_cell)["."]

        # The session is the recording from this cell. The day directory timestamp is not
        # reliable (it can be inherited from an earlier day), and the cell directory
        # timestamp can be later than its first protocol, so use the earliest of the cell
        # and protocol directory timestamps.
        protocol_timestamps = [
            self.AR.readDirIndex(currdir=Path(path_to_cell, p.strip()))["."]["__timestamp__"]
            for p in protocols
        ]
        self.session_start_ts = min([cell_index["__timestamp__"]] + protocol_timestamps)
        data_date = datetime.date.fromtimestamp(self.session_start_ts)
        # get the cell file metadata and clean up the representation.
        age, age_days, age_qualifiers = parse_age(info.get("age", None))
        if age is None:
            CP("r", f"No age found for {path_to_cell!s}; please fix in the acq4 data.")
        info["sex"] = sex_code(info.get("sex", ""))
        weight = weight_with_units(info.get("weight", None))
        if weight is None:
            CP(
                "r",
                f"No usable weight ('{info.get('weight', '')!s}') for {path_to_cell!s}; please fix in the acq4 data.",
            )
        if "species" not in info.keys() or info["species"] == "":
            info["species"] = "Mus musculus"
        dob = None
        # only compute a date of birth when the age is exact
        if age_days is not None and age_days < 1200:
            dob = datetime.datetime.combine(
                data_date - datetime.timedelta(days=round(age_days)),
                datetime.time(),
                tzlocal(),
            )
        # dobstr = dob.strftime("%d/%m/%Y%Z")
        if "mouse" in info["species"] or "Mouse" in info["species"]:
            info["species"] = "Mus musculus"
        subject_id = info.get("animal identifier", "No ID")
        if outfilename.is_file():
            if not self.overwrite:
                print(f"File {outfilename} already exists and overwrite is False. Skipping conversion.")
                return outfilename
            print(f"File {outfilename} already exists and overwrite is True; it will be replaced.")
        if self.testmode:
            print(f"Testmode is True: would write {outfilename}, but skipping conversion.")
            return outfilename



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




        dset = Path(path_to_cell).parts
        if subject_id.strip() in [None, "?", "", "NA"]:
            subject_id = str(Path(dset[-3]))  # day directory: one subject per day
        session_id = self._get_short_name(path_to_cell, rig_id=rig_id, expt_id=expt_id)

        if "type 1" not in list(cell_index.keys()):
            ctypes = "Not specified"
        else:
            ctypes = f"{cell_index['type 1']:s} and {cell_index.get('type 2', 'Not specified'):s}"

        if "notes" not in info.keys():
            info["notes"] = "  No Notes"

        # populate the NWB subject object from the acq4 metadata
  

        subject_description = info.get("strain", "Unknown")
        if len(age_qualifiers) > 0:
            subject_description += (
                f"; age recorded as '{info['age']!s}' ({', '.join(age_qualifiers):s})"
            )
        subject = NWB.file.Subject(
            age=age,
            description=subject_description,
            genotype=info.get("genotype", "Unknown"),
            sex=info.get("sex", "U"),
            species=info.get("species", "Unknown"),
            subject_id=subject_id,
            weight=weight,
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
            session_start_time=datetime.datetime.fromtimestamp(self.session_start_ts, tz=tzlocal()),
            session_id=f"{session_id:s}",
            session_description=info["description"],
            keywords=keywords,
            notes=f"Cell Type: {ctypes:s}\n" + info["notes"],
            protocol=iacuc_protocols,  # str(path_to_cell.name),
            experimenter=experimenters,
            lab="Manis Lab",
            institution="UNC Chapel Hill",
            experiment_description=experiment_description,  # grant.
            subject=subject,
        )

        self.NWBFile.add_device(device)
        self.sweep_counter = 0  # cumulative counter of traces/sweeps that are stored in this file.
        self.electrodes = {}  # one electrode per amplifier channel for this cell
        # sequential recordings (one per protocol), grouped by protocol name without the
        # repeat number, so that repeats of a protocol form one repetition.
        self.protocol_sequences = {}

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
            print(
                f"    Datafile: {path_to_protocol!s}\n{' '*14:s}Exists: {path_to_protocol.is_dir()!s}"
            )

            # get the protocol data set
            self.AR.setProtocol(Path(path_to_cell, protocol))
            proto_index = self.AR.readDirIndex(currdir=path_to_protocol)["."]
            self.get_protocol_data(
                path_to_cell=Path(path_to_cell),
                protocol=protocol,
                proto_index=proto_index,
                info=info,
                device=device,
                protonum=protonum,
            )

        # repeats of the same protocol (e.g., CCIV_long_000, CCIV_long_001) form one repetition
        for sequential_recordings in self.protocol_sequences.values():
            self.NWBFile.add_icephys_repetition(sequential_recordings=sequential_recordings)

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
        stimulus_type: str = "square",
    ):
        """Make NWB intracellular electrophysiology tables from a
        single acq4 current clamp or voltage clamp protocol

        Each sweep becomes a stimulus and a response series and a row of the intracellular
        recordings and simultaneous recordings tables; the sweeps together are one sequential
        recording, grouped with the repeats of the same protocol. Optical stimulation, if any, is
        added for each sweep.

        Parameters
        ----------
        recording_mode : str
            Either CC or VC for current or voltage clamp
        series_description : str
            Name for this  series
        path_to_cell : pathlib.Path
            full path of the data, from which the protocol is pulled
        protocol : str
            The protocol directory name.
        info : dict
            acq4 "info" block provided with this protocol
        device : pynwb.device.Device
            NWB device structure
        acq4_source_name : str, optional
            name of device in acq4. Defaults to "MultiClamp1.ma".
        electrode : int, optional
            # of recording electrode. Defaults to 1.
        protonum : int, optional
            Number of the protocol in the cell (default 0).
        stimulus_type : str, optional
            stimulus type for the sequential recordings table

        Returns
        -------
        None
            The recordings are added to ``self.NWBFile``; None is also returned early if the data
            cannot be read.
        """

        match recording_mode:
            case "CC":
                stim_name = f"{protocol:s}~Ics{electrode:d}"
                rec_name = f"{protocol:s}~Vcs{electrode:d}"
                stim_units = "amperes"
                rec_units = "volts"

            case "VC":
                stim_name = f"{protocol:s}~Vcs{electrode:d}"
                rec_name = f"{protocol:s}~Ics{electrode:d}"
                stim_units = "volts"
                rec_units = "amperes"

            case _:
                print(f"Recording mode {recording_mode:s} not implemented")
                return None

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
        nsweeps = self.AR.cmd_wave.shape[0]
        if len(self.AR.clampInfo["dirs"]) == 0 or nsweeps == 0:
            print(f"No sweeps found in protocol {protocol:s}: ", path_to_cell)
            return None
        # trace_index maps each sweep (row of data_array) to its acq4 sweep directory
        sweep_dir_index = self.AR.trace_index[:nsweeps]
        sweep_start_times = [self.AR.trace_StartTimes[i] for i in sweep_dir_index]
        # clamp settings for each sweep, from that sweep's own data file
        datainfos = [
            self.AR.getDataInfo(Path(self.AR.clampInfo["dirs"][i], self.AR.dataname))
            for i in sweep_dir_index
        ]

        # make the electrode: one per amplifier channel for this cell, shared by all protocols
        if electrode not in self.electrodes:
            self.electrodes[electrode] = NWB.icephys.IntracellularElectrode(
                name=f"electrode{electrode:d}",
                description=f"Sutter 1.5 mm patch as intracellular electrode: {self.day_index.get('internal', 'unknown'):s}",
                cell_id=f"{slicen:s}_{cell:s}",
                filtering=str(datainfos[0][1]["ClampState"]["ClampParams"]["PrimarySignalLPF"]),
                device=device,
            )
        elec1 = self.electrodes[electrode]

        step_description = f"; step from {self.AR.tstart:.4f} s to {self.AR.tend:.4f} s"
        optical = self.read_optical_stimulation(protocol=protocol)  # once per protocol

        recordings = []  # keep a list of the recordings

        print("    # sweeps: ", nsweeps)
        for sweepno in range(nsweeps):
            clamp_state = datainfos[sweepno][1]["ClampState"]
            cparams = clamp_state["ClampParams"]
            # time of the sweep relative to the session start (timestamps_reference_time)
            starting_time = float(sweep_start_times[sweepno] - self.session_start_ts)

            match recording_mode:
                case "CC":
                    stimulus = NWB.icephys.CurrentClampStimulusSeries(
                        name=f"{stim_name:s}_{sweepno:d}",
                        data=np.array(self.AR.cmd_wave[sweepno, :]),
                        electrode=elec1,
                        gain=clamp_state["extCmdScale"],  # from acq4
                        stimulus_description=protocol,
                        description=series_description + step_description,
                        comments=recording_mode,
                        unit=stim_units,
                        starting_time=starting_time,
                        rate=self.AR.sample_rate[0],
                        sweep_number=np.uint32(sweepno),
                    )

                    response = NWB.icephys.CurrentClampSeries(
                        name=f"{rec_name:s}_{sweepno:d}",
                        description=series_description,
                        data=self.AR.data_array[sweepno, :],
                        unit=rec_units,
                        electrode=elec1,
                        gain=clamp_state["primaryGain"],
                        bias_current=cparams["Holding"],
                        bridge_balance=cparams["BridgeBalResist"],
                        capacitance_compensation=cparams["NeutralizationCap"],
                        stimulus_description=protocol,
                        conversion=1.0,
                        timestamps=None,
                        starting_time=starting_time,
                        rate=self.AR.sample_rate[0],
                        comments=recording_mode,
                        sweep_number=np.uint32(sweepno),
                    )

                case "VC":
                    stimulus = NWB.icephys.VoltageClampStimulusSeries(
                        name=f"{stim_name:s}_{sweepno:d}",
                        stimulus_description=protocol,
                        description=series_description + step_description,
                        comments=recording_mode,
                        data=np.array(self.AR.cmd_wave[sweepno, :]),
                        unit=stim_units,
                        starting_time=starting_time,
                        rate=self.AR.sample_rate[0],
                        electrode=elec1,
                        gain=clamp_state["extCmdScale"],
                        sweep_number=np.uint32(sweepno),
                    )

                    response = NWB.icephys.VoltageClampSeries(
                        name=f"{rec_name:s}_{sweepno:d}",
                        # the holding potential is metadata, not an offset to the recorded current
                        description=series_description
                        + f"; holding potential {float(clamp_state['holding']):.4f} V",
                        data=self.AR.data_array[sweepno, :],
                        unit=rec_units,
                        electrode=elec1,
                        gain=clamp_state["primaryGain"],
                        stimulus_description=protocol,
                        capacitance_fast=cparams["FastCompCap"],
                        capacitance_slow=cparams["SlowCompCap"],
                        resistance_comp_bandwidth=cparams["RsCompBandwidth"],
                        resistance_comp_correction=cparams["RsCompCorrection"],
                        resistance_comp_prediction=0.0,  # not recorded in acq4, we rarely use this.
                        whole_cell_capacitance_comp=cparams["WholeCellCompCap"],
                        whole_cell_series_resistance_comp=cparams["WholeCellCompResist"],
                        resolution=np.nan,
                        conversion=1.0,
                        timestamps=None,
                        starting_time=starting_time,
                        rate=self.AR.sample_rate[0],
                        comments=recording_mode,
                        sweep_number=np.uint32(sweepno),
                    )

            IR_index = self.NWBFile.add_intracellular_recording(
                electrode=elec1,
                stimulus=stimulus,
                response=response,
                id=sweepno + self.sweep_counter,
            )
            # one electrode: each sweep is its own "simultaneous" recording
            recordings.append(self.NWBFile.add_icephys_simultaneous_recording(recordings=[IR_index]))
            self.add_optical_stimulation(
                optical=optical,
                protocol=protocol,
                sweepno=sweepno,
                dir_index=sweep_dir_index[sweepno],
                starting_time=starting_time,
            )

        self.sweep_counter += nsweeps

        # the sweeps of this protocol form one sequential recording
        IR_sequence_index = self.NWBFile.add_icephys_sequential_recording(
            simultaneous_recordings=recordings,
            stimulus_type=stimulus_type,
        )
        # group with repeats of this protocol (name without the "_000" repeat number)
        protocol_base = re.sub(r"_\d{3}$", "", protocol)
        self.protocol_sequences.setdefault(protocol_base, []).append(IR_sequence_index)

    def get_protocol_data(
        self,
        path_to_cell: Path,
        protocol: str,
        proto_index: dict,
        info: dict,
        device: NWB.device,
        protonum: int = 0,
    ):
        """Check if a protocol is one of our known protocols (see known_protocols),
        and if so, put the data into the NWB format.
        The recording mode is read from the amplifier settings in the acq4 protocol index,
        and must match the mode expected for the protocol.

        This adds the protocol information directly to the self.NWBFile object.

        Parameters
        ----------
        path_to_cell : pathlib.Path
            full path to the cell
        protocol : str
            protocol name
        proto_index : dict
            acq4 index for the protocol directory
        info : dict
            acq4 metadata dictionary
        device : pynwb.device.Device
            The amplifier device.
        protonum : int, optional
            Number of the protocol in the cell (default 0).

        Returns
        -------
        None
        """
        ptype = protocol_type(protocol)
        if ptype is None:
            print(f"    Protocol {protocol:s} is not a known protocol type; skipping")
            return None
        expected_mode, series_description, stimulus_type = ptype
        devices = proto_index.get("devices", {})
        for ampchannel in self.ampchannels:
            amp_name = f"{self.amplifier_prefix:s}{ampchannel:d}"
            if amp_name not in devices:
                print(f"    Amplifier {amp_name:s} was not used in protocol {protocol:s}; skipping")
                continue
            acq4_mode = devices[amp_name]["mode"]
            match acq4_mode:
                case "IC" | "I=0" | "CC":
                    recording_mode = "CC"
                case "VC":
                    recording_mode = "VC"
                case _:
                    print(f"    Recording mode {acq4_mode!s} is not implemented; skipping {protocol:s}")
                    continue
            if recording_mode != expected_mode:
                print(
                    f"    Protocol {protocol:s} was recorded in {acq4_mode:s} mode, "
                    f"but a {expected_mode:s} protocol was expected; skipping"
                )
                continue
            self.get_one_recording(
                recording_mode=recording_mode,
                path_to_cell=path_to_cell,
                series_description=series_description,
                protocol=protocol,
                info=info,
                device=device,
                acq4_source_name=f"{amp_name:s}.ma",
                electrode=ampchannel,
                protonum=protonum,
                stimulus_type=stimulus_type,
            )

    def read_optical_stimulation(self, protocol: str):
        """Read the laser and LED command waveforms for the current protocol.
        This is done once per protocol; add_optical_stimulation then stores each sweep.

        Optical stimulation is stored in NWB ogen structures:
            OptogeneticStimulusSite: the light source and where it was aimed. For laser
                scanning there is one site per sweep (the spot position); for the
                widefield LED there is one site per protocol.
            OptogeneticSeries: the light source command waveform for each sweep. The data
                are the command values as recorded by acq4 (volts); "conversion" scales
                the largest command in the protocol to the nominal source power in watts.
                The power is not calibrated, and assumes a linear relation between command
                and power.

        Parameters
        ----------
        protocol : str
            The protocol name.

        Returns
        -------
        dict
            one entry per light source used in this protocol ("laser", "LED"),
            holding the per-sweep waveforms and the scaling information.
        """
        optical = {}
        #  laser scanning:
        if self.AR.getLaserBlueCommand():
            if not self.AR.getScannerPositions():
                raise ValueError(
                    f"Scanner positions for {protocol:s} could not be read; "
                    f"please check the acq4 data in {self.AR.protocol!s}"
                )
            if self.laser_device is None:
                self.laser_device = self.NWBFile.create_device(
                    name="Oxxius Laser", description="450 nm 70 mW single-photon Laser"
                )
            # use the Pockels cell command if it was recorded, otherwise the shutter
            waveforms = [
                pcell if pcell is not None else shutter
                for pcell, shutter in zip(self.AR.LaserBlue_pCell_sweeps, self.AR.LaserBlue_Raw_sweeps)
            ]
            optical["laser"] = {
                "waveforms": waveforms,
                "rates": self.AR.LaserBlue_sample_rate,
                "power_W": self.laser_power_W,
                "max_command": np.max([np.max(np.abs(w)) for w in waveforms]),
                "spotsize": self.AR.scanner_spotsize,
                "light_source": "450 nm 70 mW Laser",
            }
        # widefield LED:
        if self.AR.getLEDCommand():
            light_source = "470 nm Thorlabs LED"
            if np.isnan(self.led_power_W):
                CP("y", "'led_power_W' is not set for this dataset; LED power will be recorded as NaN")
            if self.LED_device is None:
                self.LED_device = self.NWBFile.create_device(
                    name=f"{light_source:s}", description="Widefield LED"
                )
            waveforms = [np.array(w) for w in self.AR.LED_Raw]
            led_site = NWB.ogen.OptogeneticStimulusSite(
                name=f"{protocol:s}~optosite",
                device=self.LED_device,
                description="Widefield optical stimulation through the objective",
                excitation_lambda=470.0,
                location="x=0 um, y=0 um (widefield, centered on the field of view)",
            )
            self.NWBFile.add_ogen_site(led_site)
            optical["LED"] = {
                "waveforms": waveforms,
                # the sample rate from the LED file index is not reliable; use the time base
                "rates": [1.0 / np.mean(np.diff(tb)) for tb in self.AR.LED_time_base],
                "power_W": self.led_power_W,
                "max_command": np.max([np.max(np.abs(w)) for w in waveforms]),
                "spotsize": 1e-4,  # 100 microns.
                "light_source": light_source,
                "site": led_site,
            }
        return optical

    def add_optical_stimulation(
        self,
        optical: dict,
        protocol: str,
        sweepno: int,
        dir_index: int,
        starting_time: float,
    ):
        """Add the optical stimulation for one sweep to the NWB file

        Parameters
        ----------
        optical : dict
            light source information from read_optical_stimulation
        protocol : str
            protocol name
        sweepno : int
            sweep number (row in the data array)
        dir_index : int
            index of the acq4 sweep directory for this sweep
        starting_time : float
            start of the sweep relative to the session start
        """
        for source, light in optical.items():
            if source == "laser":
                if dir_index >= self.AR.scanner_positions.shape[0]:
                    raise ValueError(
                        f"Scanner position for sweep {sweepno:d} of {protocol:s} could not be read; "
                        f"please check the acq4 data in {self.AR.protocol!s}"
                    )
                # scanner positions are in meters; store them in microns
                x, y = (int(np.round(p * 1e6)) for p in self.AR.scanner_positions[dir_index])
                site = NWB.ogen.OptogeneticStimulusSite(
                    name=f"{protocol:s}~optosite_{sweepno:d}",
                    device=self.laser_device,
                    description="Laser Scanning Photostimulation site",
                    excitation_lambda=450.0,
                    location=f"x={x:d} um, y={y:d} um (scanner coordinates)",
                )
                self.NWBFile.add_ogen_site(site)
                series_name = f"{protocol:s}~OptoSeries_{sweepno:d}"
            else:
                site = light["site"]
                # keep the laser series name unique if both light sources were used
                series_name = f"{protocol:s}~OptoSeries{'LED' if 'laser' in optical else '':s}_{sweepno:d}"
            if light["max_command"] > 0:
                conversion = light["power_W"] / light["max_command"]
            else:
                conversion = np.nan
            ogenseries = NWB.ogen.OptogeneticSeries(
                name=series_name,
                description=(
                    f"Optical stimulation command waveform for {light['light_source']:s} (volts as recorded); "
                    f"conversion scales the largest command in the protocol to the nominal power of "
                    f"{light['power_W']:.4g} W (not calibrated)"
                ),
                data=light["waveforms"][dir_index],
                site=site,
                conversion=conversion,
                starting_time=starting_time,
                rate=float(light["rates"][dir_index]),
                comments=f"Spotsize: {light['spotsize']:e}",
            )
            self.NWBFile.add_stimulus(ogenseries)


# Read the data back in
def validate(testpath, nwbfile):
    """validate Compare the data in memory with what has been written to the disk file

    Parameters
    ----------
    testpath : str or pathlib.Path
        path to the data file on disk
    nwbfile : pynwb.NWBFile
        nwb structure in memory

    Raises
    ------
    AssertionError
        If the intracellular recordings or the experimental conditions differ.
    """
    with NWBHDF5IO(path=testpath, mode="r") as io:
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
    output_path: Union[str, Path] = None,
    dataset: dict = None,
    overwrite: bool = False,
    testmode: bool = False,
):
    # experiment_name = None
    """Convert one cell directory to an NWB file.

    Parameters
    ----------
    experiment_name : str
        The dataset name.
    filename : str or pathlib.Path
        The cell directory.
    protocols : list of str
        The protocols to convert.
    appendmode : bool, optional
        Not used (default False).
    output_path : str or pathlib.Path
        Directory for the NWB file (required).
    dataset : dict
        The dataset's settings from the configuration file (``device`` names the amplifier).
    overwrite : bool, optional
        Replace an existing NWB file (default False).
    testmode : bool, optional
        Do not write the file (default False).

    Returns
    -------
    pathlib.Path or None
        The NWB file written, or None.
    """
    assert output_path is not None
    A2N = ACQ4toNWB(output_path=output_path, overwrite=overwrite, testmode=testmode)
    A2N.set_amplifier_name(dataset['device'])
    NWBFile = A2N.acq4tonwb(
        experiment_name=experiment_name,
        path_to_cell=filename,
        output_path=output_path,
        protocols=protocols,
        appendmode=appendmode,
        dataset=dataset,
    )
    print("NWB file written:  ", NWBFile)
    return NWBFile


def select(protocols: list):
    """Keep the protocols that can be converted (see known_protocols)

    Parameters
    ----------
    protocols : list of str
        Protocol names.

    Returns
    -------
    list of str
        The protocols of known types; the others are listed on the console.
    """
    ok_prots = [p for p in protocols if protocol_type(p) is not None]
    skipped = [p for p in protocols if protocol_type(p) is None]
    if len(skipped) > 0:
        print(f"    Protocols not converted (not a known protocol type): {skipped!s}")
    return ok_prots


def cvt(
    row,
    experiment_name,
    output_path: Union[str, Path] = None,
    dataset: dict = None,
    exclusions: dict = None,
    rigs: Union[Tuple, List, None] = None,
    row_num: int = 0,
    overwrite: bool = False,
    testmode: bool = False,
):
    """Do the file conversion for each complete data set in a given database row

    The protocols are those in the row's "data_complete" column that are of known types, less
    those excluded for the slice or the cell.

    Parameters
    ----------
    row : pandas.Series
        A row containing information about the day (one cell)
    experiment_name : str
        The dataset name.
    output_path : str or pathlib.Path, optional
        Directory for the NWB files.
    dataset : dict, optional
        The dataset's settings from the configuration file.
    exclusions : dict, optional
        "day/slice" or "day/slice/cell" -> {"protocols": [...]} ("all" excludes everything).
    rigs : tuple or list of str, optional
        Convert only data recorded on these rigs (a rig is a directory in the data path).
    row_num : int, optional
        Row number, for messages.
    overwrite : bool, optional
        Replace existing NWB files (default False).
    testmode : bool, optional
        Do not write files (default False).

    Returns
    -------
    nwbfile : pathlib.Path or None
        The NWB file written.
    filename : pathlib.Path or None
        The cell directory.
    """

    path = row.data_directory
    day = row.date
    if rigs is not None:  # restrict by rig: the rig is a directory in the data path
        pathparts = Path(path, day).parts
        if not any(rig in pathparts for rig in rigs):
            print(f"Skipping {day!s}: not recorded on one of the rigs {rigs!s}")
            return None, None  # did not match rigs
    slice = row.slice_slice
    cell = row.cell_cell
    protocols = row.data_complete.split(",")
    protocols = [p.strip() for p in protocols if p.strip() not in [" ", "", None]]
    protocols = select(protocols)
    if len(protocols) == 0:
        print(f"{day}{slice}{cell} No protocols:")
        return None, None
    dsc = f"{day}/{slice}/{cell}"
    ds = f"{day}/{slice}"
    for exclusion_key in [ds, dsc]:  # slice-level and cell-level exclusions
        if exclusion_key not in exclusions.keys():
            continue
        excluded_protocols = exclusions[exclusion_key]["protocols"]
        if any(p in ["all", "All"] for p in excluded_protocols):
            CP("m", f"{dsc} Excluded all protocols ({exclusion_key:s})")
            return None, None
        CP("m", f"{dsc} Excluded protocols: {excluded_protocols!s}")
        protocols = [p for p in protocols if p not in excluded_protocols]
        if len(protocols) == 0:
            CP("m", f"{dsc} No protocols left after exclusion")
            return None, None

    print("\nCell id: ", row.cell_id, " index: ", row_num)
    print("     Protocols: ", protocols)
    filename = Path(path, day, slice, cell)
    NWBFile = ConvertFile(
        experiment_name=experiment_name,
        dataset=dataset,
        filename=filename,
        output_path=output_path,
        protocols=protocols,
        overwrite=overwrite,
        testmode=testmode,
    )
    # NWBFiles.append(NWBFile)
    return NWBFile, filename


def convert_dataset(
    output_path: Union[str, Path],
    experiment_name: str,
    dataset: dict,
    experiment_config: dict,
    exclusions: dict,
    overwrite: bool = False,
    testmode: bool = False,
    rigs: Union[Tuple, List, None] = None,
):
    """Read the pandas database for this set of experiments, and convert each cell.
    The pandas database is generated by dataSummary.py (from github.com/pbmanis/ephys/util),
    and its location is given in the ephys experiment configuration file
    (databasepath/directory/datasummaryFilename).
    This database is a complete record of the experiments, not filtered by any criteria.
    Each row of the database corresponds to one cell.
    The column "data_complete" holds the list of protocols that were completed and are
    potentially suitable for further analysis.

    Parameters
    ----------
    output_path : str or pathlib.Path
        Directory for the NWB files.
    experiment_name : str
        The dataset name.
    dataset : dict
        The dataset's settings from the configuration file.
    experiment_config : dict
        The experiment's entry in the ephys experiments file.
    exclusions : dict
        Protocols to leave out, by slice or cell (see :func:`cvt`).
    overwrite : bool, optional
        Replace existing NWB files (default False).
    testmode : bool, optional
        Do not write files (default False).
    rigs : tuple or list of str, optional
        Convert only data recorded on these rigs.

    Returns
    -------
    list of pathlib.Path
        The NWB files written.
    """
    import pandas as pd

    database = Path(
        experiment_config["databasepath"],
        experiment_config["directory"],
        experiment_config["datasummaryFilename"],
    )
    print(f"Converting acq4 to NWB from database: {database!s}")
    db = pd.read_pickle(database)

    NWBFiles = []
    files = []
    for row_num, (index, row) in enumerate(db.iterrows()):
        print(row.cell_id)
        nwbfile, filename = cvt(
            row,
            experiment_name,
            output_path=output_path,
            dataset=dataset,
            exclusions=exclusions,
            rigs=rigs,
            row_num=row_num,
            overwrite=overwrite,
            testmode=testmode,
        )
        if nwbfile is None:
            continue
        NWBFiles.append(nwbfile)
        files.append(filename)
        print("ACQ4 source file: ", filename)
    return NWBFiles


def check_conversion(nwbf):
    """Inspect one NWB file with nwbinspector, using the DANDI configuration
    (as for "nwbinspector <file> --config dandi").
    The file is reported as OK if there are only best practice suggestions.

    Parameters
    ----------
    nwbf : str or pathlib.Path
        The NWB file.

    Returns
    -------
    bool
        True if the file is OK
    """
    try:
        results = list(inspect_nwbfile(nwbfile_path=nwbf, config=load_config("dandi")))
    except Exception as e:
        print(f"    nwbinspector failed for file {nwbf!s}: {e!r}")
        return False
    problems = [r for r in results if r.importance != Importance.BEST_PRACTICE_SUGGESTION]
    counts = collections.Counter(r.importance.name for r in results)
    if len(problems) == 0:
        print(f"    Conversion OK for file: {nwbf!s} {dict(counts) if len(counts) > 0 else ''}")
        return True
    print(f"    Problems in conversion for file: {nwbf!s} {dict(counts)!s}")
    for r in problems:
        print(f"        {r.importance.name:s} {r.check_function_name:s} {r.object_name!s}: {r.message:s}")
    return False


def check_conversions(output_path: Union[str, Path]):
    """Run nwbinspector on every NWB file in the output directory for a dataset

    Parameters
    ----------
    output_path : str or pathlib.Path
        The directory of NWB files.
    """
    output_path = Path(output_path)
    print("=" * 80)
    print(f"Checking conversions in: {output_path!s}")
    NWBFiles = sorted(output_path.rglob("*.nwb"))
    if len(NWBFiles) == 0:
        print(f"    No NWB files found in {output_path!s}")
        return
    print(f"    {len(NWBFiles):d} NWB files found")
    nok = sum(check_conversion(nwbf) for nwbf in NWBFiles)
    print(f"    {nok:d} of {len(NWBFiles):d} files OK")


def upload_to_dandi(output_path: Path, dandi_path: Path, dataset_id: str):
    """Prepare the converted files for upload to the DANDI archive.
    Downloads the dandiset metadata (dandiset.yaml) into dandi_path/dataset_id,
    organizes the converted NWB files into it (dry run first), and validates.
    Each step must succeed before the next is run.
    The final "dandi upload" is left to be run by hand (see README.md).

    Parameters
    ----------
    output_path : pathlib.Path
        directory holding the converted NWB files
    dandi_path : pathlib.Path
        directory in which the dandiset directory is created
    dataset_id : str
        the DANDI dataset id (e.g. "001423")

    Returns
    -------
    bool
        True if all steps succeeded
    """
    import subprocess

    dandi_path.mkdir(parents=True, exist_ok=True)
    dandiset_dir = Path(dandi_path, dataset_id)
    steps = [
        (
            [
                "dandi",
                "download",
                "--download",
                "dandiset.yaml",
                "--existing",
                "refresh",
                f"https://dandiarchive.org/dandiset/{dataset_id:s}/draft",
            ],
            dandi_path,
        ),
        (["dandi", "organize", str(output_path), "-f", "dry"], dandiset_dir),  # dry run to see if it works
        (["dandi", "organize", str(output_path)], dandiset_dir),
        (["dandi", "validate", "."], dandiset_dir),
    ]
    for cmd, cwd in steps:
        print(f"Running: {' '.join(cmd):s}\n    in {cwd!s}")
        result = subprocess.run(cmd, cwd=cwd)
        if result.returncode != 0:
            print(f"'{' '.join(cmd):s}' failed (exit code {result.returncode:d}); stopping")
            return False
    print(f"Ready to upload. To upload, run:\n    cd {dandiset_dir!s}\n    dandi upload")
    return True


def main():

    """Command line entry point (``acq4_to_nwb <dataset> [-c] [-o] [-v] [-s file] [-u] [-t]``).

    Reads the dataset's entry in the configuration file and the ephys experiments file it names,
    then converts (-c), verifies with nwbinspector (-v, or -s for one file) and prepares the upload
    (-u), as requested. See ``acq4_to_nwb --help``.
    """
    parser = argparse.ArgumentParser(
        description="Convert acq4 data to NWB format",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        type=str,
        dest="dataset",
        help="Dataset name: a top-level entry in the configuration file",
    )
    parser.add_argument(
        "--config",
        type=str,
        default=str(default_config_file),
        dest="config_file",
        help="Dataset configuration file",
    )
    parser.add_argument(
        "-c",
        "--convert",
        default=False,
        action="store_true",
        help="Convert the files to nwb format",
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        default=False,
        action="store_true",
        help="Overwrite existing nwb files (if not set, will skip conversion if nwb file already exists)",
    )
    parser.add_argument(
        "-v",
        "--verify",
        default=False,
        action="store_true",
        help="Verify the converted files",
    )
    parser.add_argument(
        "-s",
        "--spotcheck",
        type=str,
        default=None,
        action="store",
        help="Spot check the converted files",
    )
    parser.add_argument(
        "-u",
        "--upload",
        default=False,
        action="store_true",
        help="Prepare the files for upload to the Dandi archive (download, organize, validate)",
    )
    parser.add_argument(
        "-t",
        "--testmode",
        default=False,
        action="store_true",
        help="Run in test mode (read the metadata but do not write nwb files, and do not run the dandi steps)",
    )
    args = parser.parse_args()
    config = read_config(args.config_file)
    id = args.dataset
    if id not in config.keys():
        print(f"Dataset {id} not found in {args.config_file!s}")
        print(f"Available datasets: {[k for k in config.keys() if k != 'datasets']!s}")
        exit(1)
    else:
        print(f"Dataset {id} found: {config[id]!s}")
    dataset_id = config[id].get("dataset_id", None)
    if dataset_id is None:
        print(f"Dataset {id} has no dataset_id in {args.config_file!s}")
        exit(1)
    dataset_id = str(dataset_id)
    experiments = pyqtgraph.configfile.readConfigFile(config[id]["configfile"])["experiments"]
    # the experiment name in the ephys experiments file, if it differs from the dataset name
    experiment_id = config[id].get("experiment", id)
    if experiment_id not in experiments:
        print(f"Experiment {experiment_id} is not in the experiments file {config[id]['configfile']!s}")
        print(f"Experiments in that file: {list(experiments.keys())!s}")
        print("Set 'experiment' in the dataset configuration to one of these.")
        exit(1)
    experiment_config = experiments[experiment_id]
    exclusions = experiment_config.get("excludeIVs", None) or {}
    # relative output directories are relative to the project directory
    output_path = Path(project_root, config[id]["output_dir"], dataset_id)
    dandi_path = Path(project_root, config[id]["dandi_output_dir"])
    print(f"Output path: {output_path!s}")
    if args.convert:
        convert_dataset(
            output_path=output_path,
            experiment_name=id,
            dataset=config[id],
            experiment_config=experiment_config,
            exclusions=exclusions,
            overwrite=args.overwrite,
            testmode=args.testmode,
        )
    if args.verify:
        check_conversions(output_path=output_path)
    if args.spotcheck is not None:
        check_conversion(Path(output_path, args.spotcheck))  # check the conversion
    if args.upload:
        if args.testmode:
            print("Testmode is True: not running the dandi download/organize/validate steps")
        else:
            upload_to_dandi(output_path=output_path, dandi_path=dandi_path, dataset_id=dataset_id)


if __name__ == "__main__":
    main()
