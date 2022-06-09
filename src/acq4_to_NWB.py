"""
Read acq4 data and convert it from ACQ4 to NWB format
"""
import datetime
import time
from collections import OrderedDict
from pathlib import Path

import ephys.ephysanalysis as EP
import numpy as np
from dateutil.tz import tzlocal

AR = EP.acq4read.Acq4Read()
import pynwb as NWB
from ephys.ephysanalysis import MatdatacRead as MDR


class ACQ4toNWB():
    def __init__(self):
        pass

    def get_short_name(self, f: Union[str, Path, None]: None):
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
        foname = str(Path(day.name[:10], 'S'+slicen[-1]+'C'+cell[-1], protocol)).replace('/', '~')
        return foname

    
    def acq4tonwb(self, protocolname, outfilename:Union[Path, str, None] = None):
        print('datafile: ', protocolname, protocolname.is_dir())
        if outfilename is None:
            outfilename = self.get_short_name(protocolname)
        print('Out file: ', outfilename)

        AR.setProtocol(protocolname)
        AR.setDataName('MultiClamp1.ma')
        dataok = AR.getData()
        sample_rate = AR.sample_rate
        data_array = AR.data_array
        cmd_wave = AR.cmd_wave
        if not dataok:
            return
        info = AR.readDirIndex(currdir=protocolname.parent.parent.parent)['.']
        slice_index = AR.readDirIndex(currdir=protocolname.parent.parent)['.']
        cell_index = AR.readDirIndex(currdir=protocolname.parent)['.']
        proto_index = AR.readDirIndex(currdir=protocolname)['.']


        print(cell_index.keys())
    #     print(proto_index)
        print(info['temperature'], info['sex'], info['age'])

        #
        # data in AR.data_array, AR.time_base, stimuli in AR.cmd_wave
        data_date = datetime.date.fromtimestamp(info['__timestamp__'])
        age = int(info['age'].replace('p', '').replace('P', '').replace('~', '').replace('+', ''))
        dob = datetime.datetime.combine(data_date - datetime.timedelta(days=age), datetime.time())
        dobstr = dob.strftime('%d/%m/%Y')
        subject = NWB.file.Subject(age=info['age'], 
                description=None, genotype='None', sex=info['sex'].upper(), species='Rat', 
                subject_id=None, weight=info['weight'], date_of_birth=dob)

    #     print(info)
        if 'type 1' not in list(cell_index.keys()):
            ctypes = "unknown"
        else:
            ctypes = f"{cell_index['type 1']:s} and {cell_index['type 2']:s}" 
        nwbfile = NWB.NWBFile('Dual Recordings', str(protocolname), datetime.datetime.now(tzlocal()),
                    experimenter='Hughes, Ben; Kasten, Michael',
                    lab='Manis Lab',
                    institution='UNC Chapel Hill',
                    experiment_description='mPFC Pairs',
                    session_id=f"{0:d}",
                    notes=f"Cell Type: {ctypes:s} Protocol: {str(protocolname):s}",
                    subject=subject)
        device = NWB.device.Device('MC700A', parent=None)
        nwbfile.add_device(device)
        elec = NWB.icephys.IntracellularElectrode(name="elec0",
                                    description='KG33 patch as intracellular electrode',
                                    device=device)
        nwbfile.add_ic_electrode(elec)  # not well documented!

        istim1 = NWB.icephys.CurrentClampStimulusSeries( name="Ics1", data=np.array(AR.cmd_wave), unit='amperes',
                starting_time=info['__timestamp__'], rate=sample_rate[0],
                electrode=elec, gain=1.0, sweep_number=1)

        vdata1 = NWB.icephys.CurrentClampSeries('Vcs1', data=AR.data_array, 
            unit='volts', electrode=elec, gain=1.0, bias_current=0.0, bridge_balance=0., 
            capacitance_compensation=0., stimulus_description='NA', 
            resolution=0.0, conversion=1.0, timestamps=None, 
            starting_time=info['__timestamp__'], rate=AR.sample_rate[0], 
            comments='no comments', description='no description', 
            control=None, control_description=None, sweep_number=1, parent=None)

        AR.setDataName('MultiClamp2.ma')
        dataok = AR.getData()
        vdata2 = NWB.icephys.CurrentClampSeries('Vcs2', data=AR.data_array, 
            unit='V', electrode=elec, gain=1.0, bias_current=0.0, bridge_balance=0., 
            capacitance_compensation=0., stimulus_description='NA', 
            resolution=0.0, conversion=1.0, timestamps=None, 
            starting_time=info['__timestamp__'], rate=AR.sample_rate[0], 
            comments='no comments', description='no description', 
            control=None, control_description=None, sweep_number=1, parent=None)
        istim2 = NWB.icephys.CurrentClampStimulusSeries( name="Ics2", data=np.array(AR.cmd_wave), unit='A',
                starting_time=info['__timestamp__'], rate=AR.sample_rate[0],
                electrode=elec, gain=1.0, sweep_number=1)


        nwbfile.add_acquisition(istim1)
        nwbfile.add_acquisition(vdata1)
        nwbfile.add_acquisition(istim2)
        nwbfile.add_acquisition(vdata2)

        outfile = Path(opath, outfilename)
        print(f"writing to: {str(outfile)+'.nwb':s}")
        with NWB.NWBHDF5IO(str(outfile)+'.nwb', 'w') as io:
            io.write(nwbfile)
