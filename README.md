Dandi tools for posting data to the Dandi Archive.

Convert from acq4 to NWB

Test for NWB form correctness

see: https://www.dandiarchive.org/handbook/10_using_dandi/

In particular:
nwbinspector <source_folder> --config dandi
or:
nwbinspector <source_folder> --config dandi --report-file-path <report_location>.txt

then validate:
dandi validate <source_folder>


To upload:
dandi download https://dandiarchive.org/dandiset/<dataset_id>/draft
cd <dataset_id>
dandi organize <source_folder> -f dry
dandi organize <source_folder>
For final:  dandi upload

or for staging:
dandi upload -i dandi-staging
