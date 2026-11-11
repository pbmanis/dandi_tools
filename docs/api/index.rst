API
===

The modules are in the ``src`` package. Not documented here: ``src.tifffile`` (a copy of
Christoph Gohlke's tifffile, with its own documentation), ``src.imageinfo`` (does not import:
it uses ``nwbinspector.inspect_nwb``, which no longer exists) and ``src.tutorial_test`` (a
pynwb tutorial script that runs when imported).

.. toctree::
   :maxdepth: 1

   acq4_to_nwb
   acq4_reader
   display_nwb
   browse_nwb
   display_acq4
   acq4_mosaic
   stack_projection
   image_kind
   mosaic_exceptions
