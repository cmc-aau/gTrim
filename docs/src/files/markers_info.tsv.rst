.. _files/markers_info.tsv:

markers_info.tsv
================

A tab-separated file listing the markers of a concatenated marker gene alignment in the order they are concatenated.
It is required by :ref:`commands/trim` (``-m/--marker-info-file``). The GTDB marker info files, e.g.
``gtdb_r226_bac120_markers_info.tsv``, and those of GTDB-Tk can be used directly.

The file must have a header line and the following columns, in any order. Column names are matched ignoring case,
so both ``Marker Id`` and ``Marker ID`` are accepted. Other columns are ignored.

* ``Marker Id``: identifier of the marker.
* ``Description``: description of the marker.
* ``Length (bp)``: number of columns of the marker in the concatenated alignment.

The marker lengths must sum to the length of the alignment.

Used by
-------

* :ref:`commands/trim`

Example
-------

.. code-block:: text

    Marker Id	Name	Description	Length (bp)
    PF00380.20	Ribosomal_S9	Ribosomal protein S9/S16	121
    PF00410.20	Ribosomal_S8	Ribosomal protein S8	126
    PF00466.21	Ribosomal_L10	Ribosomal protein L10	100
