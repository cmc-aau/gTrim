.. _files/trimmed.markers.tsv:

trimmed.markers.tsv
===================

A tab-separated file reporting how the columns of each marker were filtered and subsampled, written to
``[prefix].trimmed.markers.tsv``. Markers are listed in the order of the :ref:`marker info file <files/markers_info.tsv>`.
Each column of a marker is counted in exactly one of ``retained_columns``, ``gappy_columns``,
``conserved_columns`` and ``unconserved_columns``, so these sum to ``length``. See :ref:`about/how_trim_works` for how
columns are classified.

The columns in this file are as follows:

* ``marker_id``: identifier of the marker.
* ``description``: description of the marker.
* ``length``: number of columns of the marker in the input alignment.
* ``retained_columns``: number of columns that passed the filters.
* ``retained_columns_perc``: percentage of the marker's columns that passed the filters.
* ``sampled_columns``: number of columns sampled for the trimmed alignment.
* ``gappy_columns``: number of columns removed for having too few residues.
* ``conserved_columns``: number of columns removed for having a consensus at or above ``--max-consensus``.
* ``unconserved_columns``: number of columns removed for having a consensus below ``--min-consensus``.

Produced by
-----------

* :ref:`commands/trim`

Example
-------

.. code-block:: text

    marker_id	description	length	retained_columns	retained_columns_perc	sampled_columns	gappy_columns	conserved_columns	unconserved_columns
    PF00380.20	Ribosomal protein S9/S16	121	84	69.42	50	0	29	8
    PF00410.20	Ribosomal protein S8	126	97	76.98	50	0	17	12
    PF00466.21	Ribosomal protein L10	100	81	81.00	50	4	2	13
