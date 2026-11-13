.. _files/column_stats.tsv:

column_stats.tsv
================

A tab-separated file with the gap and residue statistics of each column of the alignment, written to
``[prefix].column_stats.tsv``. Percentages are of the number of sequences in the alignment and are given to two
decimal places.

The columns in this file are as follows:

* ``column_idx``: index of the column in the alignment, starting at 0.
* ``num_gaps``: number of sequences with a gap (``-`` or ``.``) in the column.
* ``gap_perc``: percentage of sequences with a gap in the column.
* ``num_residues``: number of sequences without a gap in the column.
* ``residue_perc``: percentage of sequences without a gap in the column.
* ``per_residue_counts``: comma-separated ``residue:count`` pairs for each character other than a gap, from most to
  least common. Empty if every sequence has a gap.
* ``per_residue_perc``: comma-separated ``residue:percentage`` pairs, in the same order as ``per_residue_counts``.

Produced by
-----------

* :ref:`commands/column-stats`

Example
-------

.. code-block:: text

    column_idx	num_gaps	gap_perc	num_residues	residue_perc	per_residue_counts	per_residue_perc
    0	10	100.00	0	0.00
    1	9	90.00	1	10.00	C:1	C:10.00
    2	8	80.00	2	20.00	G:2	G:20.00
