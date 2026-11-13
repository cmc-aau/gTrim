.. _files/pruned.tsv:

pruned.tsv
==========

A tab-separated file with the gap statistics of each column of the input alignment, written to
``[prefix].pruned.tsv``.

The columns in this file are as follows:

* ``column_idx``: index of the column in the input alignment, starting at 0.
* ``num_gaps``: number of sequences with a gap (``-`` or ``.``) in the column.
* ``gap_perc``: percentage of sequences with a gap in the column.
* ``retained``: ``true`` if the column was kept, otherwise ``false``.

Produced by
-----------

* :ref:`commands/prune-gap`

Example
-------

.. code-block:: text

    column_idx	num_gaps	gap_perc	retained
    0	10	100	false
    1	9	90	false
    2	8	80	false
    3	7	70	false
