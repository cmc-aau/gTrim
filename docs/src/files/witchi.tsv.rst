.. _files/witchi.tsv:

WitChi pruned TSV
=================

The ``*_pruned.tsv`` file written by ``witchi prune`` (`WitChi <https://github.com/stephkoest/witchi>`_), which lists
the columns WitChi removed from an alignment, one per line. It is read by :ref:`commands/trim` with ``--mask-file``
(the ``--witchi_info_file`` input of GenomeTreeTk ``msa_filter``).

The file is tab-separated with a header line. Only the ``Original Index`` column is used (the header is matched
case-insensitively): the 0-based index of a removed column in the full concatenated alignment. Each index must be
less than the summed length of the markers in the :ref:`marker info file <files/markers_info.tsv>` and may be listed
only once. A file with no rows removes no columns.

Used by
-------

* :ref:`commands/trim` (``--mask-file``)

Example
-------

.. code-block:: text

    Iteration	Original Index	Global ChiScore	Initial Significant ChiScore	ChiScore Difference	Significant taxa	Null-Pvalue	Algorithm
    1	668	1.0	1.0	0.1	3		quartic
    2	310	1.0	1.0	0.1	3		quartic
