.. _commands/column-stats:

column-stats
============

Report the number and percentage of gaps and of each residue in every column of an alignment. Gaps are the ``-`` and
``.`` characters. No alignment is written.

Arguments
---------

.. program:: gtrim column-stats

Required
^^^^^^^^

.. option:: -i, --alignment-file <ALIGNMENT_FILE>

   File containing the multiple sequence alignment.

.. option:: -o, --out-dir <OUT_DIR>

   Output directory.

Optional
^^^^^^^^

.. option:: -t, --threads <THREADS>

   Number of threads, in the range [1, 1024]. This option currently has no effect and ``column-stats`` uses a
   single thread.

   Default: ``1``


Files output
------------

* :ref:`[prefix].column_stats.tsv <files/column_stats.tsv>`
* :ref:`gtrim.log <files/gtrim.log>`


Example
-------

Input
^^^^^

.. code-block:: bash

    gtrim column-stats -i gap_test.fna.gz -o out_stats

Output
^^^^^^

.. code-block:: text

    [2026-11-14 03:16:57] INFO: gtrim v{{ gtrim_version }}
    [2026-11-14 03:16:57] INFO: gtrim column-stats -i gap_test.fna.gz -o out_stats
    [2026-11-14 03:16:57] INFO: Using 1 threads.
    [2026-11-14 03:16:57] INFO: Starting column stats workflow.
    [2026-11-14 03:16:57] INFO: Reading alignment from file: gap_test.fna.gz
    [2026-11-14 03:16:57] INFO: Alignment contains 10 sequences and 22 columns.
    [2026-11-14 03:16:57] INFO: Sequence type detected: DNA
    [2026-11-14 03:16:57] INFO: Execution time: 0.00 seconds

This writes ``gap_test.column_stats.tsv`` and ``gtrim.log`` to ``out_stats``.
