.. _commands/prune-gap:

prune-gap
=========

Remove columns from an alignment in which the percentage of sequences with a gap is greater than a threshold. Gaps
are the ``-`` and ``.`` characters. A column is kept if its percentage of gaps is less than or equal to
``--gap-threshold``, so the default of 100 keeps every column.

Arguments
---------

.. program:: gtrim prune-gap

Required
^^^^^^^^

.. option:: -i, --alignment-file <ALIGNMENT_FILE>

   File containing the multiple sequence alignment.

.. option:: -o, --out-dir <OUT_DIR>

   Output directory.

Optional
^^^^^^^^

.. option:: --gap-threshold <GAP_THRESHOLD>

   Maximum percentage of gaps to retain a column, in the range [0, 100].

   Default: ``100``

.. option:: --max-residue-pruned <MAX_RESIDUE_PRUNED>

   Maximum number of columns to prune. This option currently has no effect.

   Default: ``100000``

.. option:: -t, --threads <THREADS>

   Number of threads, in the range [1, 1024].

   Default: ``1``


Files output
------------

* :ref:`[prefix].pruned.fasta <files/alignment>`
* :ref:`[prefix].pruned.mask <files/mask>`
* :ref:`[prefix].pruned.tsv <files/pruned.tsv>`
* :ref:`gtrim.log <files/gtrim.log>`


Example
-------

Input
^^^^^

.. code-block:: bash

    gtrim prune-gap -i gap_test.fna.gz -o out_gap --gap-threshold 50

Output
^^^^^^

.. code-block:: text

    [2026-11-14 03:16:57] INFO: gtrim v{{ gtrim_version }}
    [2026-11-14 03:16:57] INFO: gtrim prune-gap -i gap_test.fna.gz -o out_gap --gap-threshold 50
    [2026-11-14 03:16:57] INFO: Using 1 threads.
    [2026-11-14 03:16:57] INFO: Starting gap pruning workflow.
    [2026-11-14 03:16:57] INFO: Reading alignment from file: gap_test.fna.gz
    [2026-11-14 03:16:57] INFO: Alignment contains 10 sequences and 22 columns.
    [2026-11-14 03:16:57] INFO: Sequence type detected: DNA
    [2026-11-14 03:16:57] INFO: Pruning gappy columns.
    [2026-11-14 03:16:57] INFO: Pruned 10 columns.
    [2026-11-14 03:16:57] INFO: Writing filtered alignment to file.
    [2026-11-14 03:16:57] INFO: Execution time: 0.01 seconds

This writes ``gap_test.pruned.fasta.gz``, ``gap_test.pruned.mask``, ``gap_test.pruned.tsv`` and ``gtrim.log`` to
``out_gap``.
