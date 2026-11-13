.. _commands/mask:

mask
====

Apply a mask to an alignment, keeping only the columns marked with a ``1``. The mask must have one character for
each column of the alignment (see :ref:`files/mask`). This can be used to apply the mask from :ref:`commands/trim`
or :ref:`commands/prune-gap` to another alignment of the same columns, e.g. one with additional sequences.

Arguments
---------

.. program:: gtrim mask

Required
^^^^^^^^

.. option:: -i, --alignment-file <ALIGNMENT_FILE>

   File containing the multiple sequence alignment.

.. option:: -m, --mask-file <MASK_FILE>

   File containing the binary mask.

.. option:: -o, --out-dir <OUT_DIR>

   Output directory.

Optional
^^^^^^^^

.. option:: -t, --threads <THREADS>

   Number of threads, in the range [1, 1024]. This option currently has no effect and ``mask`` uses a single
   thread.

   Default: ``1``


Files output
------------

* :ref:`[prefix].masked.fasta <files/alignment>`
* :ref:`gtrim.log <files/gtrim.log>`


Example
-------

Input
^^^^^

.. code-block:: bash

    gtrim mask -i trim_test.faa.gz -m out_trim/trim_test.trimmed.mask -o out_mask

Output
^^^^^^

.. code-block:: text

    [2026-11-14 03:16:57] INFO: gtrim v{{ gtrim_version }}
    [2026-11-14 03:16:57] INFO: gtrim mask -i trim_test.faa.gz -m out_trim/trim_test.trimmed.mask -o out_mask
    [2026-11-14 03:16:57] INFO: Using 1 threads.
    [2026-11-14 03:16:57] INFO: Starting mask workflow.
    [2026-11-14 03:16:57] INFO: Reading alignment from file: trim_test.faa.gz
    [2026-11-14 03:16:57] INFO: Alignment contains 200 sequences and 1792 columns.
    [2026-11-14 03:16:57] INFO: Sequence type detected: amino acid
    [2026-11-14 03:16:57] INFO: Reading mask from file: out_trim/trim_test.trimmed.mask
    [2026-11-14 03:16:57] INFO: Mask specifies 500 columns to retain and 1292 columns to filter.
    [2026-11-14 03:16:57] INFO: Writing filtered alignment to file.
    [2026-11-14 03:16:57] INFO: Execution time: 0.01 seconds

This writes ``trim_test.masked.fasta.gz`` and ``gtrim.log`` to ``out_mask``. The masked alignment is identical to
``trim_test.trimmed.fasta.gz`` written by :ref:`commands/trim`.
