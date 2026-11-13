.. _commands/trim:

trim
====

Trim a concatenated marker gene alignment, such as the bac120 or ar53 alignments used by the GTDB. Columns are
removed if they are gappy, highly conserved or poorly conserved, and the same number of columns can then be randomly
sampled from each marker to obtain an alignment of a desired length. This reimplements the ``msa_filter`` command of
`GenomeTreeTk <https://github.com/dparks1134/GenomeTreeTk>`_. See :ref:`about/how_trim_works` for a description of
the filters and the subsampling, and how ``trim`` differs from ``msa_filter``.

By default (``--msa-length 1000000``), every column that passes the filters is kept and no columns are subsampled.
Use ``--msa-length 5000`` to subsample a 5,000 column alignment, as done by ``msa_filter`` by default.

If columns were removed from the concatenated alignment before trimming, for example by
`WitChi <https://github.com/stephkoest/witchi>`_, give the columns that were removed with ``--mask-file`` so that
each column of the alignment is assigned to the correct marker.

Arguments
---------

.. program:: gtrim trim

Required
^^^^^^^^

.. option:: -i, --alignment-file <ALIGNMENT_FILE>

   File containing the concatenated marker gene alignment.

.. option:: -m, --marker-info-file <MARKER_INFO_FILE>

   Tab-separated file listing the markers in the order they are concatenated in the alignment
   (see :ref:`files/markers_info.tsv`). The marker lengths must sum to the length of the alignment, or, with
   ``--mask-file``, to the length of the full alignment before columns were removed.

.. option:: -o, --out-dir <OUT_DIR>

   Output directory.

Optional
^^^^^^^^

.. option:: --mask-file <MASK_FILE>

   Columns removed from the full alignment before trimming (e.g. by WitChi); the alignment must contain only the
   retained columns. Either a :ref:`WitChi pruned TSV <files/witchi.tsv>` giving the 0-based index of each removed
   column in an ``Original Index`` column, or a :ref:`binary mask <files/mask>` as used by :ref:`commands/mask` (a
   single line of 0s and 1s, one per column, ``1`` = retained). The format is detected from the content of the
   file. This is the ``--witchi_info_file`` input of ``msa_filter``.

   The mask must have one column for each column described by the marker info file, and the number of columns it
   retains must equal the length of the alignment. The output masks are over the full alignment.

.. option:: --msa-length <MSA_LENGTH>

   Desired length of the alignment after randomly subsampling the same number of columns from each marker. The
   smallest number of columns per marker that gives an alignment of at least this length is used. If fewer columns
   pass the filters, every column passing the filters is kept.

   Default: ``1000000``

.. option:: --min-perc-aa <MIN_PERC_AA>

   Minimum percentage of sequences with a residue (not a gap) required to retain a column.

   Default: ``50``

.. option:: --min-consensus <MIN_CONSENSUS>

   Minimum percentage of residues that must be the most common residue to retain a column (inclusive).

   Default: ``25``

.. option:: --max-consensus <MAX_CONSENSUS>

   Maximum percentage of residues that can be the most common residue to retain a column (exclusive). Must be
   greater than ``--min-consensus``.

   Default: ``95``

.. option:: --seed <SEED>

   Seed for randomly subsampling columns. The same seed gives the same columns for any number of threads.

   Default: ``42``

.. option:: -t, --threads <THREADS>

   Number of threads, in the range [1, 1024].

   Default: ``1``


Files output
------------

* :ref:`[prefix].trimmed.fasta <files/alignment>`
* :ref:`[prefix].trimmed.mask <files/mask>`
* :ref:`[prefix].trimmed.filtered.mask <files/mask>`
* :ref:`[prefix].trimmed.markers.tsv <files/trimmed.markers.tsv>`
* :ref:`[prefix].trimmed.genomes.tsv <files/trimmed.genomes.tsv>`
* :ref:`gtrim.log <files/gtrim.log>`


Example
-------

Input
^^^^^

.. code-block:: bash

    gtrim trim -i trim_test.faa.gz -m trim_test_markers_info.tsv -o out_trim --msa-length 500 -t 4

Output
^^^^^^

.. code-block:: text

    [2026-11-14 03:16:57] INFO: gtrim v{{ gtrim_version }}
    [2026-11-14 03:16:57] INFO: gtrim trim -i trim_test.faa.gz -m trim_test_markers_info.tsv -o out_trim --msa-length 500 -t 4
    [2026-11-14 03:16:57] INFO: Using 4 threads.
    [2026-11-14 03:16:57] INFO: Starting MSA trimming workflow.
    [2026-11-14 03:16:57] INFO: Reading marker info file: trim_test_markers_info.tsv
    [2026-11-14 03:16:57] INFO: Read 10 markers with a concatenated length of 1792 columns.
    [2026-11-14 03:16:57] INFO: Reading alignment from file: trim_test.faa.gz
    [2026-11-14 03:16:57] INFO: Alignment contains 200 sequences and 1792 columns.
    [2026-11-14 03:16:57] INFO: Sequence type detected: amino acid
    [2026-11-14 03:16:57] INFO: Filtering columns with < 50% residues, < 25% consensus, or >= 95% consensus.
    [2026-11-14 03:16:57] INFO: Gappy columns: 73
    [2026-11-14 03:16:57] INFO: Conserved columns: 211
    [2026-11-14 03:16:57] INFO: Unconserved columns: 213
    [2026-11-14 03:16:57] INFO: Retained 1295 columns (min: 53, max: 280 per marker) after filtering.
    [2026-11-14 03:16:57] INFO: Randomly subsampling columns to obtain an MSA of ~500 columns (seed = 42).
    [2026-11-14 03:16:57] INFO: Sampling up to 50 columns from each marker.
    [2026-11-14 03:16:57] INFO: 0 markers had fewer than 50 columns available for sampling.
    [2026-11-14 03:16:57] INFO: Trimmed MSA has 500 columns.
    [2026-11-14 03:16:57] INFO: Writing trimmed alignment to file.
    [2026-11-14 03:16:57] INFO: Execution time: 0.01 seconds

This writes ``trim_test.trimmed.fasta.gz``, ``trim_test.trimmed.mask``, ``trim_test.trimmed.filtered.mask``,
``trim_test.trimmed.markers.tsv``, ``trim_test.trimmed.genomes.tsv`` and ``gtrim.log`` to ``out_trim``. The input
files are in ``tests/input_file/`` of the gTrim repository.

Example with a mask file
------------------------

``trim_test_witchi.faa.gz`` is ``trim_test.faa.gz`` with the 150 columns listed in ``trim_test_witchi.tsv`` removed.

.. code-block:: bash

    gtrim trim -i trim_test_witchi.faa.gz -m trim_test_markers_info.tsv --mask-file trim_test_witchi.tsv -o out_trim_mask --msa-length 500 -t 4

.. code-block:: text

    [2026-11-14 04:40:12] INFO: gtrim v{{ gtrim_version }}
    [2026-11-14 04:40:12] INFO: gtrim trim -i trim_test_witchi.faa.gz -m trim_test_markers_info.tsv --mask-file trim_test_witchi.tsv -o out_trim_mask --msa-length 500 -t 4
    [2026-11-14 04:40:12] INFO: Using 4 threads.
    [2026-11-14 04:40:12] INFO: Starting MSA trimming workflow.
    [2026-11-14 04:40:12] INFO: Reading marker info file: trim_test_markers_info.tsv
    [2026-11-14 04:40:12] INFO: Read 10 markers with a concatenated length of 1792 columns.
    [2026-11-14 04:40:12] INFO: Reading mask file: trim_test_witchi.tsv
    [2026-11-14 04:40:12] INFO: Mask file (WitChi format) retains 1642 of 1792 columns; 150 columns were removed from 10 markers.
    [2026-11-14 04:40:12] INFO: Reading alignment from file: trim_test_witchi.faa.gz
    [2026-11-14 04:40:12] INFO: Alignment contains 200 sequences and 1642 columns.
    [2026-11-14 04:40:12] INFO: Sequence type detected: amino acid
    [2026-11-14 04:40:12] INFO: Filtering columns with < 50% residues, < 25% consensus, or >= 95% consensus.
    [2026-11-14 04:40:12] INFO: Gappy columns: 68
    [2026-11-14 04:40:12] INFO: Conserved columns: 194
    [2026-11-14 04:40:12] INFO: Unconserved columns: 193
    [2026-11-14 04:40:12] INFO: Retained 1187 columns (min: 44, max: 254 per marker) after filtering.
    [2026-11-14 04:40:12] INFO: Randomly subsampling columns to obtain an MSA of ~500 columns (seed = 42).
    [2026-11-14 04:40:12] INFO: Sampling up to 51 columns from each marker.
    [2026-11-14 04:40:12] INFO: 1 markers had fewer than 51 columns available for sampling.
    [2026-11-14 04:40:12] INFO: Trimmed MSA has 503 columns.
    [2026-11-14 04:40:12] INFO: Writing trimmed alignment to file.
    [2026-11-14 04:40:12] INFO: Execution time: 0.01 seconds

The masks written to ``out_trim_mask`` have 1,792 columns, one for each column of ``trim_test.faa.gz``. The same
output is obtained with the binary mask ``trim_test_witchi.mask`` in place of ``trim_test_witchi.tsv``.
