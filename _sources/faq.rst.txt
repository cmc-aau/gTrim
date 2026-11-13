.. _faq:

FAQ
===


Which alignment formats are supported?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Alignments in FASTA format with a ``.faa``, ``.fna``, ``.fasta`` or ``.gz`` extension. Gzip-compressed alignments are
read directly, and the output alignment is then also gzip-compressed. All sequences must have the same length.


Does gTrim work on nucleotide alignments?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Yes. :ref:`commands/prune-gap`, :ref:`commands/mask` and :ref:`commands/column-stats` can be used with nucleotide or
amino acid alignments. :ref:`commands/trim` is designed for concatenated protein alignments, but its filters do not
depend on the type of sequence.


Which characters are treated as gaps?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:ref:`commands/prune-gap` and :ref:`commands/column-stats` treat ``-`` and ``.`` as gaps. :ref:`commands/trim` also
treats ``_`` and ``*`` as gaps, to match ``msa_filter`` of GenomeTreeTk.


Why are my output files overwritten?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Output file names depend only on the input file name and the command, not on the parameters used. Running a command
with different parameters in the same output directory overwrites the previous results. Use a different output
directory for each set of parameters. The parameters of each run are recorded in :ref:`files/gtrim.log`.


How do I trim an alignment the same way as GTDB?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run :ref:`commands/trim` with the marker info file for the alignment's marker set. To subsample a 5,000 column
alignment as done by the default settings of ``msa_filter``, use ``--msa-length 5000``.


How do I apply the same trimming to another alignment?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use :ref:`commands/mask` with the ``.mask`` file from the first run. The other alignment must have the same columns,
e.g. additional genomes aligned to the same markers.


Is the output reproducible?
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Yes. Random subsampling in :ref:`commands/trim` is seeded with ``--seed``, and the results do not depend on the
number of threads.


How much memory does gTrim need?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

About one byte per character of the alignment, i.e. the number of sequences multiplied by the number of columns.
See :ref:`installing`.


How do I cite gTrim?
~~~~~~~~~~~~~~~~~~~~

gTrim has not yet been published. See :ref:`references`.
