.. _commands:

####################
  Commands
####################

gTrim has the following commands:

.. toctree::
   :maxdepth: 1

   prune-gap
   trim
   mask
   column-stats

Every command reads a single alignment in FASTA format (``-i/--alignment-file``), which can be gzip-compressed, and
writes its output to a directory (``-o/--out-dir``) that is created if it does not exist. Input files must end in
``.faa``, ``.fna``, ``.fasta`` or ``.gz``.

Output files are named after the input file, with the ``.gz``, ``.faa``, ``.fna`` or ``.fasta`` extension removed and
a word for the command added (see :ref:`files`). Parameter values are not included in output file names, so runs with
different parameters should use different output directories. If the input alignment is gzip-compressed, the output
alignment is also compressed.

Whether the alignment contains nucleotides or amino acids is detected automatically from the first 100,000
characters. Alignments containing characters other than the nucleotides ``ACGT`` or the 20 standard amino acids
(ignoring ``N``, ``X``, ``?``, ``-`` and ``.``) are treated as amino acid alignments and a warning is reported.
