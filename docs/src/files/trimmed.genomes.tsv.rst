.. _files/trimmed.genomes.tsv:

trimmed.genomes.tsv
===================

A tab-separated file reporting how complete each sequence is in the trimmed alignment, written to
``[prefix].trimmed.genomes.tsv``. This can be used to identify genomes with few residues in the trimmed alignment,
which may be poorly placed in a phylogeny. Residues are counted as the ASCII letters, so gaps and characters such as
``*`` and ``?`` are not counted.

The columns in this file are as follows:

* ``genome_id``: identifier of the sequence.
* ``msa_length``: number of columns in the trimmed alignment.
* ``residues``: number of residues of the sequence in the trimmed alignment.
* ``residues_perc``: percentage of columns of the trimmed alignment in which the sequence has a residue.

Produced by
-----------

* :ref:`commands/trim`

Example
-------

.. code-block:: text

    genome_id	msa_length	residues	residues_perc
    G000176175	500	446	89.20
    G001464235	500	496	99.20
    G002406055	500	487	97.40
