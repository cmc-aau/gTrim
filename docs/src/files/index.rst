.. _files:

####################
  Files
####################

Output files are named after the input alignment. The ``[prefix]`` is the name of the input file with the ``.gz``,
``.faa``, ``.fna`` or ``.fasta`` extension removed, followed by a word for the command: ``.pruned``
(:ref:`commands/prune-gap`), ``.trimmed`` (:ref:`commands/trim`), ``.masked`` (:ref:`commands/mask`) or
``.column_stats`` (:ref:`commands/column-stats`). For example, :ref:`commands/trim` on
``gtdb_r226_bac120_concatenated.faa.gz`` writes ``gtdb_r226_bac120_concatenated.trimmed.fasta.gz``.

Input files
-----------

.. toctree::
   :maxdepth: 1

   markers_info.tsv
   witchi.tsv

Output files
------------

.. toctree::
   :maxdepth: 1

   alignment
   mask
   pruned.tsv
   trimmed.markers.tsv
   trimmed.genomes.tsv
   column_stats.tsv
   gtrim.log
