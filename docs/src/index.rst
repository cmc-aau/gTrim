*****
gTrim
*****

.. warning::

   gTrim is a work in progress and in active development. Commands, options and output files may change
   between versions. Please contact Donovan Parks before using gTrim in your work.

**gTrim** is a command-line tool for removing columns from multiple sequence alignments (MSAs) before inferring a
phylogeny. It is written in Rust and designed for very large alignments, such as the concatenated marker gene
alignments of the `Genome Taxonomy Database (GTDB) <https://gtdb.ecogenomic.org>`_, with tens of thousands of columns
and hundreds of thousands of sequences.

gTrim can:

* remove columns with too many gaps (:ref:`commands/prune-gap`);
* trim a concatenated marker gene alignment by removing gappy, highly conserved and poorly conserved columns, and
  optionally subsample the same number of columns from each marker (:ref:`commands/trim`), as done for the GTDB
  reference trees;
* apply a mask of columns to keep to an alignment (:ref:`commands/mask`);
* report gap and residue statistics for each column of an alignment (:ref:`commands/column-stats`).

gTrim is open source and released under the
`GNU General Public License (Version 3) <https://www.gnu.org/licenses/gpl-3.0.en.html>`_.

Please post questions and issues related to gTrim on the
`Issues <https://github.com/cmc-aau/gTrim/issues>`_ section of the GitHub repository.


Running gTrim
=============

1. Build gTrim from source (:ref:`installing`).

2. Read the documentation for each of the :ref:`commands`, or view the program help menu: ``gtrim -h``

* Note: individual help can be accessed via the specific command, e.g.: ``gtrim trim -h``

For example, to trim a concatenated marker gene alignment:

.. code-block:: bash

    gtrim trim -i bac120_msa.faa.gz -m bac120_markers_info.tsv -o trim_out --threads 16


Citing gTrim
============

gTrim has not yet been published. If you use gTrim, please cite the GitHub repository
(https://github.com/cmc-aau/gTrim) and the version you used. If you use the :ref:`commands/trim` command, please
also cite the GTDB publications that describe this trimming approach (see :ref:`references`).


.. toctree::
   :caption: Getting started
   :maxdepth: 1

   Home <self>
   installing
   faq


.. toctree::
   :caption: Running gTrim
   :maxdepth: 1

   commands/index
   files/index


.. toctree::
   :caption: About
   :maxdepth: 1

   about/how_trim_works
   changelog
   references
