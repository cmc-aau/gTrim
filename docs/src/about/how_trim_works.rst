.. _about/how_trim_works:

How trim works
==============

:ref:`commands/trim` removes columns from a concatenated marker gene alignment that are unlikely to be informative
for inferring a phylogeny, and can then subsample the same number of columns from each marker. It reimplements the
``msa_filter`` command of `GenomeTreeTk <https://github.com/dparks1134/GenomeTreeTk>`_, which has been used to trim
the concatenated marker gene alignments of the Genome Taxonomy Database (GTDB).


Columns removed before trimming
-------------------------------

The marker info file gives the length of each marker in the full concatenated alignment. If columns were removed from
this alignment before trimming, for example by `WitChi <https://github.com/stephkoest/witchi>`_, the markers no
longer start where the marker info file says they do. ``--mask-file`` gives the columns that were removed, either as a
:ref:`WitChi pruned TSV <files/witchi.tsv>` or as a :ref:`binary mask <files/mask>`. Each marker is then reduced by the
number of its columns that were removed, and the columns of the input alignment are assigned to markers in order.

``trim`` checks that the inputs line up before doing any filtering:

* Without ``--mask-file``, the marker lengths must sum to the length of the alignment.
* With ``--mask-file``, a binary mask must have one character for each column described by the marker info file, the
  indices in a WitChi file must be less than the summed marker length and be listed only once, and the number of
  columns retained by the mask must equal the length of the alignment.

The output masks are over the full alignment: columns removed by ``--mask-file`` are ``0``, and the
``masked_columns`` column of :ref:`files/trimmed.markers.tsv` gives the number removed from each marker.


1. Filtering columns
--------------------

Each column of the alignment is classified as follows, using the ``-``, ``.``, ``_`` and ``*`` characters as gaps:

* **Gappy:** more than (100 − ``--min-perc-aa``)% of sequences have a gap in the column, or every sequence does.
  With the default of 50, a column is gappy if more than half of the sequences have a gap.
* Otherwise, the **consensus** of the column is the number of sequences with the most common residue divided by the
  number of sequences without a gap. A column is:

  * **unconserved** if its consensus is below ``--min-consensus`` (default: 25%);
  * **conserved** if its consensus is at or above ``--max-consensus`` (default: 95%);
  * **retained** otherwise.

Only retained columns pass the filters. They are written to ``[prefix].trimmed.filtered.mask``, and the number of
columns in each class is reported for each marker in :ref:`files/trimmed.markers.tsv`.


2. Subsampling columns
----------------------

To give each marker a similar influence on the phylogeny, the same number of columns, *k*, is randomly sampled from
the retained columns of each marker. *k* is the smallest number for which the sum over all markers of
min(retained columns, *k*) is at least ``--msa-length``. Markers with fewer than *k* retained columns contribute all
of their retained columns. If fewer than ``--msa-length`` columns are retained, every retained column is kept and no
subsampling is done; this is the case by default, as ``--msa-length`` is 1,000,000.

Columns are sampled with a ChaCha8 random number generator seeded with ``--seed`` (default: 42), using a single
generator for the markers in the order they are listed in the marker info file. The sampled columns therefore depend
only on the seed and the input, not on the number of threads. The sampled columns are written to
``[prefix].trimmed.mask`` and form the trimmed alignment.


Differences from msa_filter
---------------------------

``trim`` produces the same filtered columns as ``msa_filter`` and samples the same number of columns from each
marker. It differs in the following ways:

* Sampling is seeded, so results are reproducible. ``msa_filter`` is not seeded, so the specific columns sampled
  from each marker differ between the two programs and between runs of ``msa_filter``.
* The default ``--msa-length`` is 1,000,000, so by default no columns are subsampled. The default of ``msa_filter``
  is 5,000.
* If ``--msa-length`` cannot be reached, every retained column is used. ``msa_filter`` samples at most 999 columns
  from each marker.
* A run in which no columns pass the filters is reported as an error.
* ``--mask-file`` accepts a binary mask as well as a WitChi pruned TSV.
* ``msa_filter`` assigns a removed column that is the first column of a marker (other than the first marker) to the
  preceding marker, so the boundary between the two markers is shifted by one column in the input alignment. ``trim``
  assigns each removed column to the marker it belongs to, and reports in the log how many removed columns are the
  first column of a marker. When there are any, the filtered columns are the same, but the per-marker statistics and
  the number of columns sampled from the affected markers can differ from ``msa_filter``.

The ``--min-perc-aa`` option filters columns, not genomes, despite the wording of the ``msa_filter`` help text.


Validation
----------

``trim`` was compared to ``msa_filter`` on an alignment of 5,000 GTDB genomes using the bac120 markers of GTDB
release 226, with ``--msa-length 5000`` and with
``--msa-length 12000 --min-perc-aa 70 --min-consensus 10 --max-consensus 90``. In both cases, the filtered masks and
the number of retained and sampled columns of every marker were identical. ``trim`` took 2 seconds compared to 26
seconds for ``msa_filter``.

With ``--mask-file``, ``trim`` was compared to ``msa_filter`` given a WitChi pruned TSV removing 2,500 random
columns from this alignment, none of which was the first column of a marker, with ``--msa-length 5000``: the filtered
masks, per-marker statistics and number of columns sampled from each marker were identical, and a binary mask of the
same columns gave identical output. With a second WitChi file removing 2,504 columns, 9 of which were the first column of a
marker, the filtered masks were still identical but the statistics of the affected markers differed, as described
above.

The integration tests of gTrim include ``trim`` tests on the first 200 genomes and first 10 bac120 markers of this
alignment, with and without a mask file, whose expected output was checked against ``msa_filter``.
