.. _files/mask:

Mask (.mask)
============

A mask of the columns kept from an alignment. The file contains a single line with one character for each column of
the alignment: ``1`` if the column was kept and ``0`` if it was removed. A mask can be applied to an alignment with
:ref:`commands/mask`, and can be given to :ref:`commands/trim` with ``--mask-file`` to indicate the columns removed
from a concatenated alignment before trimming.

The following masks are written:

* ``[prefix].pruned.mask`` (:ref:`commands/prune-gap`): columns with a percentage of gaps less than or equal to
  ``--gap-threshold``.
* ``[prefix].trimmed.mask`` (:ref:`commands/trim`): columns in the trimmed alignment, i.e. after subsampling.
* ``[prefix].trimmed.filtered.mask`` (:ref:`commands/trim`): columns that passed the filters, before subsampling.

The masks written by :ref:`commands/trim` have one character for each column described by the marker info file. With
``--mask-file``, this is the full alignment before columns were removed, so columns removed by the mask file are
``0``.

Produced by
-----------

* :ref:`commands/prune-gap`
* :ref:`commands/trim`

Used by
-------

* :ref:`commands/mask`
* :ref:`commands/trim` (``--mask-file``)

Example
-------

``gap_test.pruned.mask`` for an alignment of 22 columns:

.. code-block:: text

    0000011111100000111111
