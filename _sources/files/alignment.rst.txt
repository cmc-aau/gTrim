.. _files/alignment:

Alignment (.fasta)
==================

The alignment with the removed columns excluded, in FASTA format. Each sequence is written on a single line, in the
same order and with the same identifiers as in the input alignment. If the input alignment is gzip-compressed, the
output alignment is also compressed and ``.gz`` is appended to its name.

The file is named:

* ``[prefix].pruned.fasta`` by :ref:`commands/prune-gap`;
* ``[prefix].trimmed.fasta`` by :ref:`commands/trim`;
* ``[prefix].masked.fasta`` by :ref:`commands/mask`.

Produced by
-----------

* :ref:`commands/prune-gap`
* :ref:`commands/trim`
* :ref:`commands/mask`

Example
-------

The first two sequences of ``gap_test.fna.gz``:

.. code-block:: text

    >seq1
    ----------G-ACGTACGTAC
    >seq2
    ---------CG--CGTACGTAC

and of ``gap_test.pruned.fasta.gz``, after removing columns with more than 50% gaps:

.. code-block:: text

    >seq1
    -----GACGTAC
    >seq2
    ----CGACGTAC
