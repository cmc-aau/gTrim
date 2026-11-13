.. _installing:

Installing gTrim
================

gTrim is not yet available from Bioconda or crates.io and must be built from source.


Building from source
--------------------

gTrim requires a `Rust toolchain <https://www.rust-lang.org/tools/install>`_ that supports the 2024 edition of Rust
(Rust 1.85 or later). Rust can be installed with ``rustup``:

.. code-block:: bash

    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

Clone the repository and build an optimised binary:

.. code-block:: bash

    git clone https://github.com/cmc-aau/gTrim.git
    cd gTrim
    cargo build --release

The binary is written to ``target/release/gtrim``. Copy it to a directory on your ``PATH``, or install it into
``~/.cargo/bin`` with:

.. code-block:: bash

    cargo install --path .

Check the installation:

.. code-block:: bash

    gtrim --version


Testing
-------

The unit and integration tests can be run from the repository with:

.. code-block:: bash

    cargo test

The integration tests run each command on the small alignments in ``tests/input_file/`` and compare the output
to the expected results in ``tests/ground_truth/``.


Hardware requirements
---------------------

gTrim reads the entire alignment into memory, using one byte per residue. For example, an alignment of
10,000 sequences and 41,084 columns needs about 0.4 GB of memory. The ``prune-gap`` and ``trim`` commands
calculate column statistics in parallel; use ``--threads`` to set the number of threads.

On an alignment of 10,000 sequences and 41,084 columns, each command takes 2 to 3 seconds with 16 threads. Most of
this time is spent decompressing the gzipped input alignment.


Rust libraries
--------------

gTrim uses the following libraries, which are downloaded and compiled automatically by ``cargo``:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Library
     - Purpose
   * - `needletail <https://github.com/onecodex/needletail>`_
     - Reading FASTA files.
   * - `ndarray <https://github.com/rust-ndarray/ndarray>`_
     - Storing the alignment as a two-dimensional array.
   * - `rayon <https://github.com/rayon-rs/rayon>`_
     - Parallel computation of column statistics.
   * - `gzp <https://github.com/sstadick/gzp>`_
     - Multithreaded gzip compression of output alignments.
   * - `rand <https://github.com/rust-random/rand>`_ and `rand_chacha <https://github.com/rust-random/rand>`_
     - Seeded random subsampling of columns.
   * - `clap <https://github.com/clap-rs/clap>`_
     - Command line interface.
   * - `log4rs <https://github.com/estk/log4rs>`_
     - Logging to the console and ``gtrim.log``.
