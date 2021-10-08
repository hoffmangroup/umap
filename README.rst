.. image:: https://readthedocs.org/projects/umap/badge/?version=latest


Umap and Bismap: quantifying genome and methylome mappability
=============================================================


Introduction
------------

**The free umap software package efficiently identifies uniquely mappable regions of any genome.
Its Bismap extension identifies mappability of the bisulfite converted genome (methylome).**

the mappability of a genome for a given read length *k*.
First, it generates all possible *k*-mers of the genome.
Second, it maps these unique *k*-mers to the genome with `Bowtie <http://bowtie-bio.sourceforge.net/index.shtml>`_ version 1.1.0.
Third, Umap marks the start position of each *k*-mer that aligns to only one region in the genome.
Umap repeats these steps for a range of different *k*-mers and stores the data of each chromosome
in a binary vector *X* with the same length as the chromosome's sequence.
For read length *k*, :math:`X_i = 1` means that the sequence starting at :math:`X_i` and ending
at :math:`X_{i+k}` is uniquely mappable on the forward strand.
Since we align to both strands of the genome, the reverse complement of this same sequence
starting at :math:`X_{i+k}` on the reverse strand is also uniquely mappable.
:math:`X_i = 0` means that the sequence starting at :math:`X_i` and ending at
:math:`X_{i+k}` can be mapped to at least two different regions in the genome. 


Mappability of the bisulfite-converted genome
---------------------------------------------

To identify the single-read mappability of a bisulfite-converted genome,
we create two altered genome sequences.
In the first sequence, we convert all cytosines to thymine (C :math:`\rightarrow` T).
In the other sequence we convert all guanines to adenine (G :math:`\rightarrow` A).
Our approach follows those of `Bismark <http://www.bioinformatics.babraham.ac.uk/projects/bismark/>`_
and `BWA-meth <https://github.com/brentp/bwa-meth>`_.
We convert the genome sequence this way because bisulfite
treatment converts un-methylated cytosine to uracil which is read as thymine.
Similarly the guanine that is base-pairing with the un-methylated cytosine in
reverse strand converts to adenine. However these two conversions never occur
at the same time on the same read. We identify the uniquely mappable regions
of these two genomes separately, and then combine the data to represent the
single-read mappability of the forward and reverse strands in the bisulfite-converted genome.


Bismap requires special handling of reverse complementation of C :math:`\rightarrow` T
or G :math:`\rightarrow` A converted genomes.
Conversion of C :math:`\rightarrow` T on the sequence AATTCCGG produces AATT **TT** GG.
In the Bowtie index, the reverse complement would be CCAAAATT.
However for the purpose of identifying the mappability of the bisulfite-converted genome,
we expect the reverse complement to be TTGGAA **TT**. The reason is that both forward and reverse
strands undergo bisulfite treatment simultaneously. There is no DNA replication after bisulfite treatment.
To handle this issue, Bismap creates its own reverse complemented chromosomes and suppresses Bowtie's usual reverse complement mapping.

Umap and Bismap each take approximately 200 CPU hours to run for a given read length. This can be parallelized in a computing cluster over 400 cores to take only 30 minutes. 


Measures of mappability
:::::::::::::::::::::::

Umap efficiently identifies the single-read mappability of any genome for a
range of sequencing read lengths. The single-read mappability of a genomic
region is a fraction of that region which overlaps with at least one uniquely
mappable *k*-mer. The Bismap extension of Umap produces the single-read mappability
of a bisulfite-converted genome. Both Umap and Bismap produce an integer vector for
each chromosome that efficiently defines the mappability for any region and can be
converted to a browser extensible data (BED) file. In addition to single-read mappability,
we can measure the mappability of a genomic region by another approach. To quantify
the single-read mappability of a given genomic region, we measure the fraction of
potential uniquely mappable reads in that region. A region, however, can have 100% single-read
mappability, but in practice require a high coverage sequencing to properly identify that region.
For example, a 1 kbp region with 100% single-read mappability can be mappable due to a
minimum of 10 unique 100-mers that none of them overlap or a maximum of 1100 unique 100-mers
that highly overlap. Therefore, we define the multi-read mappability, the probability that a
randomly selected read of length k in a given region is uniquely mappable. For the genomic
region :math:`G_{i:j}` starting at *i* and ending at *j*, there are :math:`j - i + k + 1`
*k*-mers that overlap with :math:`G_{i:j}`.
The multi-read mappability of :math:`G_{i:j}` is the fraction of those *k*-mers that are uniquely mappable. 


News
----

* Version 1.2.0: Fixed the issue with 1-based inclusive interval BED files. As of this
  version, the BED files are 0-based and they don't contain overlapping intervals.
* Version 1.1.1: In addition to wiggle, you can create bedGraph of multi-read mappability.
  New expectation: order of chromosomes in chromosome-size file must match genome FASTA.
  In this version, *unify_bowtie* filters out any unexpected unique mapping of a *k*-mer
  to unintended chromosomes.
  This version released on November 24th 2017.
* Version 1.1.0: We released Umap version 1.1.0 on May 2nd 2017.
  As of this version, uint8 non-zero values correspond to the start position of
  a k-mer that starts at that position and ends in k nucleotides downstream, and is unique.
  Before this version, non-zero values in the uint8 files corresponded to all of nucleotides
  covered by that unique k-mer (not only the start position).
  Any multi-read mappability measure based on previous versions,
  therefore, is over-estimated.
* Version 1.0.0: We released Umap version 1.0.0 on December 20th 2016.
  *Deprecated, do not use!*



Quick start
-----------

We have tested Umap installation on a CentOS system using python 2.7.11.
Bismap requires numpy and pandas and it uses other python modules such as:

* gzip
* os
* re
* subprocess

Umap uses mercurial version control. Make sure that mercurial (hg) is installed.
Download Umap to the directory of your python packages using::

    git clone https://github.com/hoffmangroup/umap.git
    cd umap
    python setup.py install


Now we will run a test using a wrapper in the umap directory called ubismap.py
and a toy genome stored under umap/data ::

    cd umap
    python ubismap.py -h
    python ubismap.py data/genome.fa data/chrsize.tsv data/TestGenomeMappability all.q $BOWTIEDIR/bowtie-build --kmer 8 12 -write_script test_run.sh
    sh test_run.sh


**Important: Order of chromosomes in your genome FASTA must match order
of chromosomes in your chromosome-size file. Chromosome-size file should not
contain header (2 column file with chromosome name and length of the chromosome).**

    
The scripts that are produced by **ubismap.py** assume that you are using a Sun Grid Engine computing cluster.
You can use parameters of this script to adjust it to your own system. You may need to manually edit this file
because many of the SGE settings are very different than other computing clusters.
However, all of the Umap modules accept *-job_id* which allows you to use the modules without a cluster or if
your cluster does not support job arrays.

Basically, the job array saves an environmental variable (defined by (*-var_id*) that Umap uses for paralellizing processes.
You can run the modules in a for loop and set the *-job_id* manually.
For example, in order to find *k*-mers of a genome with 10 million base pairs, the get_kmers
and run_bowtie modules each need to be executed with -job_ids ranging between 1 to 10.


Get *k*-mers
------------

.. module:: umap
    :platform: Unix

.. currentmodule:: get_kmers
.. autoclass:: GetKmers


The first step of finding the mappability of the genome for a given read length,
is to create all of the possible sequences of the genome with the read length of interest.
GetKmers requires an index reference and an index ID to do this step for one chunk of
the genome at a time. The index ID can be specified by -job_id, or if GetKmers is
submitted as a job array, GetKmers will use the variable name set by -var_id to obtain
the environmental variable for the job array ID.



**Important: The chromosome names in the fasta file should not contain underscore.
Underscore is used in Bismap to differentiate reverse complement chromosomes.**



Run Bowtie
----------

.. currentmodule:: run_bowtie
.. autoclass:: BowtieWrapper


When all of the possible *k*-mers of the genome are created, BowtieWrapper keeps record
of all of the k-mers that are unique. Similar to GetKmers, the parallel process relies on
an index reference, or -job_id.


Merge bowtie outputs
--------------------


.. currentmodule:: unify_bowtie
.. autoclass:: UnifyBowtie


For each *k*-mer file, there will be a bowtie file with information of unique *k*-mers.
For each chromosome, a binary vector with length of the chromosome will be created.
If a *k*-mer starting at a given position is unique, the value will be 1.


**Important: Name of you chromosome should start with "chr" and not contain underscore.
Order of chromosomes in your genome FASTA must match order
of chromosomes in your chromosome size file.**


Merge data of various *k*-mers
------------------------------


.. currentmodule:: combine_umaps
.. autoclass:: CombineUmaps

If the above steps are repeated for various *k*-mers, these data can be efficiently
stored in a numeric vector for each chromosome. For each position, the length of the smallest
*k*-mer that uniquely maps to that position is specified.


Convert numeric vectors to BED and Wiggle
-----------------------------------------


.. currentmodule:: uint8_to_bed
.. autoclass:: Int8Handler
.. method:: Int8Handler.write_beds
.. method:: Int8Handler.write_as_wig

To visualize binary and numeric vectors that are produced by Umap, you can use Int8Handler.


Special instructions for Bismap
-------------------------------

For Bismap, you must create the mappability of C :math:`\rightarrow` T
and G :math:`\rightarrow` A by specifying -Bismap, -C2T, and -G2A each time.
After creating these two mappability, for each genome you will have a
kmers/globalmap_k<min>tok<max> folder with normal and reverse complemented chromosome mappabilities.


The first step is to merge the mappability of normal and reverse complemented chromosomes for each
genome. This can be done through uint8_to_bed.py ::

    python uint8_to_bed.py <MergedKmerDir> <OutDirC2T> <LabelForOutputFiles> -C2T -chrsize_path <ChromSizeFile> -WriteUnique
    python uint8_to_bed.py <MergedKmerDir> <OutDirG2A> <LabelForOutputFiles> -G2A -chrsize_path <ChromSizeFile> -WriteUnique


When you have created the merged mappability of each chromosome once for C :math:`\rightarrow` T genome
and once for G :math:`\rightarrow` A genome, you should use combine_umaps.py and specify -kmer_dir_2 ::

    qsub <...> -t 1-<NumberOfChromosomes> python combine_umaps.py <OutDirC2T> <ChromSizePath> -out_dir <OutDir> -kmer_dir_2 <OutDirG2A>


Requesting Genomes
------------------

In case you need these data for other genomes and do not have access to a Sun Grid Engine computing cluster,
we may accept to do this for you.
Please contact the software maintainer by email and we will do our best to assist you as soon as possible.



Contact, support and questions
------------------------------

For support of Umap, please user our `mailing list <https://groups.google.com/forum/#!forum/ubismap>`_.
Specifically, if you want to report a bug or request a feature,
please do so using
the `Umap issue tracker <https://github.com/hoffmangroup/umap/issues>`_.
We are interested in all comments on the package,
and the ease of use of installation and documentation.


Credits
-------


Umap was originally developed by Anshul Kundaje and was written in MATLAB.
The original repository is available `here <https://sites.google.com/site/anshulkundaje/projects/mappability>`_.
The version of Umap that is available in this repository, is a python reimplementation of the original Umap, and
is initially written by Mehran Karimzadeh.
