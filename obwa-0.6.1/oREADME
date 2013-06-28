Introduction
------------

This README file is about bwa add-ons 'obwase' and 'obwape', further referred as 'the programs'.

The programs attempt to optimise parts of bwa by merging the functionalities
of aligning ('bwa aln') and alignment generation ('bwa samse' or 'bwa sampe')
to a single step. Traditionally, alignment and alignment generation
steps in bwa are either executed in pipe or as separate consequentive steps.
The aim in this optimisation is to reduce the footprint of RAM and disk I/O.
The RAM is reduced by not having to read two copies of the indexes while running
alignment and alignment generation steps for a single end alignment in a pipe.
The disk I/O is reduced by not writing the alignment results (sai files) to disk
in consecutively executed alignment and alignment generation steps both in single end
and in paired end alignment.

The programs are expected to behave and produce the same output as 'bwa aln',
'bwa samse' and 'bwa sampe' in combined. The output sam file is expected to be
identical to the one created by bwa.

The source code for the programs is provided together with the original source code of the bwa.
The bwa source code is copied unchanged from the author(s) of the bwa (http://bio-bwa.sourceforge.net/)
except for modified Makefile and files added to create obwase and obwape. None of the original files
except the Makefile of the bwa has been edited. The Makefile has not been changed for making the bwa,
only additions for making obwase and obwape have been added. All of the files added to the
bwa source code and none of the original files start with prefix o. The versions
of the programs can be found by executing the programs without any parameters.
The versioning of obwa follows the versioning of bwa. For example: obwa-0.6.1 copies the code
unchanged from bwa-0.6.1.

Comments and suggestions (but no requests for help in any form) are welcome to <margus.lukk@cruk.cam.ac.uk>.

The code is released under the licence of GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007,
as is the original code of bwa.


Installation and instructions
-----------------------------

While in directory of the source code, type:

make

This compiles and links bwa as well as obwase and obwape.

For instructions of how to use the programs, for information about passing the 'bwa aln',
'bwa samse' or 'bwa sampe' options to the programs; and for examples, the programs
should be run without any parameters.
