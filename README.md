Getting started
================================

git clone https://github.com/shiquan/bedutils.git

cd bedutils; make

Introduction
----------
Bedutils is a software package for handling three columns bed
files. It consist of six algrithm till now. See roadmap.txt for
my further plans.

Usage
------
merge several bed files:

	bedutils merge [option] <in1.bed> [in2.bed ...]

find the uniq regions of several bed files:

	bedutils uniq [option] <in1.bed> <in2.bed> [...]

find the difference between the first bed file and others:

	bedutils diff <in1.bed> <in2.bed> [...]

trim a bed file:

	bedutils trim -r <right_offset> -l <left_offset> <in1.bed>

compare the first bed file with others:

	bedutils comp <in1.bed> <in2.bed> [in3.bed ..]

count the overlap regions of bed files:

	bedutils count <in1.bed> [in2.bed...]


Global options
------------

**-o**  set the output file

**-r**  the right offset, default is flank

**-l**  the left offset, default is flank


Demo
-----

## format in.bed, merge the overlap regions into one region

	bedutils merge in.bed

## flank 100bp of in.bed

	bedutils merge -r 100 -l 100 in.bed

## find the uniq regions of in1.bed and in2.bed

	bedutils uniq in1.bed in2.bed

