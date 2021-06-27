# tspOrder2
Using traveling salesman problem (TSP) software Concorde or LKH2 to fine-order markers in a genetic map

# Introduction

I have been using [ASMap](https://cran.r-project.org/web/packages/ASMap/index.html) R package, which uses [MSTMap](http://mstmap.org/), to construct genetic maps for wheat (14 or 21 pairs of chromosomes) that are genotyped with 9K or 90K SNPs. It can make very good genetic maps for thousands of markers in a short time. But I have no ways to check whether the marker order is optimal. Later I found another R package [TSPmap](https://github.com/mckaylab/TSPmap), a tool making use of traveling salesperson problem solvers in the efficient and accurate construction of high-density genetic linkage maps. I tested it on wheat with 90K markers, but it did not work well. In wheat at least, it is common to get partial linkage groups due to either marker distribution or closeness of parents. So I often need to try different grouping stringencies and use a higher stringency to break them into more groups in the chromosome numbers, then merge linkage groups from the same chromosomes. **TSPmap** needs you to give a chromosome number, which is not practical for wheat. Although I cannot use **TSPmap** to make genetic maps from scratch, the "tspOrder" function is a very good tool to check the order of a constructed map.

**TSPmap** is not in the R package repository, so it is not very easy to install it. **tspOrder** require another R package **TSP**. To make it easier to use and no requires other packages, I took some functions out of **TSPmap** and added some new functions to make an R script that is enough to run **tspOrder**. 

In addition, I modified the **tspOrder** function so it can also use [LKH2](http://akira.ruc.dk/~keld/research/LKH/) to check marker orders. LKH2 can get the same shortest map length as [concorde](https://www.math.uwaterloo.ca/tsp/concorde.html), a recommended software by **TSPmap**, but LKH2 is a little faster than **concorde**.

Another problem for high-density maps is that multiple optimal orders can get the shortest maps. I would like the one that is more consistent with the physical maps. Using random seeds, tsp software can give different optimal orders (same length but a little different orders in some regions). So I just run multiple times of **tspOrder** to improve regions that can have multiple orders. Those regions usually include a group of markers with equal recombination frequency to each other, so their order can change.
Multiple runs of **tspOrder** are not the most resources-efficient way, but it is the easiest way to do. I may improve the code later when I get time.

# Usage
Try the [examples.R](examples.R) file.

These R functions only help you check whether your marker order is optimal. So you should already have genetic maps made for your population.

If you want to improve your marker order to be more consistent with the physical map, please add physical positions (chromosome names and position in Mb, separated by "_") to your marker names, such as "IWB3087_1A_9.16". The physical chromosome names should be the same as your genetic chromosome names, do not use "ch1A" for the physical map, but "1A" for genetic maps.

# Software needed
[Concorde](https://www.math.uwaterloo.ca/tsp/concorde/downloads/downloads.htm) OR [LKH2](http://akira.ruc.dk/~keld/research/LKH/) are needed for these R scripts. Although my scripts use a GPL3 license, [Concorde](https://www.math.uwaterloo.ca/tsp/concorde/downloads/downloads.htm) and [LKH2](http://akira.ruc.dk/~keld/research/LKH/) are for **academic use only**. Please check their website for details. I downloaded the Linux and Windows version of [Concorde](https://www.math.uwaterloo.ca/tsp/concorde/downloads/downloads.htm) for my only use. The windows version needs the **cygwin1.dll** (32-bit, no need to install cygwin-32bit) and I have put them together. LKH2 windows version can be [downloaded](http://akira.ruc.dk/~keld/research/LKH/LKH-2.exe) from their website. I compiled the Linux version for my use only.