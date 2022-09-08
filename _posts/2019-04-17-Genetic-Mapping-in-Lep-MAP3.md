---
layout: post
title: "Genetic map construction in LepMap3"
tags: [Lep-MAP3, Genetic Mapping, Mapchart]
image: /image/lepmap/cropped_nih_geneticmap.jpg
share-img: /image/lepmap/cropped_nih_geneticmap.jpg
---

__Building genetic maps__ can be challenging and sometimes quite stressful, especially, when dealing with thousands or even millions of markers. In this post, I am hoping to help anyone who would like to get started to build a decent genetic map in an open software 
<a href="https://sourceforge.net/projects/lep-map3/"> Lep-MAP3 </a>, and finally, evaluating the accuarcy of the genetic map and plotting it.

# Table of contents


- [QC analysis](#-quality-control-analysis-)
- [Install LepMap3](#-installing-lep-map3-)
- [Running LepMap modules](#-running-lep-map3-)
- [Convert Phased data into Gnenotyes](#-30-converting-phased-output-data-from-ordermarkers2-to-genotypes-)
- [4way-Cross rQTL](#-covnvert-genotype-data-into-4-way-cross-rqtl-input-data-)
- [Validate Genetic Map](#-40-validate-the-genetic-map-by-conducting-qtl-analysis-)
- [Plot in MapChart](#-50-graphical-presentation-of-linkage-maps-in-mapchart-)


{: .box-note}
<i class="fa fa-commenting" aria-hidden="true"></i> **Note:** If you have an amplicon sequencing (`AmpSeq` or `rhAmpSeq`) haplotype data, you can convert the data into a psuedo VCF file using <a href="https://github.com/avinashkarn/analyze_amplicon/blob/master/haplotype_to_VCF.pl"> Haplotype to VCF </a> PERL script. 

<h2> Quality control analysis </h2>
Prior to building genetic maps - I strongly advise to perform QC analysis on your genetic data.
There are two QC tests that i usually perform: 1) `Multidimensional scaling (MDS)` and 2) `Check for Mendelian Error`

<p><a href="https://avikarn.com/2019-05-06-MDS/">Click here</a> for Multidimensional scaling (MDS) tutorial.</p>
<p><a href="https://avikarn.com/2019-08-04-Mendel/">Click here</a> for Mendelian Error Check tutorial. </p>

<br>

<h2>How to transfer files using FileZilla </h2>
Please watch below video to: download , install and configure `FileZilla`. It will show you how to upload files and folders to your server.

<iframe width="100%" height="300" src="https://www.youtube-nocookie.com/embed/adxmlHDim6c" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

<br>


<center><h2> Installing Lep-MAP3 </h2></center>
<strong>The Lep-MAP3 software is built in Linux and one has to have some experience in working in command-line environment.</strong>

Downloand and install Lep-Map3 on your computer following below steps:

<center><img src="/image/lepmap/install_lepmap.gif"></center>


<center><h2> Running Lep-MAP3 </h2></center>

The steps invloved in the genetic mapping process in Lep-MAP3 are shown in the flow chart below. 

<center><img src="/image/lepmap/lepmap_flow.png"></center>

<hr>

<h3>Step 1.1. File Preparation</h3>
<strong> Important - </strong> Correctly install the Lep-MAP3 software on your computer, and please make sure its the latest version. 
There are two files that are required as input files: 
<li>(1) <strong> genotype file</strong> in a <strong>VCF</strong> format, and </li>
<a href="/image/lepmap/sample_VCF.vcf" target="_blank">Download sample genotype file (VCF) here.</a>

<li>(2) <strong> pedigree</strong> file in <strong>.txt</strong> format. A snippet of the pedigree file showing relationship between all individuals in a family or population is shown below. It is important that the pedigree file is formatted exactly as shown in the below figure: </li>

<a href="/image/lepmap/sample_pedigree.txt" target="_blank">Download sample Pedigree file (.txt) here.</a>
<br>
<center><img src="/image/lepmap/ped.png"></center>
<hr>

<h3>Step 1.2. Parent Call</h3>
The parental genotypes are called using the `ParentCall2` module, using the below command:

```powershell
$ java -cp [path]/Lep-MAP3/bin ParentCall2 data = pedigree.txt  vcfFile = File.vcf > p.call
```
{: .box-note}
<i class="fa fa-commenting" aria-hidden="true"></i> **Note:** `path` is the directory where Lep-MAP3 is located on your computer.

<hr>
<center><img src="/image/lepmap/runlepmap3.gif"></center>

<hr>

<h3>Step 1.3. Filtering </h3>
This an optional step - However, One may use the `Filtering2` module to remove `non-informative markers` (Markers that are monomorphic or homozygous in both parents), and `distorted markers` (markers segregating in a non-Mendelian fashion) using the below command line:

```powershell
 $ java -cp /path/Lep-MAP3/bin Filtering2 data=p.call  removeNonInformative=1 dataTolerance=0.0000001  > p_fil.call
```

{: .box-note}
<i class="fa fa-commenting" aria-hidden="true"></i> **Note:** Use the parameter `removeNonInformative` to remove markers that are homozygous/monomorphic, and `dataTolerance` to remove distorted markers at __given__ p-value threshold.

<h3>Step 1.4. Separate Chromosomes </h3>

In this step, `SeparateChromosomes2` module is used to categorize markers into linkage groups (LGs) using the below command:

```powershell
$ java -cp /path/Lep-MAP3/bin  SeparateChromosomes2 data=p_fil.call lodLimit=5 > map.txt
```

__Note__: One can use parameters such as `lodLimit` and `theta` to split the linkage groups. 

One can check the number of markers in the in `map` file using the below command:
```powershell
$ sort map.txt|uniq -c|sort -n 
```
<br>

<h3>Step 1.5. Order Markers </h3>

In this step, markers separated into their corresponding linkage groups are ordered using `OrderMarkers2` module using the below command:

```powershell
$ java -cp /path/Lep-MAP3/bin  OrderMarkers2 data=p_fil.call map=map.txt > order.txt
```

	
One may use the parameter `sexAveraged`  to calculate sex-averaged map distances (by default male and female genetic maps are curated), also `numMergeIterations` parameter can be used to adjust number of iterations (by deafault its 6 iterations per linkage group). 

<h2> 2.0 Checking the accuracy of the marker order </h2>

If the physical positons of the markers in the curated genetic map curation are known, then one may use that information to evaluate the quality of the marker order in the genetic map, especially markers that <strong> inflate </strong> the chromosome length, by making a correlation plot of the genetic and physical positions of the markers for each chromosome or linkage group. 

<strong>Note:</strong> It is a common scenario to see the marker orders are flipped relative to their physical positions. There is nothing to panic about, one may fix it by manually sorting it.

Command to obtain the marker information using `cut`:

```powershell
cut -f1,2 p.call > cut_pcall.txt

```
{: .box-note}
<i class="fa fa-commenting" aria-hidden="true"></i> **Note:** Please make sure to use the p.call file that you used in the ordering step!

Follow the below steps to perform the correlation analysis:

<center><img src="/image/lepmap/order_correlation.gif"></center>

<br>
<center><img src="/image/lepmap/corr_geneticmap.png"></center>

<hr>

<h2> 3.0 Converting phased output data from OrderMarkers2 to genotypes </h2>

The phased data from `OrderMarkers2` step can be converted to fully informative "genotype" data by using `map2gentypes.awk` script and command below: 
<a href="/image/lepmap/map2genotypes.awk" target="_blank">Download the map2gentypes.awk script here.</a>
<br>

Next, run the `map2genotypes.awk` script by following the command shown below.
```powershell
	$ awk -vfullData=1 -f map2genotypes.awk order.txt > genotypes.txt
```

Snippet of the `map2gentypes.awk` output:

<center><img src="/image/lepmap/orderOutput.png"></center>
<hr>

One may convert the genotypes in <strong> 1 1 => A, 2 2 => B, 1 2 or 2 1 => H format </strong> (See below figure) in `MS Excel` using `Find` and `Replace` function, which can be then be loaded in `R/Qtl` for QTL mapping.

LepMap3 imputes and phase the genotype calls, therefore, A and B allele represent major and minor allele frequencies and it will change from parent/phase. They do not represent one specific parent and information depends on the parent and the phase of the marker. 

<center><img src="/image/lepmap/rqtlFormat.png">

</center>
<hr>

<h2> Covnvert genotype data into 4-way cross RQtl input data </h2>
In case both parents are heterozygotes, the cross is a 4way cross, also know as `AB x CD` phased output data from OrderMarkers2 __(1 1, 12, 21, 22)__ into Rqtl 4way code that __1, 2, 3 and 4__ as shown in the below table.  Also, please remember, in LepMap pedigree file: `male parent = 1 , female parent  = 2`,  and first digit of the phased genotypes is inherited from paternal parent and the second from maternal parent.

Read more on genotype file formatting here: https://github.com/kbroman/qtl/blob/main/man/read.cross.Rd

```console
LepMap	RQtl-4way-code	RQtl genotype
1 1	1	AC
1 2	2	BC
2 1	3	AD
2 2	4	BD
```
Finally import the converted the 4-waycross genetic map in RQTL using the below command: 
 ```powershell
 GenoData = read.cross(format = "csv", file = "geneticMap.csv", genotypes=NULL,
                     estimate.map = F, crosstype="4way")
```

<h2> 4.0 Validate the genetic map by conducting QTL analysis </h2>
<p>It is a good QC step to perform a QTL analysis of a well studied trait to check if expected QTL region is observed in the curated genetic mmap </p> 

<p><a href="https://avikarn.com/2019-04-23-composite-interval-mapping/">Click here</a> for QTL mapping tutorial.</p>

<center><img src="/image/lepmap/qtl.png"></center>


<hr>	


<h2> 5.0 Graphical presentation of linkage maps in Mapchart </h2>

Software Mapchart can be downloaded from below link.
https://www.wur.nl/en/show/Mapchart.htm

<p><a href="https://avikarn.com/2019-05-19-mapchart/">Click here</a> for Mapchart tutorial.</p>

<center><img src="/image/lepmap/mapChart.png"></center>


__Thank you__ for reading this tutorial. I really hope these steps will get you started in genetic map construction in Lep-MAP3. The key is to <strong> PRACTISE. </strong>. If you have any comments or suggestions, please let comment below or send me an email. 

Happy mapping !


<h3> Bibliography </h3>
Rastas, Pasi. "Lep-MAP3: robust linkage mapping even for low-coverage whole genome sequencing data." Bioinformatics 33.23 (2017): 3726-3732.
