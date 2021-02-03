# bsa
Bulked-Segregant Analysis using vcf file with or without parents(you can call it vcfbsa), you don't need to construct parent's reference to call snp and don't need to polarize alleles first to make direction of delta-snpindex meaningful.

It's memory effecient and very fast, you can complete the analysis in about one hour  with only < 4Gb RAM when there are ~10 millions markers in your vcf file.

The most memory cosuming step is the slidding windows step, because the slidewindow.pl script need to read all avalible sites in index/ed/gst/fet file into RAM.

Can also be used for mutmap analysis(just give the same name for -b1 and -p1 when calcualte index).

# Dependences
You need to install these dependence perl modules from cpan first

use Cwd

use Getopt::Long

use Data::Dumper

use File::Basename

use Text::NSP::Measures::2D::Fisher::twotailed

use FindBin


R (only >= v3.60 were tested) is needed, and [gtools](https://cran.r-project.org/web/packages/gtools/index.html) package is needed if you want to use point_line_plot.pl script.

# BSA calculations
This script can be used to calculate snp-index/delta snp-index(index), g-statistic(gst), euclidean distance(ed) and Fisher-exact test(fet) and some other methods (coming soon).

For bi-parents populations (f2, ril/dh, bc), only aaxbb sites are used when two parents are avalible, only aax?? or ??xbb sites are used when only one parent is avalible, and all sites are used when two parents are absence. For outcross population (f1/cp), only lmxll nnxnp and hkxhk sites are used when two parents are avalible, and all sites are used if one parent is absence.

For delta snp-index, the direction is meaningful for bi-parents populations (f2, ril/dh, bc) only at least one parent is used. And the direction is not meaningful for outcross population (f1/cp) in all case.

It is wise to set -ab T when the direction of delta snp-index is not meaningful.

For mutmap, please give the same sample name for -b1 and -p1(wild parent), the index of bulk1(wild parent actually) will not be calculated.

Reference:

SNP-Index and Delta SNP-index, [Fekih R, Takagi H, Tamiru M, et al. MutMap+: genetic mapping and mutant identification without crossing in rice[J]. PloS one, 2013, 8(7): e68529.](https://doi.org/10.1371/journal.pone.0068529)

G-statistic, [Magwene P M, Willis J H, Kelly J K. The statistics of bulk segregant analysis using next generation sequencing[J]. PLoS Comput Biol, 2011, 7(11): e1002255.](https://doi.org/10.1371/journal.pcbi.1002255)

Euclidean distance, [Hill J T, Demarest B L, Bisgrove B W, et al. MMAPPR: mutation mapping analysis pipeline for pooled RNA-seq[J]. Genome research, 2013, 23(4): 687-697.](https://genome.cshlp.org/content/23/4/687.short)

## step 1 simulation
### get help infomation
perl simulation_v2.pl
```
Function: simulation and get delta snpindex confidence intervals
	for given population type(f2, ril, bc1), depth and bulk size.

	-k   <str>    output prefix                  [force]
	-od  <dir>    output directory               [force]
	-pt  <str>    pop type (f2, ril, bc)         [force]
	-sd  <dir>    shell dir                      [od]
	-s1  <str>    bulk1(dominant/wild) size      [30]
	-s2  <str>    bulk2(recessive/mutant) size   [30]
	-rp  <int>    replications number            [10000]
	-ci  <int>    confidence intervals           [95,99]
	-md  <int>    min depth                      [10]
	-xd  <int>    max depth                      [500]
	-mi  <float>  min snp index                  [0]
	-rb  <bin>    Rscript bin                    [/Bio/bin/Rscript-3.6.0]

	-h            help document
```
### command line
```
perl simulation_v2.pl -k bsa -od bsa -pt ril -s1 30 -s2 30 -rp 10000 -ci 95,99 -md 10 -xd 500 -mi 0 -rb /path/to/Rscript
```

### results
\*.cisim.xls

## step 2 bsa calculation
### get help infomation
perl bsaindex.pl
```
Function: calculate (snpindex, g-statistic, ed, fisher exact test(fet))
	  of two bulks using vcf file with/without parent(s)

	-v   <file>    vcf file                       [force]
	-k   <str>     output prefix                  [force]
	-od  <dir>     output directory               [force]
	-b1  <str>     bulk1(dominant/wild)           [force]
	-b2  <str>     bulk2(recessive/mutant)        [force]
	-p1  <str>     dominant/wild parent           [optional]
	-p2  <str>     recessive/mutant parent        [optional]
	-pt  <str>     pop type (f1/cp, f2, ril, bc)  [ril]
	-md  <int>     min depth                      [10]
	-xd  <int>     max depth                      [500]
	-pd  <int>     parent min depth               [5]
	-mt  <str>     method(index,gst,ed,fet,all)   [all]
	-ep  <int>     ed power                       [5]
	-mi  <float>   min snp index                  [0.2]
	-sf  <file>    simulation result file         [optional]
	-sd  <str/int> min,max,other number in sf     [max]
	-id  <T/F>     use(T) indel or not(F)         [F]
	-ab  <T/F>     output absolutely delta(T)     [F]
	-rb  <bin>     Rscript bin path               [Rscript]
	-h             help document
```
### command line
```
perl bsaindex.pl -v snp.indel.vcf -k snp -od bsa -b1 bulk1 -b2 bulk2 -p1 parent1 -p2 parent2 -pt ril -ep 2 -mt all -sf bsa/snp.cisim.xls
```
if you do not care about direction of delta snpindex, -sf can be ignored and set -ab T 

### results
├── bsa.all.xls ## all index file

├── bsa.basic.xls ## used sites informations

├── bsa.drop.vcf ## dropped sites

├── bsa.ed.xls ## ed results

├── bsa.fet.xls ## fisher exact test results

├── bsa.gst.xls ## g-statistic results

└── bsa.index.xls ## snp-index results

# sliding windows and ploting

## step 1 sliding windows

### get help infomation
perl slidewindow.pl
```
Function: sliding window calculate smooth(mean) of snp index

        -i   <file>   snpindex file                    [force]
        -k   <str>    output prefix                    [force]
        -o   <dir>    output directory                 [force]
        -f   <file>   chromosome fai file              [force]
        -w   <int>    window length(kb)                [1000]
        -s   <int>    step length(kb)                  [100]
        -cp  <str>    chr pos columns                  [1,2]
        -cv  <str>    values ( index && ci ) columns   [3]
        -ms  <int>    min snp number in a window       [10]
        -lg  <T/F>    get -log10 value of mean value   [F]

        -h            help document
```
### command line
```
perl slidewindow.pl -i bsa/bsa.index.xls -k bsa.index -o bsa/ -f har.fa.fai -w 2000 -s 10 -cp 1,2 -cv 3-13 -ms 10
perl slidewindow.pl -i bsa/bsa.ed.xls -k bsa.ed -o bsa/ -f har.fa.fai -w 2000 -s 10 -cp 1,2 -cv 3 -ms 10
perl slidewindow.pl -i bsa/bsa.gst.xls -k bsa.gst -o bsa/ -f har.fa.fai -w 2000 -s 10 -cp 1,2 -cv 3 -ms 10
perl slidewindow.pl -i bsa/bsa.fet.xls -k bsa.fet -o bsa/ -f har.fa.fai -w 2000 -s 10 -cp 1,2 -cv 3,5 -ms 10 -lg T
```
fai file can be obtained from samtools fasta index or just add chromosome length manully:
```
chr01	18757413
chr02	15481669
chr03	14254934
chr04	13967672
chr05	13931063
chr06	13863567
chr07	13341776
chr08	13100779
chr09	13064641
chr10	12996920
```

### results
├── bsa.ed.mean.xls ## ed sliding window results

├── bsa.fet.mean.xls ## fisher exact test sliding window results

├── bsa.gst.mean.xls ## g-statistic sliding window results

└── bsa.index.mean.xls ## snp-index sliding window results

## step 2 plotting
bsaplot.r, plot.scanone.r (copy from [R/qtl](https://github.com/kbroman/qtl)), and point_line_plot.pl must be in the same directory.

### get help infomation
perl point_line_plot.pl

Don't be afraid of seeing so many options, just some of them are required!

```
Function: point and/or line plot.

	-k   <str>    output prefix                       [force]
	-f   <file>   genome fai file                     [force]
	-od  <dir>    output directory                    [force]
	-sd  <dir>    shell directory                     [od]
	-ch  <str>    chr name(s) use for plotting        [optional]
	-cm  <int>    chr start position                  [1]
	-lf  <file>   line plot file                      [NULL]
	-lp  <ints>   line file chr,pos columns           [1,2]
	-lv  <ints>   line file values columns            [3,4,5,6,7]
	-lc  <color>  line plot colors (len(lc)==len(lv)) [black,blue,red,blue,red]
	-pf  <file>   point plot file                     [NULL]
	-pp  <ints>   point plot chr,pos columns          [1,2]
	-pv  <ints>   point plot values columns           [3]
	-pc  <colors> point colors                        [#33A02C,#FF7F00]
	-hl  <float>  horizon line position               [0]
	-hc  <color>  horizon line color                  [grey]
	-ht  <int>    horizon line type(R)                [1]
	-cg  <int>    gap between chromosomes             [0]
	-ym  <float>  ymin                                [1*-1]
	-yx  <float>  ymax                                [1]
	-xl  <str>    xlab                                [Chromosome]
	-yl  <str>    ylab                                [expression(Delta(SNPindex))]
	-ca  <float>  axis cex                            [1.5]
	-cl  <float>  labels cex                          [2]
	-cp  <float>  point cex                           [0.3]
	-pch <int>    point type(R)                       [19]
	-lwd <float>  line width                          [2.5]
	-las <1/2>    parallele(1)/vertical(2) x-axis lab [1]
	-bdc <color>  band color(for distinguish chrs)    [grey90]
	-alt <T/F>    alt chr labels(avoid overlaping)    [F]
	-width  <int> plot width in inches                [15]
	-height <int> plot height in inches               [6]
	-res    <int> png resolution ratio                [300]
	-rb  <bin>    Rscript bin                         [/Bio/bin/Rscript-3.6.0]

	-h            help document
```
### command line
```
perl point_line_plot.pl -k bsa.index -f har.fa.fai -od bsa/ -lf bsa/bsa.index.mean.xls -lp 1,4 -lv 7,8,9,12,13 -lc "black,#2121D9,#2121D9,#DF0101,#DF0101" -pf bsa/bsa.index.xls -pp 1,2 -pv 5 -ym 1*-1 -yx 1 -xl "" -yl "expression(Delta(SNPindex))" -las 2
perl point_line_plot.pl -k bulk1.index -f har.fa.fai -od bsa/ -lf bsa/bsa.index.mean.xls -lp 1,4 -lv 5,10,14 -lc "black,#2121D9,#DF0101" -pf bsa/bsa.index.xls -pp 1,2 -pv 3 -ym NULL -yx 1 -xl "" -yl "SNPindex(bulk1)" -las 2
perl point_line_plot.pl -k bulk2.index -f har.fa.fai -od bsa/ -lf bsa/bsa.index.mean.xls -lp 1,4 -lv 6,11,15 -lc "black,#2121D9,#DF0101" -pf bsa/bsa.index.xls -pp 1,2 -pv 4 -ym NULL -yx 1 -xl "" -yl "SNPindex(bulk2)" -las 2
perl point_line_plot.pl -k bsa.ed -f har.fa.fai -od bsa/ -lf bsa/bsa.ed.mean.xls -lp 1,4 -lv 5 -lc "black" -pf bsa/bsa.ed.xls -pp 1,2 -pv 3 -hl 0.149,0.753 -hc "#2121D9,#DF0101" -ht 1,1 -ym 0 -yx NULL -xl "" -yl "ED2" -las 2
perl point_line_plot.pl -k bsa.gst -f har.fa.fai -od bsa/ -lf bsa/bsa.gst.mean.xls -lp 1,4 -lv 5 -lc "black" -pf bsa/bsa.gst.xls -pp 1,2 -pv 3 -hl 4.786,6.639 -hc "#2121D9,#DF0101" -ht 1,1 -ym 0 -yx NULL -xl "" -yl "G" -las 2
perl point_line_plot.pl -k bsa.fet -f har.fa.fai -od bsa -lf bsa/bsa.fet.mean.xls -lp 1,4 -lv 6 -lc "black" -pf bsa/bsa.fet.xls -pp 1,2 -pv 4 -hl 1.976,2.833 -hc "#2121D9,#DF0101" -ht 1,1 -ym 0 -yx NULL -xl "" -yl "expression(-log[10](italic(p)))" -las 2
```

for -hl option, you need to calculate these values using quantile/percentile or other method first, or you just run qtl_region.pl first, and the corresponding value will be printed in the standard output.

fai file here ony contain chr/scaffolds use for plotting.

#  get qtl region
### get help infomation
perl qtl_region.pl

```
Function: get significant region

	-i      <str>    input file                            [force]
	-k      <str>    output prefix                         [force]
	-o      <dir>    output directory                      [force]
	-chr    <int>    chr column                            [force]
	-value  <int>    value column                          [force]
	-start  <int>    window start column                   [undef]
	-end    <int>    window end column                     [undef]
	-lower  <int>    lower bound column                    [undef]
	-upper  <int>    upper bound column                    [undef]
	-alpha  <fload>  alpha for p-value bf correction       [undef]
	-uthres <fload>  upper bound hard threshold            [undef]
	-lthres <fload>  lower bound hard threshold            [undef]
	-quant  <int>    percentile to defined threshold       [undef]
	-type   <str>    significant type smaller/larger/both  [smaller]
	-h              help document

priority:
lower/upper > thres > alpha > quant
```
### command line
```
perl qtl_region.pl -i bsa/bsa.index.mean.xls -k bsa.index95 -o bsa/ -chr 1 -value 7 -start 2 -end 3 -lower 8 -upper 9
perl qtl_region.pl -i bsa/bsa.index.mean.xls -k bsa.index99 -o bsa/ -chr 1 -value 7 -start 2 -end 3 -lower 12 -upper 13
perl qtl_region.pl -i bsa/bsa.ed.mean.xls -k bsa.ed95 -o bsa/ -chr 1 -start 2 -end 3 -value 5 -quant 95 -type larger
perl qtl_region.pl -i bsa/bsa.ed.mean.xls -k bsa.ed99 -o bsa/ -chr 1 -start 2 -end 3 -value 5 -quant 99 -type larger
perl qtl_region.pl -i bsa/bsa.gst.mean.xls -k bsa.gst95 -o bsa/ -chr 1 -start 2 -end 3 -value 5 -quant 95 -type larger
perl qtl_region.pl -i bsa/bsa.gst.mean.xls -k bsa.gst99 -o bsa/ -chr 1 -start 2 -end 3 -value 5 -quant 99 -type larger
perl qtl_region.pl -i bsa/bsa.fet.mean.xls -k bsa.fet0.05 -o bsa/-chr 1 -start 2 -end 3 -value 7 -lthres 0.05 -type smaller
perl qtl_region.pl -i bsa/bsa.fet.mean.xls -k bsa.fet0.01 -o bsa/-chr 1 -start 2 -end 3 -value 7 -lthres 0.01 -type smaller
```
### results
*.sig.xls ## possible qtl region of each index
