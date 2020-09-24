# bsa
Bulked-Segregant Analysis using vcf file with or without parents

dependences, you need to install these dependence module from cpan first

use Cwd

use Getopt::Long

use Data::Dumper

use File::Basename

use Text::NSP::Measures::2D::Fisher::twotailed

use FindBin

# step 1 do simulation
simulation_v2.pl and simulation_v2.r must be in the same directory

if you do not care about the direction of delta snpindex, you can skip this step
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
perl simulation_v2.pl -k snp -od bsa -pt ril -s1 50 -s2 50 -md 10 -xd 500 -rb /path/to/your/Rscript
```
### results
snp.cisim.xls is the main result file

# step 2 index calculation
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
	-p1  <str>     parent of bulk1                [optional]
	-p2  <str>     parent of bulk2                [optional]
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

# step 3 slidding windows
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
perl slidewindow.pl bsa/bsa.index.xls -k bsa.index -o bsa/ -f har.fa.fai -w 2000 -s 10 -cp 1,2 -cv 3-13 -ms 10
perl slidewindow.pl bsa/bsa.ed.xls -k bsa.ed -o bsa/ -f har.fa.fai -w 2000 -s 10 -cp 1,2 -cv 3 -ms 10
perl slidewindow.pl bsa/bsa.gst.xls -k bsa.gst -o bsa/ -f har.fa.fai -w 2000 -s 10 -cp 1,2 -cv 3 -ms 10
perl slidewindow.pl bsa/bsa.fet.xls -k bsa.fet -o bsa/ -f har.fa.fai -w 2000 -s 10 -cp 1,2 -cv 3,5 -ms 10 -lg T
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

after this step, you can use cmplot(https://github.com/YinLiLin/CMplot) to plot window results

# step 4 get qtl region
### get help infomation
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
### commad line
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
