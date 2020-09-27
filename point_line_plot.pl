#!/usr/bin/perl -w
use strict;
use Cwd; ## get current directory
use Getopt::Long qw(:config no_ignore_case);
use Data::Dumper;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);
my $version = "v1.0";
my $writer = "xiek";
my $programe_date = "2020-2-19";
###### get options and help documents ######
my ($key, $outdir, $fai, $chrmin, $chr, $shdir);
my ($linef, $linep, $linev, $linec);
my ($pointf, $pointp, $pointv, $pointc);
my ($hline, $hcolor, $hlty);
my ($ymin, $ymax);
my ($xlab, $ylab);
my ($cexa, $cexl, $cexp);
my ($pch, $lwd, $las, $bandc, $alt);
my ($width, $height, $res);
my $cgap;
my $rbin;
GetOptions(
	"h|?"=>\&help,"k:s"=>\$key,"f:s"=>\$fai,"od:s"=>\$outdir,"ch:s"=>\$chr, "sd:s"=>\$shdir,
	"cm:s"=>\$chrmin,"lf:s"=>\$linef,"lp:s"=>\$linep,"lv:s"=>\$linev,"lc:s"=>\$linec,
	"pf:s"=>\$pointf,"pp:s"=>\$pointp,"pv:s"=>\$pointv,"pc:s"=>\$pointc,
	"hl:s"=>\$hline,"hc:s"=>\$hcolor,"ht:s"=>\$hlty,"cg:s"=>\$cgap,
	"ym:s"=>\$ymin,"yx:s"=>\$ymax,"xl:s"=>\$xlab,"yl:s"=>\$ylab,
	"ca:s"=>\$cexa,"cl:s"=>\$cexl,"cp:s"=>\$cexp,
	"pch:s"=>\$pch,"lwd:s"=>\$lwd,"las:s"=>\$las,"bdc:s"=>\$bandc,"alt:s"=>\$alt,
	"width:s"=>\$width,"height:s"=>\$height,"res:s"=>\$res,
	"rb:s"=>\$rbin,
) || &help;
&help unless ($key && $fai && $outdir);

sub help
{
	print "Function: point and/or line plot.

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
";
	exit;
}
############# start time #################
my $current_T = &date_format(localtime());
print "Programe start: $current_T\n\n";
my $begin_time = time();
##########################################
$fai = &abs_dir($fai);
die("ERROR: one of linef and pointf must be given!\n") if(!defined $linef && !defined $pointf);
if(defined $linef){
	$linef = &abs_dir($linef);
	$linef = "\"$linef\"";
}else{
	$linef = "NULL";
}
if(defined $pointf){
	$pointf = &abs_dir($pointf);
	$pointf = "\"$pointf\"";
}else{
	$pointf = "NULL";
}

$chrmin ||= 1;
$linep ||= "1,2";
$linev ||= "3,4";
$linec ||= "black,blue,red,blue,red";
$linec = join("\",\"", split(/,/, $linec));
$pointp ||= "1,2";
$pointv ||= 3;
$pointc ||= "#33A02C,#FF7F00";
$pointc = join("\",\"", split(/,/, $pointc));
$hline ||= "NULL";
$hcolor ||= "grey";
$hlty ||= 1;
if($hline =~ /,/){
	my @h_value = split(/\,/, $hline);
	my @h_color = split(/\,/, $hcolor);
	my @h_type = split(/\,/, $hlty);
	if(@h_value != @h_color){
		die("ERROR: length of -hl was not equal to -hc\n");
	}
	if(@h_value != @h_type){
		die("ERROR: length of -hl was not equal to -ht\n");
	}
	$hline = "c(" . join(",", @h_value) . ")";
	my @new_color;
	foreach my $c (@h_color) {
		push(@new_color, "\"$c\"");
	}
	$hcolor = "c(" . join(",", @new_color) . ")";
	$hlty = "c(" . join(",", @h_type) . ")";
}else{
	$hcolor = "\"$hcolor\"";
}
$cgap //= 0;
$ymin //= -1;
$ymax //= 1;
$xlab //= "Chromosome";
$ylab //= "expression(Delta(SNPindex))";
if($ylab !~ /expression/){
	$ylab = "\"$ylab\"";
}
$cexa ||= 1.5;
$cexl ||= 2;
$cexp ||= 0.3;
$pch ||= 19;
$lwd ||= 3;
$las ||= 1;
$bandc ||= "grey90";
$alt ||= "F";
$width ||= 15;
$height ||= 6;
$res ||= 300;
$rbin ||= "/Bio/bin/Rscript-3.6.0";
if(!-d $outdir){
	mkdir $outdir;
}
$outdir = &abs_dir($outdir);
$shdir ||= $outdir;
if(!-d $shdir){
	mkdir $shdir;
}
$shdir = &abs_dir($shdir);

my %chr_used;
if(defined $chr){
	$chr =~ s/\s+//g;
	my @chrs = split(/\,/, $chr);
	foreach my $c (@chrs){
		$chr_used{$c} = 1;
	}
}

# chrid and length
my $chrom;
my $chrmax;
open(IN, $fai);
while (<IN>) {
	chomp;
	next if(/^$/);
	my @line = split(/\s+/);
	next if(defined $chr && !exists $chr_used{$line[0]});
	$line[1] = $line[1] + $cgap;
	if(defined $chrom){
		$chrom .= "\",\"$line[0]";
		$chrmax .= ",$line[1]";
	}else{
		$chrom = $line[0];
		$chrmax = $line[1];
	}
}
close(IN);
# r command
my $rscript = "$shdir/$key\.plplot.r";
open(R, ">$rscript");
print R "setwd(\"$outdir\")
source(\"$Bin/plot.scanone.r\")
source(\"$Bin/bsaplot.r\")
bsaplot(
  key = \"$key\",
  chrom = c(\"$chrom\"),
  chrmax = c($chrmax),
  chrmin = $chrmin,
  linef = $linef,
  linep = c($linep),
  linev = c($linev),
  linec = c(\"$linec\"),
  pointf = $pointf,
  pointp = c($pointp),
  pointv = c($pointv),
  pointc = c(\"$pointc\"),
  hline = $hline,
  hcolor = $hcolor,
  hlty = $hlty,
  ymin = $ymin,
  ymax = $ymax,
  xlab = \"$xlab\",
  ylab = $ylab,
  cexa = $cexa,
  cexl = $cexl,
  cexp = $cexp,
  pch = $pch,
  lwd = $lwd,
  las = $las,
  bandc = \"$bandc\",
  alt = $alt,
  width = $width,
  height = $height,
  res = $res
)
";
close(R);

&run_or_die("$rbin $rscript");

############# end time ###################
$current_T = &date_format(localtime());
print "Programe end: $current_T\n\n";
&runtime($begin_time);
##########################################
# sub date format
sub date_format()
{
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year + 1900, $mon + 1, $day, $hour, $min, $sec);
}

# sub Runtime
sub runtime()
{
	my ($begin_time) = @_;
	my $now_time = time();
	my $total_time = $now_time - $begin_time;
	my $sec = 0; my $minu = 0; my $hour = 0;
	if($total_time >= 3600){
		$hour = int($total_time/3600);
		my $left = $total_time % 3600;
		if($left >= 60){
			$minu = int($left/60);
			$sec = $left % 60;
			$total_time = $hour."h\-".$minu."m\-".$sec."s";
		} else {
			$minu = 0;
			$sec = $left;
			$total_time = $hour."h\-".$minu."m\-".$sec."s";
		}
	} else {
		if($total_time >= 60){
			$minu = int($total_time/60);
			$sec = $total_time % 60;
			$total_time = $minu."m\-".$sec."s";
		} else {
			$sec = $total_time;
			$total_time = $sec."s";
		}
	}
	print "Total elapsed time [$total_time]\n\n";
}

# sub absolutely directory
sub abs_dir()
{
	my ($in) = @_;
	my $current_dir = `pwd`;
	chomp($current_dir);
	my $return_dir = "";
	if(-f $in){
		my $in_dir = dirname($in);
		my $in_file = basename($in);
		chdir $in_dir;
		$in_dir = `pwd`;
		chomp($in_dir);
		$return_dir = "$in_dir/$in_file";
	} elsif(-d $in) {
		chdir $in;
		my $in_dir = `pwd`;
		chomp($in_dir);
		$return_dir = $in_dir;
	} else {
		die("ERROR: there is no file or dir called [$in], please check!\n");
	}
	chdir $current_dir;
	return $return_dir;
}

# show log
sub show_log()
{
	my ($text) = @_;
	my $current_time = &date_format(localtime());
	print "$current_time: $text\n";
}

# run or die
sub run_or_die()
{
	my ($cmd) = @_;
	&show_log($cmd);
	my $flag = system($cmd);
	if($flag != 0){
		&show_log("ERROR: command $cmd");
		exit(1);
	}
	&show_log("done\n");
}

# qsub
sub qsub()
{
	my ($sh, $vf, $maxjob, $queue) = @_;
	$vf ||= "2G";
	$maxjob ||= 20;
	$queue ||= "all.q,avx.q";
	my $cmd_qsub = "/Bio/bin/qsub-sge.pl --convert no --queue $queue --interval 1 --maxjob $maxjob --resource vf=$vf $sh";
	&run_or_die($cmd_qsub);
}
