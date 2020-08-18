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
my ($key, $outdir, $shdir);
my ($repeat, $ci_interval);
my ($b1size, $b2size, $popt);
my ($min_depth, $max_depth, $min_index);
my $rbin;
GetOptions(
	"h|?"=>\&help,
	"k:s"=>\$key,
	"od:s"=>\$outdir,
	"pt:s"=>\$popt,
	"sd:s"=>\$shdir,
	"s1:s"=>\$b1size,
	"s2:s"=>\$b2size,
	"rp:s"=>\$repeat,
	"ci:s"=>\$ci_interval,
	"md:s"=>\$min_depth,
	"xd:s"=>\$max_depth,
	"mi:s"=>\$min_index,
	"rb:s"=>\$rbin,
) || &help;
&help unless ($key && $popt && $outdir);

sub help
{
	print "Function: simulation and get delta snpindex confidence intervals
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
";
	exit;
}
############# start time #################
my $current_T = &date_format(localtime());
print "Programe start: $current_T\n\n";
my $begin_time = time();
##########################################
$popt ||= "ril";
$popt = lc($popt);
$popt =~ s/\s+//g;
$min_depth ||= 10;
$max_depth ||= 500;
$min_index //= 0;
$b1size ||= 30;
$b2size ||= 30;
$repeat ||= 10000;
$ci_interval ||= "95,99";
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
# r command
my $rscript = "$shdir/$key\.cisim.r";
my $outfile = "$outdir/$key\.cisim.xls";
open(R, ">$rscript");
print R "setwd(\"$outdir\")
bulkSize <- c\($b1size,$b2size\)
replications <- $repeat
fvalue <- $min_index
popStruc <- \"$popt\"
depth <- $min_depth:$max_depth
interv <- c\($ci_interval\)
outfile <- \"$outfile\"
source(\"$Bin/simulation_v2.r\")
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

