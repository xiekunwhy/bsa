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
my ($infile, $key, $outdir);
my ($cchr, $cstart, $cend);
my ($cvalue, $clower, $cupper);
my ($uthres, $lthres, $alpha, $quantile);
my $type;
GetOptions(
	"h|?"     =>\&help,
	"i:s"     =>\$infile,
	"k:s"     =>\$key,
	"o:s"     =>\$outdir,
	"chr:s"   =>\$cchr,
	"value:s" =>\$cvalue,
	"start:s" =>\$cstart,
	"end:s"   =>\$cend,
	"lower:s" =>\$clower,
	"upper:s" =>\$cupper,
	"alpha:s" =>\$alpha,
	"uthres:s"=>\$uthres,
	"lthres:s"=>\$lthres,
	"quant:s" =>\$quantile,
	"type:s"  =>\$type
) || &help;
&help unless ($infile && $key && $outdir && $cchr && $cvalue);

sub help
{
	print "Function: get significant region

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

";
	exit;
}
############# start time #################
my $current_T = &date_format(localtime());
print $0;
print "\nPrograme start: $current_T\n\n";
my $begin_time = time();
##########################################
$infile = &abs_dir($infile);
if(!-d $outdir){
	mkdir $outdir;
}
$outdir = &abs_dir($outdir);

## position columns
$cchr = $cchr - 1;
$cvalue = $cvalue - 1;

if(!defined $cstart && !defined $cend){
	die("ERROR: one of start and end must specified!\n");
}elsif(defined $cstart && !defined $cend){
	$cstart = $cstart - 1;
	$cend = $cstart;
}elsif(!defined $cstart && defined $cend){
	$cend = $cend - 1;
	$cstart = $cend;
}else{
	$cstart = $cstart - 1;
	$cend = $cend - 1;
}
$type ||= "smaller";
## threshold
# confidence
if(defined $clower || defined $cupper){
	undef $alpha;
	undef $uthres;
	undef $lthres;
	undef $quantile;
	$clower = $clower - 1 if defined $clower;
	$cupper = $cupper - 1 if defined $cupper;
}

# hard threshold
if(defined $uthres || defined $lthres){
	undef $alpha;
	undef $quantile;
}

# alpha or quantile
if(defined $alpha || defined $quantile){
	my @all_values;
	my $snpn = 0;
	open(IN, $infile);
	<IN>;
	while (<IN>) {
		chomp;
		next if(/^$/);
		my @line = split(/\t/, $_);
		push(@all_values, $line[$cvalue]);
		$snpn++;
	}
	close(IN);
	if(defined $alpha){
		$lthres = $alpha/$snpn if defined $alpha;
		undef $quantile;
	}
	if(defined $quantile){
		@all_values = sort{$a <=> $b} @all_values;
		my $percentile_pos1 = int($quantile/100 * $snpn);
		if($type eq "both"){
			$uthres = $all_values[$percentile_pos1];
			my $percentile_pos2 = int((1 - $quantile/100) * $snpn);
			$lthres = $all_values[$percentile_pos2];
		}elsif($type eq "larger"){
			$uthres = $all_values[$percentile_pos1];
		}elsif($type eq "smaller"){
			$lthres = $all_values[$percentile_pos1];
		}
	}
	undef @all_values;
}
if(defined $uthres){
	if(abs($uthres) < 1){
		$uthres = sprintf("%.3g", $uthres);
	}else{
		$uthres = sprintf("%.3f", $uthres);
	}
	print "uthres $uthres\n";
}

if(defined $lthres){
	if(abs($lthres) < 1){
		$lthres = sprintf("%.3g", $lthres);
	}else{
		$lthres = sprintf("%.3f", $lthres);
	}
	print "lthres $lthres\n";
}

if(defined $lthres && defined $uthres){
	my $tmp = $lthres;
	if($lthres > $uthres){
		$lthres = $uthres;
		$uthres = $tmp;
	}
}
## state defined
my %state;
open(IN, $infile);
my $head = <IN>;
chomp($head);
while (<IN>) {
	chomp;
	my @line = split(/\t/, $_);
	my $chr = $line[$cchr];
	my $pos = $line[$cstart];
	my $value = $line[$cvalue];
	$state{$chr}{$pos}{"end"} = $line[$cend];
	$state{$chr}{$pos}{"val"} = $value;
	if(defined $lthres && defined $uthres){
		if($value <= $lthres){
			$state{$chr}{$pos}{"sig"} = "l";
		}elsif($value >= $uthres){
			$state{$chr}{$pos}{"sig"} = "u";
		}else{
			$state{$chr}{$pos}{"sig"} = "n";
		}
	}elsif(defined $lthres){
		if($value <= $lthres){
			$state{$chr}{$pos}{"sig"} = "l";
		}else{
			$state{$chr}{$pos}{"sig"} = "n";
		}
	}elsif(defined $uthres){
		if($value >= $uthres){
			$state{$chr}{$pos}{"sig"} = "u";
		}else{
			$state{$chr}{$pos}{"sig"} = "n";
		}
	}else{
		if(defined $cupper && defined $clower){
			if($value >= $line[$cupper]){
				$state{$chr}{$pos}{"sig"} = "u";
			}elsif($value <= $line[$clower]){
				$state{$chr}{$pos}{"sig"} = "l";
			}else{
				$state{$chr}{$pos}{"sig"} = "n";
			}
		}elsif(defined $cupper){
			if($value >= $line[$cupper]){
				$state{$chr}{$pos}{"sig"} = "u";
			}else{
				$state{$chr}{$pos}{"sig"} = "n";
			}
		}elsif(defined $clower){
			if($value <= $line[$clower]){
				$state{$chr}{$pos}{"sig"} = "l";
			}else{
				$state{$chr}{$pos}{"sig"} = "n";
			}
		}
	}
}
close(IN);

## get chromosomes
my @chrom = &nsort(keys %state);

## continuity region defined
my %region;
foreach my $chr (@chrom) {
	my @tmp_pos;
	my @tmp_val;
	my $state_former;
	foreach my $start (sort{$a <=> $b} keys %{$state{$chr}}) {
		my $end = $state{$chr}{$start}{end};
		my $tmp_state = $state{$chr}{$start}{sig};
		my $tmp_value = $state{$chr}{$start}{val};
		if(!defined $state_former){
			$state_former = $tmp_state;
			push(@tmp_pos, $start);
			if($start != $end){
				push(@tmp_pos, $end);
			}
			push(@tmp_val, $tmp_value);
		}else{
			if($tmp_state eq $state_former){
				push(@tmp_pos, $start);
				if($start != $end){
					push(@tmp_pos, $end);
				}
				push(@tmp_val, $tmp_value);
				next;
			}else{
				if($state_former ne "n"){
					my @pos_sort = sort{$a <=> $b} @tmp_pos;
					my @val_sort = sort{$a <=> $b} @tmp_val;
					$region{$chr}{$pos_sort[0]}{end} = $pos_sort[-1];
					$region{$chr}{$pos_sort[0]}{sig} = $state_former;
					if($state_former eq "u"){
						$region{$chr}{$pos_sort[0]}{val} = $val_sort[-1];
					}elsif($state_former eq "l"){
						$region{$chr}{$pos_sort[0]}{val} = $val_sort[0];
					}
				}
				@tmp_pos = ();
				@tmp_val = ();
				push(@tmp_pos, $start);
				push(@tmp_pos, $end);
				push(@tmp_val, $tmp_value);
				$state_former = $tmp_state;
			}
		}
	}
	## last possible qtl
	if($state_former ne "n"){
		my @pos_sort = sort{$a <=> $b} @tmp_pos;
		my @val_sort = sort{$a <=> $b} @tmp_val;
		$region{$chr}{$pos_sort[0]}{end} = $pos_sort[-1];
		$region{$chr}{$pos_sort[0]}{sig} = $state_former;
		if($state_former eq "u"){
			$region{$chr}{$pos_sort[0]}{val} = $val_sort[-1];
		}elsif($state_former eq "l"){
			$region{$chr}{$pos_sort[0]}{val} = $val_sort[0];
		}
	}
}

## merge overlap region
my %merge;
foreach my $chr (@chrom) {
	if(exists $region{$chr}){
		my $start_former;
		my $end_former;
		my $val_former;
		my $state_former;
		foreach my $start (sort {$a <=> $b} keys %{$region{$chr}}) {
			my $end_tmp = $region{$chr}{$start}{end};
			my $val_tmp = $region{$chr}{$start}{val};
			my $state_tmp = $region{$chr}{$start}{sig};
			if(!defined $start_former){
				$start_former = $start;
				$end_former = $end_tmp;
				$val_former = $val_tmp;
				$state_former = $state_tmp;
			}else{
				if(($start <= $end_former) && ($state_tmp eq $state_former)){
					$end_former = $end_tmp;
					if(($state_tmp eq "u") && ($val_tmp > $val_former)){
						$val_former = $val_tmp;
					}elsif(($state_tmp eq "l") && ($val_tmp < $val_former)){
						$val_former = $val_tmp;
					}
				}else{
					$merge{$chr}{$start_former}{end} = $end_former;
					$merge{$chr}{$start_former}{val} = $val_former;
					$merge{$chr}{$start_former}{sig} = $state_former;
					$start_former = $start;
					$end_former = $end_tmp;
					$val_former = $val_tmp;
				}
				$state_former = $state_tmp;
			}
		}
		## last region
		if(defined $start_former){
			$merge{$chr}{$start_former}{end} = $end_former;
			$merge{$chr}{$start_former}{val} = $val_former;
			$merge{$chr}{$start_former}{sig} = $state_former;
		}
	}
}

## output
open(OUT, ">$outdir/$key\.sig.xls");
print OUT "chr\tstart\tend\tdirection\tpeak\n";
foreach my $chr (@chrom) {
	if(exists $merge{$chr}){
		foreach my $start (sort {$a <=> $b} keys %{$merge{$chr}}) {
			print OUT "$chr\t$start\t$merge{$chr}{$start}{end}";
			if($merge{$chr}{$start}{sig} eq "u"){
				print OUT "\tupper";
			}elsif($merge{$chr}{$start}{sig} eq "l"){
				print OUT "\tlower";
			}
			print OUT "\t$merge{$chr}{$start}{val}\n";
		}
	}
}
close(OUT);

############# end time ###################
$current_T = &date_format(localtime());
print $0;
print "\nPrograme end: $current_T\n\n";
&runtime($begin_time);
##########################################
# natural sort
sub nsort {
	my($cmp, $lc);
	return @_ if @_ < 2;   # Just to be CLEVER.
	my($x, $i);  # scratch vars
	map
		$_->[0],
	sort {
		# Uses $i as the index variable, $x as the result.
		$x = 0;
		$i = 1;
		while($i < @$a and $i < @$b) {
			last if ($x = ($a->[$i] cmp $b->[$i])); # lexicographic
			++$i;
			last if ($x = ($a->[$i] <=> $b->[$i])); # numeric
			++$i;
		}
		$x || (@$a <=> @$b) || ($a->[0] cmp $b->[0]);
	}
	map {
		my @bit = ($x = defined($_) ? $_ : '');
		if($x =~ m/^[+-]?(?=\d|\.\d)\d*(?:\.\d*)?(?:[Ee](?:[+-]?\d+))?\z/s) {
			# It's entirely purely numeric, so treat it specially:
			push @bit, '', $x;
		} else {
			# Consume the string.
			while(length $x) {
				push @bit, ($x =~ s/^(\D+)//s) ? lc($1) : '';
				push @bit, ($x =~ s/^(\d+)//s) ?    $1  :  0;
			}
		}
		\@bit;
	}
	@_;
}

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
