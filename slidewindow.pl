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
my ($snpindex, $key, $fai, $outdir);
my ($window, $step, $nlog10);
my ($cpos, $cvalue, $min_snp);
GetOptions(
	"h|?"=>\&help,
	"i:s"=>\$snpindex,
	"k:s"=>\$key,
	"o:s"=>\$outdir,
	"f:s"=>\$fai,
	"w:s"=>\$window,
	"s:s"=>\$step,
	"cp:s"=>\$cpos,
	"cv:s"=>\$cvalue,
	"ms:s"=>\$min_snp,
	"lg:s"=>\$nlog10,
) || &help;
&help unless ($snpindex && $key && $fai && $outdir);

sub help
{
	print "Function: sliding window calculate smooth(mean) of snp index

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
";
	exit;
}
############# start time #################
my $current_T = &date_format(localtime());
print "Programe start: $current_T\n\n";
my $begin_time = time();
##########################################
$snpindex = &abs_dir($snpindex);
$fai = &abs_dir($fai);
$window ||= 1000;
$step ||= 100;
$window = $window * 1000;
$step = $step * 1000;
$min_snp ||= 10;
$cpos ||= "1,2";
$cvalue ||= "12,13,14,15,16,17,18";
$cpos =~ s/\s+//g;
$cvalue =~ s/\s+//g;
$nlog10 ||= "F";
## get chr pos and values column index
my ($cindex, $pindex) = split(/,/, $cpos);
$cindex = $cindex - 1;
$pindex = $pindex - 1;
my @vtmp = split(/,/, $cvalue);
my @vcolumns;
foreach my $v (@vtmp) {
	if($v =~ /-/){
		my ($vstart, $vend) = split(/-/, $v);
		for (my $i = $vstart; $i <= $vend; $i++) {
			push(@vcolumns, $i);
		}
	}else{
		push(@vcolumns, $v);
	}
}
my @vindex;
for (my $i = 0; $i < @vcolumns; $i++) {
	push(@vindex, $vcolumns[$i] - 1)
}

## get chr length from fai file
my %chromh;
open(IN, $fai);
while (<IN>) {
	chomp;
	next if(/^$/);
	my @line = split(/\s+/, $_);
	$chromh{$line[0]} = $line[1];
}
close(IN);
my @chroma = &nsort(keys %chromh);

## read datas;
my %data;
my $vhead = "chr\tstart\tend\tmid";
my $row = 0;
open(IN, $snpindex);
while (<IN>) {
	chomp;
	next if(/^$/);
	my @line = split(/\t/, $_);
	if($row == 0){
		for (my $i = 0; $i < @vindex; $i++) {
			$vhead .= "\t$line[$vindex[$i]]";
			if($nlog10 eq "T"){
				$vhead .= "\t$line[$vindex[$i]](Nlog10)";
			}
		}
		$row++;
		next;
	}
	my $chr = $line[$cindex];
	if(!exists $chromh{$chr}){
		$row++;
		next;
	}
	my $pos = $line[$pindex];
	for (my $i = 0; $i < @vindex; $i++) {
		push(@{$data{$chr}{$pos}}, $line[$vindex[$i]]);
	}
	$row++;
}
close(IN);

print "There are ";
print (scalar keys %data);
print " chromosomes(scaffolds) need to be slided.\n";

open(OUT, ">$outdir/$key\.mean.xls");
print OUT "$vhead\n";
#my %output;
my %win_count;
foreach my $chr (@chroma) {
	&show_log("sliding chr $chr ...");
	$win_count{$chr}{e} = 0;
	$win_count{$chr}{t} = 0;
	my @pos = keys %{$data{$chr}};
	@pos = sort{$a <=> $b} @pos;
	my $posi = 0;
	my @ieffect;
	my $last_start = 0;
	my $nsnp = 0;
	my @vsum;
	## i-th window
	for (my $i = 0; ($step * $i + $window) < $chromh{$chr}; $i++) { # i = i-th window 
		foreach my $j (@pos){
			if($j >= $step * $i && $j < $step * $i + $window){
				$nsnp++;
				my @line = @{$data{$chr}{$j}};
				if(@line > @vindex){
					die("ERROR: $chr $j value count larger than value index number!\n");
				}
				for (my $v = 0; $v < @line; $v++) {
					if(defined $vsum[$v]){
						$vsum[$v] += $line[$v];
					}else{
						$vsum[$v] = $line[$v];
					}
				}
			}
			if($j < $step * ($i + 1)){
				shift @pos;
			}
		}
		push (@ieffect, $i);
		$win_count{$chr}{t}++;
		my @vmean;
		if($nsnp >= $min_snp){
			$win_count{$chr}{e}++;
			for (my $v = 0; $v < @vsum; $v++) {
				my $avg = $vsum[$v]/$nsnp;
				push(@vmean, $avg);
				if($nlog10 eq "T"){
					my $vlog10;
					if($avg <= 0){
						$vlog10 = "NA";
					}else{
						$vlog10 = &log10($avg) * -1;
					}
					push(@vmean, $vlog10);
				}
			}
			my $wstart = $step * $ieffect[-1] + 1;
			my $wend = $wstart + $window - 1;
			my $wmid = int(($wstart + $wend)/2);
			print OUT "$chr\t$wstart\t$wend\t$wmid\t", join("\t", @vmean), "\n";
			$nsnp = 0;
			@vsum = ();
		}else{
			# merge into next window
		}
		$last_start = $step * $ieffect[-1] + 1;
	}

	if($last_start < $chromh{$chr} && @pos > 0){ ## for last window
		if($nsnp >= $min_snp){
			$nsnp = 0;
			@vsum = ();
		}
		$win_count{$chr}{t}++;
		foreach my $j (@pos){
			$nsnp++;
			my @line = @{$data{$chr}{$j}};
			if(@line > @vindex){
				die("ERROR: $chr $j value count larger than value index number!\n");
			}
			for (my $v = 0; $v < @line; $v++) {
				if(defined $vsum[$v]){
					$vsum[$v] += $line[$v];
				}else{
					$vsum[$v] = $line[$v];
				}
			}
		}
		my @vmean;
		if($nsnp >= $min_snp){
			$win_count{$chr}{e}++;
			for (my $v = 0; $v < @vsum; $v++) {
				my $avg = $vsum[$v]/$nsnp;
				push(@vmean, $avg);
				if($nlog10 eq "T"){
					my $vlog10;
					if($avg <= 0){
						$vlog10 = "NA";
					}else{
						$vlog10 = &log10($avg) * -1;
					}
					push(@vmean, $vlog10);
				}
			}
			my $wstart;
			if(defined $ieffect[-1]){
				$wstart = $step * ($ieffect[-1] + 1) + 1;
			}else{
				$wstart = 1;
			}
			my $wend = $chromh{$chr};
			my $wmid = int(($wstart + $wend)/2);
			my $jmean = join("\t", @vmean);
			print OUT "$chr\t$wstart\t$wend\t$wmid\t", join("\t", @vmean), "\n";
		}else{
			# ignore
		}
	}
	@vsum = ();
	$nsnp = 0;
	&show_log("done.");
}
close(OUT);
open(STAT, ">$outdir/$key\.stat.xls");
print STAT "chr\ttotal\tavailable\tpercent(\%)\n";
my $total_all = 0;
my $total_eff = 0;
foreach my $chr (@chroma) {
	my $percent = 0;
	if($win_count{$chr}{t} > 0){
		$percent = sprintf("%.4f", $win_count{$chr}{e}/$win_count{$chr}{t});
	}
	$percent = $percent * 100;
	print STAT "$chr\t$win_count{$chr}{t}\t$win_count{$chr}{e}\t$percent\n";
	$total_all += $win_count{$chr}{t};
	$total_eff += $win_count{$chr}{e};
}
my $tpercent = sprintf("%.4f", $total_eff/$total_all);
$tpercent = $tpercent * 100;
print STAT "Total\t$total_all\t$total_eff\t$tpercent\n";
close(STAT);

############# end time ###################
$current_T = &date_format(localtime());
print "Programe end: $current_T\n\n";
&runtime($begin_time);
##########################################
sub log10()
{
	my $n = shift;
	return log($n)/log(10);
}
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
