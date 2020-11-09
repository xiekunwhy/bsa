#!/usr/bin/perl -w
use strict;
use Cwd; ## get current directory
use Getopt::Long qw(:config no_ignore_case);
use Data::Dumper;
use File::Basename qw(basename dirname);
use Text::NSP::Measures::2D::Fisher::twotailed;
use FindBin qw($Bin $Script);
my $version = "v1.0";
my $writer = "xiek";
my $programe_date = "2020-2-19";
###### get options and help documents ######
my ($vcf, $parent1, $parent2, $bulk1, $bulk2, $popt);
my ($min_depth, $max_depth, $min_index, $pmin_depth, $indel);
my ($outdir, $key, $abs);
my ($method, $edpower);
my ($sim_file, $sim_depth);
my $rbin;
GetOptions(
	"h|?"=>\&help,
	"v:s"=>\$vcf,
	"k:s"=>\$key,
	"od:s"=>\$outdir,
	"b1:s"=>\$bulk1,
	"b2:s"=>\$bulk2,
	"p1:s"=>\$parent1,
	"p2:s"=>\$parent2,
	"pt:s"=>\$popt,
	"md:s"=>\$min_depth,
	"xd:s"=>\$max_depth,
	"pd:s"=>\$pmin_depth,
	"mt:s"=>\$method,
	"ep:s"=>\$edpower,
	"mi:s"=>\$min_index,
	"sf:s"=>\$sim_file,
	"sd:s"=>\$sim_depth,
	"id:s"=>\$indel,
	"ab:s"=>\$abs,
	"rb:s"=>\$rbin,
) || &help;
&help unless ($vcf && $bulk1 && $bulk2 && $key && $outdir);

sub help
{
	print "Function: calculate (snpindex, g-statistic, ed, fisher exact test(fet))
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
";
	exit;
}
############# start time #################
my $current_T = &date_format(localtime());
print "Programe start: $current_T\n\n";
my $begin_time = time();
##########################################
$vcf = &abs_dir($vcf);
my %sim;
my $simhead;
my @mutsim;
#depth	b1size	b2size	CI_95	CI_99
if(defined $sim_file){
	$sim_file = &abs_dir($sim_file);
	open(IN, $sim_file);
	while (<IN>) {
		chomp;
		my ($depth, $b1size, $b2size, @cis) = split(/\t/, $_);
		if(/depth/){
			$simhead = join("\t", @cis);
			next;
		}else{
			push(@{$sim{$depth}}, @cis);
		}
	}
	close(IN);
	if(defined $sim_depth && $sim_depth ne "min" && $sim_depth ne "max"){
		if(!exists $sim{$sim_depth}){
			$sim_depth = "max";
		}
	}else{
		$sim_depth ||= "max";
	}
	my @simh = split(/\t/, $simhead);
	for (my $i = 0; $i < @simh; $i++) {
		my @nindex = split(/\_/, $simh[$i]);
		my $lasti = scalar(@nindex) - 2;
		my $tmp = join("_", @nindex[0..$lasti]);
		if($tmp eq "A_Index"){
			$simh[$i] = "$bulk1\_Index\_$nindex[$lasti+1]";
		}elsif($tmp eq "B_Index"){
			$simh[$i] = "$bulk2\_Index\_$nindex[$lasti+1]";
		}
	}
	$simhead = join("\t", @simh);
	if(defined $parent1 && $parent1 eq $bulk1){ ## for mutmap
		my @muthead;
		my @cis = split(/\t/, $simhead);
		for (my $i = 0; $i < @cis; $i++) {
			my @nindex = split(/\_/, $cis[$i]);
			my $lasti = scalar(@nindex) - 2;
			my $tmp = join("_", @nindex[0..$lasti]);
			if($tmp eq "$bulk2\_Index"){
				push(@mutsim, $i);
				push(@muthead, $cis[$i]);
			}
		}
		$simhead = join("\t", @muthead);
	}
}
$sim_depth ||= "max";

if(!-d $outdir){
	mkdir $outdir;
}
$outdir = &abs_dir($outdir);

$popt ||= "ril";
$popt = lc($popt);
$min_depth ||= 10;
$max_depth ||= 500;
$min_index //= 0.2;
$pmin_depth ||= 5;
$indel ||= "F";
$abs ||= "F";
$method ||= "all";
$edpower ||= 5;
$rbin ||= "Rscript";
####
my $b1i;
my $b2i;
my $p1i;
my $p2i;
my $pflag = 0;
my %posh;
my $head = "chr\tpos\tref\talt\ttype";
my $outfile = "$outdir/$key\.all.xls";
my $fin;
if($vcf =~ /gz$/){
	open($fin, "gunzip -dc $vcf|");
}else{
	open($fin, $vcf);
}
open(DROP, ">$outdir/$key\.drop.vcf");
open(INDEX, ">$outfile");
open(BASIC, ">$outdir/$key\.basic.xls");
my $oindex;
if($method eq "all" || $method =~ /index/){
	open($oindex, ">$outdir/$key\.index.xls");
	print $oindex "chr\tpos";
}
my $ogst;
if($method eq "all" || $method =~ /gst/){
	open($ogst, ">$outdir/$key\.gst.xls");
	print $ogst "chr\tpos\tG\n";
}
my $oed;
if($method eq "all" || $method =~ /ed/){
	open($oed, ">$outdir/$key\.ed.xls");
	print $oed "chr\tpos\tED\n";
}
my $ofet;
if($method eq "all" || $method =~ /fet/){
	open($ofet, ">$outdir/$key\.fet.xls");
	print $ofet "chr\tpos\tP\t-log10\(P\)\n";
}
while (<$fin>) {
	chomp;
	if(/^##/){
		print DROP "$_\n";
		next;
	}
	my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split(/\t/, $_);
	my $sample_drop = join("\t", @samples);
	if(/^#CHROM/){
		print DROP "$_\n";
		for (my $i = 0; $i < @samples; $i++) {
			$b1i = $i if($samples[$i] eq $bulk1);
			$b2i = $i if($samples[$i] eq $bulk2);
			$p1i = $i if(defined $parent1 && $samples[$i] eq $parent1);
			$p2i = $i if(defined $parent2 && $samples[$i] eq $parent2);
		}
		if(!defined $b1i){
			die("ERROR: bulk1 $bulk1 was not found in vcf file, please check!\n");
		}
		if(!defined $b2i){
			die("ERROR: bulk2 $bulk2 was not found in vcf file, please check!\n");
		}
		if(defined $parent1 && !defined $p1i){
			die("ERROR: parent1 $parent1 was not found in vcf file, please check!\n");
		}
		if(defined $parent2 && !defined $p2i){
			die("ERROR: parent2 $parent2 was not found in vcf file, please check!\n");
		}

		if(defined $p1i){
			$pflag = 1;
		}
		if(defined $p2i){
			$pflag = 2;
		}
		if(defined $p1i && defined $p2i){
			$pflag = 3;
		}
		if(defined $parent1 && $parent1 eq $bulk1){ ## for mutmap
			$head .= "\t$samples[$p1i]";
		}else{ ## bsa
			if($pflag == 1){
				$head .= "\t$samples[$p1i]";
			}elsif($pflag == 2){
				$head .= "\t$samples[$p2i]";
			}elsif($pflag == 3){
				$head .= "\t$samples[$p1i]\t$samples[$p2i]";
			}
			$head .= "\t$samples[$b1i]\_GT\t$samples[$b1i]\_AD";
		}
		$head .= "\t$samples[$b2i]\_GT\t$samples[$b2i]\_AD";
		print BASIC "$head\n";
		if($method eq "all" || $method =~ /index/){
			if(defined $parent1 && $parent1 eq $bulk1){ ## for mutmap
				$head .= "\t$samples[$b2i]\_snpindex";
				print $oindex "\t$samples[$b2i]\_snpindex";
			}else{ ## for bsa
				$head .= "\t$samples[$b1i]\_snpindex\t$samples[$b2i]\_snpindex";
				$head .= "\tDindex";
				print $oindex "\t$samples[$b1i]\_snpindex\t$samples[$b2i]\_snpindex";
				print $oindex "\tDindex";
			}
			if(defined $sim_file){
				$head .= "\t$simhead";
				print $oindex "\t$simhead\n";
			}else{
				print $oindex "\n";
			}
		}
		if($method eq "all" || $method =~ /gst/){
			$head .= "\tG";
		}
		if($method eq "all" || $method =~ /ed/){
			$head .= "\tED";
		}
		if($method eq "all" || $method =~ /fet/){
			$head .= "\tP";
			$head .= "\t-log10(P)";
		}
		print INDEX "$head\n";
		next;
	}
	## remove duplicate sites
	if(exists $posh{"$chr-$pos"}){
		$filter = "posDup";
		print DROP "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$format\t$sample_drop\n";
		next;
	}else{
		$posh{"$chr-$pos"} = 1;
	}

	## remove none bi-allelic
	if($alt =~ /,/){
		$filter = "Mallelic";
		print DROP "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$format\t$sample_drop\n";
		next;
	}

	## remove or keep indels
	my $type = &judge_type($ref, $alt);
	if($indel eq "F" && $type ne "SNP"){
		print DROP "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$type\t$info\t$format\t$sample_drop\n";
		next;
	}

	## get format
	my ($gti, $adi) = &field_format($format);

	## get info
	my %info;
	&get_info($gti, $adi, \@samples, \%info);

	## filter missing
	my $miss_state = &filter_missing(\%info);
	if($miss_state ne "pass"){
		print DROP "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$miss_state\t$info\t$format\t$sample_drop\n";
		next;
	}

	## filter depth
	my $depth_state = &filter_depth(\%info);
	if($depth_state ne "pass"){
		print DROP "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$depth_state\t$info\t$format\t$sample_drop\n";
		next;
	}

	## filter genotype
	my $genotype_state = &filter_genotype(\%info);
	if($genotype_state ne "pass"){
		print DROP "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$genotype_state\t$info\t$format\t$sample_drop\n";
		next;
	}

	## get allele
	my ($b1r, $b2r, $b1sum, $b2sum) = &get_allele(\%info);
	## snp index calculate and filter
	my ($minDP, $maxDP, $b1index, $b2index, $delta, $index_state) = &snpindex($b1r, $b2r, $b1sum, $b2sum);
	if($index_state ne "pass"){
		print DROP "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$index_state\t$info\t$format\t$sample_drop\n";
		next;
	}
	## out info
	my $outline = "$chr\t$pos\t$ref\t$alt\t$type";
	if(defined $parent1 && $parent1 eq $bulk1){ ## for mutmap
		$outline .= "\t$info{p1}{g}\t$info{b2}{g}\t$info{b2}{a}";
	}else{ ## for bsa
		if($pflag == 0){
			$outline .= "\t$info{b1}{g}\t$info{b1}{a}\t$info{b2}{g}\t$info{b2}{a}";
		}elsif($pflag == 1){
			$outline .= "\t$info{p1}{g}\t$info{b1}{g}\t$info{b1}{a}";
			$outline .= "\t$info{b2}{g}\t$info{b2}{a}";
		}elsif($pflag == 2){
			$outline .= "\t$info{p2}{g}\t$info{b1}{g}\t$info{b1}{a}";
			$outline .= "\t$info{b2}{g}\t$info{b2}{a}";
		}elsif($pflag == 3){
			$outline .= "\t$info{p1}{g}\t$info{p2}{g}";
			$outline .= "\t$info{b1}{g}\t$info{b1}{a}\t$info{b2}{g}\t$info{b2}{a}";
		}
	}

	$outline =~ s/\//\|/g;
	print BASIC "$outline\n";
	## keep snp index or not
	if($method eq "all" || $method =~ /index/){
		if($abs eq "T"){
			$delta = abs($delta);
		}
		if(defined $parent1 && $parent1 eq $bulk1){ ## for mutmap
			$outline .= "\t$b2index";
			print $oindex "$chr\t$pos\t$b2index";
		}else{ ## for bsa
			$outline .= "\t$b1index\t$b2index\t$delta";
			print $oindex "$chr\t$pos\t$b1index\t$b2index\t$delta";
		}
		if(defined $sim_file){
			my $simout;
			if($sim_depth eq "min"){
				$simout = join("\t", @{$sim{$minDP}});
			}elsif($sim_depth eq "max"){
				$simout = join("\t", @{$sim{$maxDP}});
			}else{
				$simout = join("\t", @{$sim{$sim_depth}});
			}
			if(defined $parent1 && $parent1 eq $bulk1){ ## for mutmap
				my @cis = split(/\t/, $simout);
				my @cisout;
				for (my $i = 0; $i < @mutsim; $i++) {
					push(@cisout, $cis[$mutsim[$i]]);
				}
				$simout = join("\t", @cisout);
			}
			$outline .= "\t$simout";
			print $oindex "\t$simout\n";
		}else{
			print $oindex "\n";
		}
	}
	## g statistic
	if($method eq "all" || $method =~ /gst/){
		my $gstat = &get_g($b2r, $b1r, $b2sum - $b2r, $b1sum - $b1r);
		$outline .= "\t$gstat";
		print $ogst "$chr\t$pos\t$gstat\n";
	}
	## ed value
	if($method eq "all" || $method =~ /ed/){
		my $edstat = &get_ed($b1r, $b2r, $b1sum, $b2sum);
		$outline .= "\t$edstat";
		print $oed "$chr\t$pos\t$edstat\n";
	}
	## fisher test
	if($method eq "all" || $method =~ /fet/){
		#sub fisher_exact_test { # Fisher's Exact test can just for 2 x 2 table.
		#          word2   ~word2
		#  word1    n11      n12 | n1p
		# ~word1    n21      n22 | n2p
		#           --------------
		#           np1      np2   npp
		# my ( $n11, $n1p, $np1, $npp ) = @_;
		# return ($pValue, $statisticName);
		my $n11 = $b1r;
		my $np1 = $b1r + $b2r;
		my $n1p = $b1sum;
		my $npp = $b1sum + $b2sum;
		my ($fetp, $stat_name) = &fisher_exact_test($n11, $n1p, $np1, $npp);
		$outline .= "\t$fetp";
		my $log10p = -log10($fetp);
		if($log10p < 0){
			$log10p = 0;
		}
		$outline .= "\t$log10p";
		print $ofet "$chr\t$pos\t$fetp\t$log10p\n";
	}
	## output all
	print INDEX "$outline\n";
}
close($fin);
close(DROP);
close(INDEX);
close(BASIC);

if($method eq "all" || $method =~ /index/){
	close($oindex);
}
if($method eq "all" || $method =~ /gst/){
	close($ogst);
}
if($method eq "all" || $method =~ /ed/){
	close($oed);
}
if($method eq "all" || $method =~ /fet/){
	close($ofet);
	my $fetr = "$outdir/$key\.fet.r";
	open(RS, ">$fetr");
	print RS"setwd(\"$outdir\")
myfet <- read.table(\"$outdir/$key\.fet.xls\", header = T, check.names = F)
myfet\$Q <- p.adjust(myfet\[,3\], \"BH\")
write.table(myfet, \"$outdir/$key\.fet.xls\", row.names = F, col.names = T, quote = F, sep = \"\\t\")
";
	close(RS);
	`$rbin $fetr`;
}

############# end time ###################
$current_T = &date_format(localtime());
print "Programe end: $current_T\n\n";
&runtime($begin_time);
##########################################
# log10
sub log10 {
	my $n = shift;
	return log($n)/log(10);
}
# judge type
sub judge_type()
{
	my ($ref, $alt) = @_;
	my $type = "SNP";
	if($ref =~ /\w{2,}|\-/ || $alt =~ /\w{2,}|\-/){
		my $rlen = length($ref);
		my $alen = length($alt);
		if($ref =~ /\-/ || $rlen < $alen){
			$type = "Insertion";
		}else{
			$type = "Deletion";
		}
	}
	return $type;
}
# gti and adi
sub field_format()
{
	my ($fmat) = @_;
	my ($gti, $adi);
	my @fm = split(/:/, $fmat);
	for (my $i = 0; $i < @fm; $i++) {
		if($fm[$i] eq "GT") {
			$gti = $i;
			next;
		}
		if($fm[$i] eq "AD"){
			$adi = $i;
		}
	}
	if(!defined $gti){
		&show_log("ERROR: GT is not included in FORMAT, please check!");
		die;
	}
	if(!defined $adi){
		&show_log("ERROR: AD is not included in FORMAT, please check!");
		die;
	}
	return ($gti, $adi);
}

# get info
sub get_info()
{
	my ($gti, $adi, $samples, $info) = @_;
	my @b1info = split(/:/, $$samples[$b1i]);
	my @b2info = split(/:/, $$samples[$b2i]);
	$info->{b1}->{g} = $b1info[$gti];
	$info->{b1}->{a} = $b1info[$adi];
	$info->{b2}->{g} = $b2info[$gti];
	$info->{b2}->{a} = $b2info[$adi];
	if($pflag == 1 || $pflag == 3){
		my @p1info = split(/:/, $$samples[$p1i]);
		$info->{p1}->{g} = $p1info[$gti];
		$info->{p1}->{a} = $p1info[$adi];
	}
	if($pflag == 2 || $pflag == 3){
		my @p2info = split(/:/, $$samples[$p2i]);
		$info->{p2}->{g} = $p2info[$gti];
		$info->{p2}->{a} = $p2info[$adi];
	}
}
# missing filter
sub filter_missing()
{
	my ($info) = @_;
	my $state = "pass";
	if($info->{b1}->{g} =~ /\./){
		$state = "B1miss";
		return $state;
	}
	if($info->{b2}->{g} =~ /\./){
		$state = "B2miss";
		return $state;
	}
	if($pflag == 1){
		if($info->{p1}->{g} =~ /\./){
			$state = "P1miss";
		}
	}elsif($pflag == 2){
		if($info->{p2}->{g} =~ /\./){
			$state = "P2miss";
		}
	}elsif($pflag == 3){
		if($info->{p1}->{g} =~ /\./){
			$state = "P1miss";
			return($state);
		}
		if($info->{p2}->{g} =~ /\./){
			$state = "P2miss";
			return($state);
		}
	}
	return $state;
}

# filter depth
sub filter_depth()
{
	my ($info) = @_;
	my @b1ad = split(/,/, $info->{b1}->{a});
	my @b2ad = split(/,/, $info->{b2}->{a});
	my $b1sum = $b1ad[0];
	my $b2sum = $b2ad[0];
	$b1sum += $b1ad[1];
	$b2sum += $b2ad[1];
	my $state = "pass";
	## bulk depth filter
	if(defined $parent1 && $parent1 eq $bulk1){ ## for mutmap
		# nothing	
	}else{ ## for bsa
		if($b1sum < $min_depth){
			$state = "B1depthL";
			return ($state);
		}elsif($b1sum > $max_depth){
			$state = "B1depthH";
			return ($state);
		}
	}

	if($b2sum < $min_depth){
		$state = "B2depthL";
		return ($state);
	}elsif($b2sum > $max_depth){
		$state = "B2depthH";
		return ($state);
	}
	## parents depth filter
	my ($p1depth, $p2depth);
	if($pflag == 0){
		$p1depth = $pmin_depth * 2;
		$p2depth = $p1depth;
	}elsif($pflag == 1){
		my @p1ad = split(/,/, $info->{p1}->{a});
		$p1depth = $p1ad[0] + $p1ad[1];
		$p2depth = $pmin_depth * 2;
	}elsif($pflag == 2){
		my @p2ad = split(/,/, $info->{p2}->{a});
		$p2depth = $p2ad[0] + $p2ad[1];
		$p1depth = $pmin_depth * 2;
	}elsif($pflag == 3){
		my @p1ad = split(/,/, $info->{p1}->{a});
		$p1depth = $p1ad[0] + $p1ad[1];
		my @p2ad = split(/,/, $info->{p2}->{a});
		$p2depth = $p2ad[0] + $p2ad[1];
	}
	if($p1depth < $pmin_depth){
		$state = "P1depthL";
		return($state);
	}elsif($p1depth > $max_depth){
		$state = "P1depthH";
		return($state);
	}
	if($p2depth < $pmin_depth){
		$state = "P2depthL";
		return($state);
	}elsif($p2depth > $max_depth){
		$state = "P2depthH";
		return($state);
	}
	return ($state);
}

# filter genotype
sub filter_genotype()
{
	my ($info) = @_;
	my $state = "pass";
	my @b1gt = split(/\/|\|/, $info->{b1}->{g});
	my @b2gt = split(/\/|\|/, $info->{b2}->{g});
	my %allele;
	$allele{$b1gt[0]}++;
	$allele{$b1gt[1]}++;
	$allele{$b2gt[0]}++;
	$allele{$b2gt[1]}++;
	my @geno = sort{$allele{$a} <=> $allele{$b}} keys %allele;
	if(@geno != 2){
		$state = "Oallelic";
		return $state;
	}
	## parents genotype
	my ($p1, $p2, $p3, $p4);
	if($pflag == 1){
		($p1, $p2) = split(/\/|\|/, $info->{p1}->{g});
		if($p1 ne $geno[0]){
			$p3 = $geno[0];
			$p4 = $geno[0];
		}else{
			$p3 = $geno[1];
			$p4 = $geno[1];
		}
	}elsif($pflag == 2){
		($p3, $p4) = split(/\/|\|/, $info->{p2}->{g});
		if($p3 ne $geno[0]){
			$p1 = $geno[0];
			$p2 = $geno[0];
		}else{
			$p1 = $geno[1];
			$p2 = $geno[1];
		}
	}elsif($pflag == 3){
		($p1, $p2) = split(/\/|\|/, $info->{p1}->{g});
		($p3, $p4) = split(/\/|\|/, $info->{p2}->{g});
	}
	if($pflag > 0){
		if($popt ne "f1" && $popt ne "cp"){
			if($p1 eq $p3 || $p1 ne $p2 || $p3 ne $p4){
				$state = "Pgenotype";
			}
		}else{
			if($p1 eq $p2 && $p3 eq $p4){
				$state = "Pgenotype";
			}
		}
	}
	return $state;
}
# sub get alleles
sub get_allele()
{
	my ($info) = @_;
	my @b1ad = split(/,/, $info->{b1}->{a});
	my @b2ad = split(/,/, $info->{b2}->{a});
	my $b1sum = $b1ad[0] + $b1ad[1];
	my $b2sum = $b2ad[0] + $b2ad[1];
	my $recessive_pos = 1; ## if no parents, alts are recessive alleles
	if($popt ne "f1" && $popt ne "cp"){
		if($pflag == 1){
			if($info->{p1}->{g} eq "1/1" || $info->{p1}->{g} eq "1|1"){
				$recessive_pos = 0;
			}
		}elsif($pflag > 1){
			if($info->{p2}->{g} eq "0/0" || $info->{p2}->{g} eq "0|0"){
				$recessive_pos = 0;
			}
		}
	}else{
		if($pflag == 1){
			if($info->{p1}->{g} eq "1/1" || $info->{p1}->{g} eq "1|1"){# lmxll ?
				$recessive_pos = 0;
			} # else (nnxnp & hkxhk ?)
		}elsif($pflag == 2){
			if($info->{p2}->{g} eq "0/0" || $info->{p2}->{g} eq "0|0"){# lmxll ?
				$recessive_pos = 0;
			} #else (nnxnp & hkxhk ?)
		}elsif($pflag == 3){
			my ($p1, $p2) = split(/\/|\|/, $info->{p1}->{g});
			my ($p3, $p4) = split(/\/|\|/, $info->{p2}->{g});
			if($p1 eq $p2 && $p1 eq $p4){ ## nnxnp 1/1 0/1
				$recessive_pos = 0;
			}elsif($p1 ne $p2 && $p1 eq $p3){ ## lmxll 0/1 0/0
				$recessive_pos = 0;
			}
		}
	}
	my $b1r = $b1ad[$recessive_pos];
	my $b2r = $b2ad[$recessive_pos];
	return($b1r, $b2r, $b1sum, $b2sum);
	# r: recessive
	# d: dominant
}

# snp index calculate
sub snpindex()
{
	my ($b1r, $b2r, $b1sum, $b2sum) = @_;
	my $state = "pass";
	my ($b1index, $b2index, $delta);
	if($b1sum > 0 && $b2sum > 0){
		$b1index = $b1r/$b1sum;
		$b2index = $b2r/$b2sum;
	}else{
		$b1index = 0;
		$b2index = 0;
	}
	$delta = $b2index - $b1index;
	if($b1index < $min_index && $b2index < $min_index){
		$state = "IndexL";
	}elsif($b1index == 1 && $b2index == 1){
		$state = "Index1";
	}elsif($b1index > 1 - $min_index && $b2index > 1 - $min_index){
		if(defined $parent2 && $parent2 eq $bulk2){ ## for mutmap
			#nothing
		}else{ ## for bsa
			$state = "IndexH";
		}
	}
	my ($minDP, $maxDP);
	if($b1sum > $b2sum){
		$minDP = $b2sum;
		$maxDP = $b1sum;
	}else{
		$minDP = $b1sum;
		$maxDP = $b2sum;
	}
	return ($minDP, $maxDP, $b1index, $b2index, $delta, $state);
}
# Compute G value
# Return: ($G)
sub get_g()
{
	# r: recessive
	# d: dominant
	my ($b2r, $b1r, $b2d, $b1d) = @_;
	my $sum = $b2r + $b1r + $b2d + $b1d;
	$sum > 0 or return('NA');
	my $np1 = ($b2r+$b1r)*($b2r+$b2d)/$sum;
	my $np2 = ($b2r+$b1r)*($b1r+$b1d)/$sum;
	my $np3 = ($b2d+$b1d)*($b2r+$b2d)/$sum;
	my $np4 = ($b1d+$b2d)*($b1d+$b1r)/$sum;
	$np1 > 0 or return('NA');
	$np2 > 0 or return('NA');
	$np3 > 0 or return('NA');
	$np4 > 0 or return('NA');
	my $v1 = ( $b2r > 0 ) ? $b2r * log($b2r/$np1) : 0;
	my $v2 = ( $b1r > 0 ) ? $b1r * log($b1r/$np2) : 0;
	my $v3 = ( $b2d > 0 ) ? $b2d * log($b2d/$np3) : 0;
	my $v4 = ( $b1d > 0 ) ? $b1d * log($b1d/$np4) : 0;
	#my $g = 2 * ($v1+$v2+$v3+$v4)/$sum;
	my $g = 2 * ($v1+$v2+$v3+$v4);
	return($g);
}# getG() 

sub get_ed()
{
	my ($b1r, $b2r, $b1sum, $b2sum) = @_; #
	my $b1rf = $b1r/$b1sum;
	my $b2rf = $b2r/$b2sum;
	my $b1df = 1 - $b1rf;
	my $b2df = 1 - $b2rf;
	my $edv = sqrt(($b1rf-$b2rf) ** 2 + ($b1df-$b2df) ** 2);
	$edv = $edv ** $edpower;
	return($edv);
}

# fisher exact test
sub fisher_exact_test { # Fisher's Exact test can just for 2 x 2 table.
#          word2   ~word2
#  word1    n11      n12 | n1p
# ~word1    n21      n22 | n2p
#           --------------
#           np1      np2   npp

	my ( $n11, $n1p, $np1, $npp ) = @_;
	my ($pValue, $statisticName, $errorCode); 

	# ($n11, $n1p, $np1, $npp) all have to be interge when Fisher's exact test
	$pValue = Text::NSP::Measures::2D::Fisher::twotailed::calculateStatistic(
		n11=>$n11,n1p=>$n1p,np1=>$np1,npp=>$npp);
	if( ($errorCode = Text::NSP::Measures::2D::Fisher::twotailed::getErrorCode()) ) {
		print STDERR $errorCode." - ".Text::NSP::Measures::2D::Fisher::twotailed::getErrorMessage();
		return ( 1, "-" ); # Return this value when Fail Test.
	}
	$statisticName = Text::NSP::Measures::2D::Fisher::twotailed::getStatisticName; 
	$statisticName = join "",(split(/\s+/,$statisticName));

	return ($pValue, $statisticName);
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
