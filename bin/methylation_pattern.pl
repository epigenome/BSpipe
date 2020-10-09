#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Author: Jeong-Hyeon Choi
## Date: 2010

## Description

## Usage:
## CHECK usage with option  -h

## $Id: perl_template,v 1.2 2007-09-21 22:51:23 jeochoi Exp $

my %RC=(
			"A"=>"T",
			"T"=>"A",
			"C"=>"G",
			"G"=>"C",
			"M"=>"K",
			"K"=>"M",
			"R"=>"Y",
			"Y"=>"R",
			"W"=>"S",
			"S"=>"W",
			"B"=>"V",
			"D"=>"H",
			"H"=>"D",
			"V"=>"B",
			"N"=>"N"
);
			
use Getopt::Long qw(:config no_ignore_case);

my ($helpFlag, $inFile, $outFile, $verbose, $quiet, $samFlag, $posFile, $bedFlag, $quality, $depth, $minCpgs, $readLen, $pattern, $minMet, $maxMet);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"sam"				=> \$samFlag,
	"length=i"		=> \$readLen,
	"p|pos=s"		=> \$posFile,
	"bed"				=> \$bedFlag,
	"q|quality=i"	=> \$quality,
	"depth=i"		=> \$depth,
	"s|site=i"		=> \$minCpgs,
	"c|context=s"	=> \$pattern,
	"min=f"			=> \$minMet,
	"max=f"			=> \$maxMet,
) || die "\n";

checkOptions();

my (%tars, %wat, %wbeg, %wend, %cre, %cbeg, %cend, %reads, $out, $mode, @buf);
my (%pats);
my ($lowQuality) = (2);
my ($rno, $ref, $start) = (0, '', -1);

my @bps = context($pattern);
load($posFile);
process($inFile);


#-------------------------------------------------------------------------------

=format
HWI-EAS407_0120:3:30:8203:2363#0/1      0       chr10   63995   255     51M     *       0       0       CGGGGTAGGGTAGTTATTTTATTTGTGGTTTGAATTGGTCGTTTTA
GTTTG   GGGGGBGGGG@CFEFGGGGGDGGGGCFGCGGGGGGHHDEDGDGGG>GGGGH     XA:i:0  MD:Z:51 ZD:Z:51 NM:i:0  ZE:i:0  ZT:Z:C  ZM:i:2  ZU:i:11 YM:Z:1,40       YU:Z:1
1,14,15,17,19,23,39,42,44,48,49
=cut
sub process
{
	my ($fileName) = @_;

	print "Processing ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput(!$fileName || $samFlag ? $fileName : "samtools view $fileName |");
	$out = openOutput($outFile);
	print $out "#Ref\tStart\tEnd\tStrand\tNcpgs\tMethylation\tNreads\tNpats\tNchange\tPattern\n";

	while (<$in>)
	{
		next if /^(#|\@)/ || /^\s*$/;
		@buf = split(/\t/);
		my (@mp, @up);
		$start = $buf[3] if $start == -1;

		die "Unknowm mapping mode: $_" if !/ZT:Z:(\w)/;
		$mode = $1 eq 'G' ? 0 : 1;

		if ($ref && ($ref ne $buf[2] || ($start != -1 && $start+$readLen < $buf[3])))
		{
			die if $start == -1;
			print "====================================\n" if $verbose;
			print "START=$start CURRENT=$buf[3]\n" if $verbose;
			my $ls = onWindow(\%wat, \%wbeg, \%wend, '+', $ref eq $buf[2]);
			my $rs = onWindow(\%cre, \%cbeg, \%cend, '-', $ref eq $buf[2]);
			$start = $ls == -1 ? $rs : $rs == -1 ? $ls : $ls < $rs ? $ls : $rs;
		}

		$rno++;
		my @indels = indel($buf[5]);
		add(1, \@indels, split(/,/, $1)) if /YM:Z:(\S+)/;
		add(0, \@indels, split(/,/, $1)) if /YU:Z:(\S+)/;

		if (exists $reads{$rno})
		{
			my @b = sort {$a <=> $b} keys %{$reads{$rno}};
			if ($mode) { push(@{$wbeg{$b[0]}}, $rno); push(@{$wend{$b[-1]}}, $rno); }
			else       { push(@{$cbeg{$b[0]}}, $rno); push(@{$cend{$b[-1]}}, $rno); }
		}

		$ref = $buf[2];
	}

		onWindow(\%wat, \%wbeg, \%wend, '+', 0) while (%wat);
		onWindow(\%cre, \%cbeg, \%cend, '-', 0) while (%cre);

	close($in) if defined $fileName;
	close($out) if defined $outFile;
}

=info
         C  G
CG  CG   0  1
CH  DG   0  1
CHG CDG  0  2
CHH DDG  0  2
GCH DGC  1  1
HCG CGD  1  1
GCG CGC  1  1
=cut
sub context
{
	my ($pat) = @_;
	$pat = uc $pat;
	my $prc = join('', map { die "Unknown IUPAC base: $_\n" if !exists $RC{$_}; $RC{$_} } reverse(split(//, $pat)));
	my $i = index($pat, 'C');
	my $j = index($prc, 'G');
	print STDERR "Pattern: [$pat] $i   [$prc] $j\n";
	die "Not C in $pat\n" if $i == -1;
	die "Not G in $pat\n" if $j == -1;
	return (-$j, -$i);
}

sub indel
{
	my ($s) = @_;
	my $chars = "MDISHN=X";
	my (@t);

	die "Unsupported cigar character $& in $s\n" if $s =~ /[^$chars\d]/;
	while ($s =~ /(\d*)([$chars])/g) { push(@t, $1, $2); }

	my ($p, $q) = (0, 0);
	my @d = (0, 0);
	for (my $i = 1; $i < @t; $i += 2)
	{
		if    ($t[$i] eq 'H'                 ) { }
		elsif ($t[$i] eq 'I' || $t[$i] eq 'S') { push(@d, $p,  $t[$i-1]); }
		elsif ($t[$i] eq 'D' || $t[$i] eq 'N') { push(@d, $p, -$t[$i-1]); $q += $t[$i-1]; }
		else                                   { $p += $t[$i-1]; }
	}
	push(@d, $p+$q, '');

	for (my $i = 3; $i < $#d; $i += 2) { $d[$i] += $d[$i-2]; }
	print "$s (@d) [$p $q]\n" if $verbose;
	return @d;
}

sub add
{
	my ($met, $rindels, @ps) = @_;
	my $i = 2; # for indels

	for (my $j = 0; $j < @ps; $j++)
	{
		while ($ps[$j] > $rindels->[$i]) { $i += 2; }
		my $p = $buf[3] + $ps[$j] - 1; # for GCH, WGC HCG GCC CGG
#		if (exists $tars{$buf[2]}{$mode ? $p : $p-1})
		if (exists $tars{$buf[2]}{$p + $bps[$mode]})
		{
			print "@buf$buf[10], $ps[$j]+$rindels->[$i-1]-1, 1\n" if length($buf[10]) <= $ps[$j]+$rindels->[$i-1]-1;
			my $m = ord(substr($buf[10], $ps[$j]+$rindels->[$i-1]-1, 1)) < $quality ? $lowQuality : $met;
			if ($mode) { push(@{$wat{$p}[$m]}, $rno); }
			else       { push(@{$cre{$p}[$m]}, $rno); }
			$reads{$rno}{$p} = $m;
		}
	}
}

sub onWindow
{
	my ($rsit, $rbeg, $rend, $strand, $flag) = @_;
	my $last = $flag ? $buf[3] : 1000000000; # process all remain sites if flag is 0
	my $_start;

	foreach my $p (sort {$a<=> $b} keys %$rbeg)
	{
		if ($p + $readLen >= $last)
		{
			$_start = $p;
			last;
		}

		foreach my $q (sort {$a<=>$b} keys %$rend)
		{
			next if $p == $q;
			last if $p + $readLen < $q;
			onInterval($rsit, $p, $q, $strand);
		}
	}

	map { delete $rsit->{$_} if $_+$readLen < $last } keys %$rsit;
	map { delete $rbeg->{$_} if $_+$readLen < $last } keys %$rbeg;
	map {
		if ($_+$readLen < $last)
		{
			map { delete $reads{$_} } @{$rend->{$_}};
			delete $rend->{$_}
		}
	} keys %$rend;

	print "$start => ", $_start || -1, "\n\n" if $verbose;
	return $_start || -1;
}

sub onInterval
{
	my ($rsit, $beg, $end, $strand) = @_;

	my (@cpgs, @ps, %rids);
	map {push(@cpgs, $_) if $beg <= $_ && $_ <= $end} keys %$rsit;
	print "$ref $beg $end ", @cpgs+0, "\n" if $verbose;
	return if @cpgs < $minCpgs;

	@cpgs = sort {$a<=>$b} @cpgs;
	print "I: [$beg $end] @cpgs\n" if $verbose;
	my ($msum, $mcnt) = (0, 0); # mcnt for the number of CpGs with minMet <= met <= maxMet

	foreach my $p (@cpgs)
	{
		my $u = $rsit->{$p}[0] ? @{$rsit->{$p}[0]} : 0;
		my $m = $rsit->{$p}[1] ? @{$rsit->{$p}[1]} : 0;
		if ($u+$m >= $depth)
		{
			push(@ps, $p);
			my $beta = $m / ($u+$m);
			$msum += $beta;
			$mcnt ++ if $minMet <= $beta && $beta <= $maxMet;
		}
		print "\t$p $u $m $msum\n" if $verbose;
	}

	print "PS: @ps\n" if $verbose;
	return if (@ps < $minCpgs);
	
	map { map { $rids{$_}++ } @$_ } (@{$rsit->{$beg}}[0,1], @{$rsit->{$end}}[0,1]);
	map { delete $rids{$_} if $rids{$_} != 2 } keys %rids;
	if (keys %rids < $depth)
	{
		print "\tless reads\n" if $verbose;
		return;
	}

	%pats = ();
	my ($nread, $npat, $nsit) = entropy(\@ps, \%rids);
	my $s = sprintf "%s\t%d\t%d\t%s\t%d\t%.3f\t%d\t%d\t%d", $ref, $beg, $end, $strand, @ps+0, $msum/@ps, $nread, $npat, $nsit;
	print $out $s, "\t", join(":", map{"$_=$pats{$_}"} keys %pats), "\n" if keys %pats;
}

sub entropy
{
	my ($rcpgs, $rrids) = @_;
	my (%sites, %gpos);
	my ($count, $gap) = (0, '-');

	foreach my $r (keys %$rrids)
	{
		my ($p, $valid) = ('', 0);
		map {
			if    (!exists $reads{$r}{$_}       ) { $gpos{length($p)}++; $p .= $gap; }
			elsif ($reads{$r}{$_} == $lowQuality) { $gpos{length($p)}++; $p .= $gap; }
			else                                  { $valid++;            $p .= $reads{$r}{$_}; }
			$sites{$_}{substr($p, length($p)-1, 1)}++;
		} @$rcpgs;
		next if $valid < $minCpgs;

		$pats{$p}++;
		$count++;
	}

	my $change = 0;
	map {my $n = keys %$_; $n-- if exists $_->{$gap}; $change++ if $n > 1} values %sites;

	return ($count, keys(%pats)+0, $change);
}

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my @a = split;
		$tars{$a[0]}{$bedFlag ? $a[1]+1 : $a[1]}++ if !$a[3] || $a[3] >= $depth;
	}

	close($in) if defined $fileName;
}

#-------------------------------------------------------------------------------

sub openInput
{
	my ($fileName) = @_;

	return *STDIN unless defined $fileName;

	my ($fd);
	open($fd, $fileName =~ /.gz(ip)?$/ ? "zcat $fileName |" : $fileName =~ /.bz(ip)?2$/ ? "bzcat $fileName |" : $fileName) || die("Open error: $fileName");
	return $fd;
}

sub openOutput
{
	my ($fileName) = @_;

	return *STDOUT unless defined $fileName;

	my ($fd);
	open($fd, $fileName =~ /.gz$/ ? "| gzip -c > $fileName" : $fileName =~ /.bz(ip)?2$/ ? "| bzip2 -z -c > $fileName" : ">$fileName") || die("Open error: $fileName");
	return $fd;
}

sub checkOptions
{
	$inFile  = shift(@ARGV) if !defined $inFile  && @ARGV > 0;
	$outFile = shift(@ARGV) if !defined $outFile && @ARGV > 0;
	die "The output file name is the same as the input file name\n" if defined $inFile && defined $outFile && $inFile eq $outFile;

	if ($helpFlag || (!$posFile && !$inFile) || !$pattern)
	{
		die("Arguments: -c pattern [-bed] -p pos_file [-min min_met] [-max max_met] [-q base_quality=0] [-d depath=5] [-c cpg_count=3] [-l read_length] [-sam] [[-i] in_file] [[-o] out_file] [-cr change_rate -pr pattern_rate -pat pattern_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}

	$depth = 5 if !$depth;
	$minCpgs = 3 if !$minCpgs;
	$quality = 0 if !$quality;
	$readLen = 500 if !$readLen;
	$verbose = 0 if !$verbose;
	$minMet = 0 if !defined $minMet;
	$maxMet = 1 if !defined $maxMet;
}
