#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Author: Jeong-Hyeon Choi
## Date: 2010

## Description

## Usage:
## CHECK usage with option  -h

## $Id: perl_template,v 1.2 2007-09-21 22:51:23 jeochoi Exp $


use Getopt::Long qw(:config no_ignore_case);

my ($helpFlag, $inFile, $outFile, $verbose, $quiet, $samFlag, $window, $step, @posFiles, $depth, $minSite, @csites, @gsites);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"sam"				=> \$samFlag,
	"window=i"		=> \$window,
	"s|step=i"		=> \$step,
	"pos=s"			=> \@posFiles,
	"depth=i"		=> \$depth,
	"number=i"		=> \$minSite,
	"csite=i"		=> \@csites,
	"gsite=i"		=> \@gsites,
) || die "\n";

checkOptions();

my (%tars, %sites, %reads, @ps, $out, $mode, @a, $length);
my ($rno, $ref, $start) = (0, '', 1);

map {load($posFiles[$_], $_)} 0..$#posFiles;
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
	print $out "#Ref\tStart\tEnd\tNsites\tNreads\tNpats\tPosIndices\tPositions\tPattern=Count\n";

	while (<$in>)
	{
		next if /^(#|\@)/ || /^\s*$/;
		@a = split(/\t/);
		my (@mp, @up);

		die "Unknowm mapping mode: $_" if !/ZT:Z:(\w)/;
		$mode = $1;

		onWindow(1) while ($ref && ($ref ne $a[2] || $start+$window+$length < $a[3]));

		$rno++;
		add(1, split(/,/, $1)) if /YM:Z:(\S+)/;
		add(0, split(/,/, $1)) if /YU:Z:(\S+)/;

		$ref = $a[2];
	}

		onWindow(0) while (%sites);
	close($in) if defined $fileName;
}

sub add
{
	my $flag = shift;
	map {
		for (my $i = 0; $i < @csites; $i++)
		{
			my $p = $a[3] + $_ - (uc($mode) eq 'G' ? $gsites[$i] : $csites[$i]) - 2;
			if (exists $tars{$a[2]}{$p})
			{
				die "Unexpected case: $_\n" if $tars{$a[2]}{$p} != $i;
				print "$a[0] $p=$a[3]+$_-", (uc($mode) eq 'G' ? $gsites[$i] : $csites[$i]), "-2 $flag\n" if $verbose >= 2;
				my $index = $tars{$a[2]}{$p}*2 + $flag;
				$sites{$p}[$flag]++;
				$reads{$rno}{$p} = $index;
			}
		}
	} @_;
}

sub onWindow
{
	my ($flag) = @_;

	print "$ref:$start  =>  $_" if $verbose;
	@ps = ();
	my @sites = sort {$a<=>$b} keys %sites;
	print "Sites: @sites\n" if $verbose;
	my @counts;

	foreach my $p (@sites)
	{
		my $u = $sites{$p}[0]||0;
		my $m = $sites{$p}[1]||0;
		if ($p < $start+$window)
		{
			print "\t$p $tars{$ref}{$p} $u+$m\n" if $verbose > 1;
			if ($u+$m < $depth) { die "Unexpected case"; }
			else { push(@ps, $p); $counts[$tars{$ref}{$p}]++; }
			#print "$p $u $m $msum\n" if $verbose;
		}
	}

	if (@ps && @ps >= $minSite)
	{
		print "PS: @ps\n" if $verbose;
		print $out join("\t", $ref, $start, $start+$window-1, join(',', map {$_||0} @counts));
		pattern();
	}

	$start += $step;

	if (@ps)
	{
		map {delete $sites{$_} if $_ < $start} (@ps);
		map {my $i = $_;
			map {delete $reads{$i}{$_} if $_ < $start} (@ps);
			delete $reads{$i} if !keys %{$reads{$i}};
		} keys %reads;
	}
	elsif (@sites)
	{
		$start = int(($sites[0]-1)/$window)*$window+1;
		print "==========> skip to $start : @sites\n" if $verbose;
	}

	if (@sites == 0 && $flag) # || $ref ne $a[2])
	{
		$start = int(($a[3]-1)/$window)*$window+1;
		$ref = $a[2] if $ref ne $a[2];
		print "==========> skip to $start\n" if $verbose;
	}
}

sub pattern
{
	my ($count, %pats) = (0);

	foreach my $r (values %reads)
	{
		my $p = '';
		my $f;
		map {
			if (!exists $r->{$_}) { $p .=  '-'; }
			else                  { $p .= $r->{$_}; $f = 1;}
		} @ps;
		next if !$f;
		#print "$count $p\n" if $verbose;
		$pats{$p}++;
		$count++;
	}

	print $out "\t", join("\t", $count, keys(%pats)+0);
	print $out "\t", join(',', map {$tars{$ref}{$_}} @ps);
	print $out "\t", join(',', @ps);
	print $out "\t", join(',', map {"$_=$pats{$_}"} sort keys %pats), "\n";
}

sub load
{
	my ($fileName, $index) = @_;

	print "Loading $index: ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my @a = split;
		$length = $a[2]-$a[1] if !$length || $length < $a[2]-$a[1];
		$tars{$a[0]}{$a[1]} = $index if !$a[3] || $a[3] >= $depth;
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

	if ($helpFlag || !$window || !@csites || !@gsites || (!@posFiles && !$inFile))
	{
		die("Arguments: -p pos_file1 -csite pos1 -gsite pos1 [-p pos_file2 -csite pos2 -gsite pos2 ...] [-n min_sites] [-d depth] -w window [-s step] [-sam] [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}

	$step = $window if !$step;
	$verbose = 0 if !$verbose;
	$minSite = 0 if !$minSite;
	$depth = 0 if !$depth;

	die "Different number of parameter:', (@posFiles) (@csites)" if @posFiles != @csites;
	die "Different number of parameter:', (@posFiles) (@gsites)" if @posFiles != @gsites;
}
