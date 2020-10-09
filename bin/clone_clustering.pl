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
use File::Basename;
use IPC::Open2;

my $rscript = "/app/redhat6/R/R-3.6.1/bin/Rscript " . dirname($0) . "/hclust.r";
my ($helpFlag, $inFile, $accFile, $outFile, $verbose, $quiet);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
) || die "\n";

checkOptions();

process($inFile);


#-------------------------------------------------------------------------------

=format
#Ref  Start End   Nsites   Nreads   Npats PosIndices  Positions   Pattern=Count
chr1  10401 10600 3,1   12 4  0,0,0,1  10550,10565,10581,10587 -112=2,0002=4,0003=4,1112=2
chr1  13201 13400 4  19 13 0,0,0,0  13308,13320,13338,13353 00--=2,000-=1,0000=2,0011=1,01--=1,0111=1,10--=2,1000=3,11--=1,1101=1,111-=1,1110=2,1111=1
=cut
sub process
{
	my ($fileName) = @_;

	print "Processing ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	my $out = openOutput($outFile);
	print $out "#Ref\tStart\tEnd\tNreads\tNpatterns\tNclusters\tClusterSizes\tPatterns=Counts:ClusterNumbers\n";


	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my @a = split(/[\t\r\n]/);
		my @poss = split(/,/, $a[7]);
		my ($count, %pats) = (0);
		map {my @b = split(/=/); $pats{$b[0]} = $b[1]; $count += $b[1]} split(/,/, $a[8]);

		print "@a[0..2]\n" if $verbose;
		my (%distMat, $avgDist, %clus, @csizes);

		die "Different number of reads: $a[4] $count\n$_" if $a[4] != $count;
		distance_matrix(\%pats, \%distMat);
		$avgDist = average_distance($count, \%distMat);
		print_distance_matrix(\*STDOUT, \%distMat) if $verbose;
		clustering(\%pats, \%distMat, \%clus, \@csizes);
		print $out join("\t", @a[0..2,4,3], $#csizes, join(',', @csizes[1..$#csizes]), join(',', map {"$_=$pats{$_}:$clus{$_}"} sort keys %pats)), "\n";
	}

	close($in) if defined $fileName;
	close($out) if defined $outFile;
}

sub distance_matrix
{
	my ($rhash, $rmat) = @_;

	my @pats =  sort keys %$rhash;
	print join("\n", map {"\t$_\t$rhash->{$_}"} @pats), "\n" if $verbose;

	for (my $i = 0; $i < $#pats; $i++)
	{
		for (my $j = $i+1; $j <= $#pats; $j++)
		{
			$rmat->{$pats[$i]}{$pats[$j]}[0] = $rhash->{$pats[$i]} * $rhash->{$pats[$j]}; # number of pairs
			$rmat->{$pats[$i]}{$pats[$j]}[1] = distance($pats[$i], $pats[$j]);
			print "\t\tdistance: $i $j $pats[$i] $pats[$j] @{$rmat->{$pats[$i]}{$pats[$j]}}\n" if $verbose > 1;
		}
	}
}

sub average_distance
{
	my ($cnt, $rmat) = @_;
	return undef if !$cnt;

	my $sum = 0;

	foreach my $p1 (keys %$rmat)
	{
		foreach my $p2 (keys %{$rmat->{$p1}})
		{
			$sum += $rmat->{$p1}{$p2}[0] * $rmat->{$p1}{$p2}[1];
		}
	}

	return $sum / ($cnt * ($cnt-1) / 2);
}

sub add_distance_matrix
{
	my ($rres, $rmet, $racc) = @_;
	my $rhash = keys %$rmet ? $rmet : $racc;

	foreach my $p1 (keys %$rhash)
	{
		$rres->{$p1}{$p1} = 0;
		foreach my $p2 (keys %{$rhash->{$p1}})
		{
			$rres->{$p2}{$p2} = 0;
			$rres->{$p1}{$p2} = ($rmet->{$p1}{$p2}[1]||0) + ($racc->{$p1}{$p2}[1]||0);
			$rres->{$p2}{$p1} = $rres->{$p1}{$p2};
		}
	}
}

sub clustering
{
	my ($rpat, $rmat, $rclu, $rsize) = @_;

	my $pid = open2(my $childOut, my $childIn, $rscript);
	#print "pid: $pid\n";

	my (@ids, @clusters);
	foreach my $id (sort keys %$rpat)
	{
		map {push(@ids, $id)} 1..$rpat->{$id};
	}

	print $childIn join("\t", @ids), "\n";
	for (my $i = 0; $i < @ids; $i++)
	{
		print $childIn join("\t", map {$ids[$i] eq $ids[$_] ? 0 : $ids[$i] lt $ids[$_] ? $rmat->{$ids[$i]}{$ids[$_]}[1] : $rmat->{$ids[$_]}{$ids[$i]}[1]} 0..$#ids), "\n";
	}

	close($childIn);

	while (<$childOut>)
	{
		print "cut: $_" if $verbose;
		@clusters = split(/\s+/);
	}

	close($childOut);
	waitpid( $pid, 0 );
	my $code = $? >> 8;
	print "child: $code\n" if $verbose;

	map {$rclu->{$ids[$_]} = $clusters[$_]; $rsize->[$clusters[$_]]++} 0..$#clusters;
}

sub distance
{
	my ($p1, $p2) = @_;

	my @s1 = split(//, $p1);
	my @s2 = split(//, $p2);
	my ($cnt, $sum) = (0, 0);

	for (my $i = 0; $i < @s1; $i++)
	{
		next if $s1[$i] eq '-' || $s2[$i] eq '-';
		$sum += abs($s1[$i] - $s2[$i]);
		$cnt++;
	}

	return $cnt ? $sum / $cnt : 0;
}

sub print_distance_matrix
{
	my ($out, $rmat) = @_;
	my @ids = sort keys %$rmat;
	#print $out "ID: @ids\n";

	foreach my $p1 (@ids)
	{
		print $out join("\t", $p1, map {!exists $rmat->{$p1}{$_} ? '-' : ref($rmat->{$p1}{$_}) eq 'ARRAY' ? "(@{$rmat->{$p1}{$_}})" : $rmat->{$p1}{$_}} @ids), "\n";
	}
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

	if ($helpFlag)
	{
		die("Arguments: [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}

	$verbose = 0 if !$verbose;
}
