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
my ($helpFlag, $inFile1, $inFile2, $outFile, $verbose, $quiet);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"i|input1=s"	=> \$inFile1,		## input file
	"j|input2=s"	=> \$inFile2,
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
) || die "\n";

checkOptions();

my (%wins); # ref=>, start=>, 0=file1, (0=positions, 1={pattern=>count}, 2=nreads, 3=npatterns, 4=nsites, 5=methylation, 6=entropy)
            #                 1=file2, (0=positions, 1={pattern=>count}, 2=nreads, 3=npatterns, 4=nsites, 5=methylation, 6=entropy)
my (@refs, $length);

load($inFile1, 0);
load($inFile2, 1);
process();


#-------------------------------------------------------------------------------

=format
#Ref    Start   End     Ncpgs   AP      Nreads  Npats   AE      Positions       Pattern=Count
chr1    13201   13400   4       0.456   19      13      0.894   13308,13320,13338,13353 00--=2,000-=1,0000=2,0011=1,01--=1,0111=1,10--=2,1000=3,11--=1,1101=1,111-=1,1110=2,1111=1
chr1    16201   16400   8       0.370   31      23      0.544   16207,16211,16216,16234,16244,16248,16263,16276 -------0=1,-------1=1,----0000=1,----0111=1,----1111=2,---00000=2,--000000=1,0000000-=3,00000000=4,000001--=2,000010--=1,000111--=1,0110000-=1,101011--=1,110011--=1,1101----=1,11100010=1,111001--=1,1110010-=1,11111---=1,111111--=1,1111111-=1,11111110=1
=cut
sub load
{
	my ($fileName, $index) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	<$in>;

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my @a = split(/[\t\r\n]/);
		$length = $a[2]-$a[1] if !$length;
		push(@refs, $a[0]) if !exists $wins{$a[0]};
		@{$wins{$a[0]}{$a[1]}[$index][0]} = split(/,/, $a[8]); # positions
		map {($pat, $cnt) = split(/=/); $wins{$a[0]}{$a[1]}[$index][1]{$pat} = $cnt} split(/,/, $a[9]); # patterns
		$wins{$a[0]}{$a[1]}[$index][2] = $a[5];
		$wins{$a[0]}{$a[1]}[$index][3] = $a[6];
		$wins{$a[0]}{$a[1]}[$index][4] = $a[3];
		$wins{$a[0]}{$a[1]}[$index][5] = $a[4];
		$wins{$a[0]}{$a[1]}[$index][6] = $a[7];
		#print "Pos=", join(',', @{$wins{$a[0]}{$a[1]}[$index][0]}), "\n";
		#print "Pat=", join(',', keys %{$wins{$a[0]}{$a[1]}[$index][1]}), "\n";
		#print "Nreads=$wins{$a[0]}{$a[1]}[$index][2], Npats=$wins{$a[0]}{$a[1]}[$index][3]\n";
		#print "Nsites=$wins{$a[0]}{$a[1]}[$index][4], Methylation=$wins{$a[0]}{$a[1]}[$index][5], Entropy=$wins{$a[0]}{$a[1]}[$index][5]\n";
	}

	close($in) if defined $fileName;
}

sub process
{
	print "Processing ...\n" if !$quiet;
	$out = openOutput($outFile);
	print $out "#Ref\tStart\tEnd\tNreads\tNsites\tNpats\tCommonSites\tCommonPats\tEntropy\n";

	foreach my $ref (@refs)
	{
		foreach my $pos (sort {$a<=>$b} keys %{$wins{$ref}})
		{
			print "$ref $pos\n" if $verbose;
			onWindow($ref, $pos) if defined $wins{$ref}{$pos}[0] && defined $wins{$ref}{$pos}[1];
		}
	}

	close($out) if defined $outFile;
}

sub onWindow
{
	my ($ref, $start) = @_;
	my (%poss, %common);
	map {$poss{$_}++} @{$wins{$ref}{$start}[0][0]};
	map {$poss{$_}++} @{$wins{$ref}{$start}[1][0]};
	map {$common{$_}++ if $poss{$_} > 1} keys %poss;

	modify_pattern($wins{$ref}{$start}[0][0], \%common, $wins{$ref}{$start}[0][1]) if $wins{$ref}{$start}[0][4] != keys %common;
	modify_pattern($wins{$ref}{$start}[1][0], \%common, $wins{$ref}{$start}[1][1]) if $wins{$ref}{$start}[1][4] != keys %common;

	my ($comPat, $entropy) = relative_entropy(@{$wins{$ref}{$start}[0]}[1,2], @{$wins{$ref}{$start}[1]}[1,2]);

	print $out join("\t", $ref, $start, $start+$length,
		join(',', map {$wins{$ref}{$start}[$_][2]} 0..1),
		join(',', map {$wins{$ref}{$start}[$_][4]} 0..1),
		join(',', map {$wins{$ref}{$start}[$_][3]} 0..1),
		scalar keys %common, $comPat,
		sprintf("%.3f", $entropy)), "\n" if $comPat;
}

sub modify_pattern
{
	my ($rposs, $rcom, $rpats) = @_;
	my @com = sort {$a<=>$b} keys %$rcom;
	my @deletes;
	for (my $i = 0; $i < @$rposs; $i++)
	{
		push(@deletes, $i) if !exists $rcom->{$rposs->[$i]};
	}
	print "\tpositions @$rposs\n" if $verbose > 1;
	print "\tcommon    @com\n" if $verbose > 1;
	print "\tdeletes   @deletes\n" if $verbose > 1;

	my @pats = sort keys %$rpats;
	foreach my $pat (@pats)
	{
		my @s = split(//, $pat);
		map {$s[$_] = ''} @deletes;
		my $ns = join('', @s);
		print "\t$pat -> $ns\n" if $verbose;
		$rpats->{$ns} = $rpats->{$pat} if $ns =~ /\d/;
		delete $rpats->{$pat};
	}
}

sub relative_entropy
{
	my ($ra, $na, $rb, $nb) = @_;
	my ($cnt, $sum) = (0, 0);

	print "\t", join(' ', sort keys %$ra), "\n" if $verbose;
	print "\t", join(' ', sort keys %$rb), "\n" if $verbose;
	foreach my $pat (keys %$ra)
	{
		print "\t$pat is not in Q\n" if !exists $rb->{$pat} && $verbose > 1;
		next if !exists $rb->{$pat};
		$cnt++;
		$sum += $ra->{$pat}/$na * log($ra->{$pat}/$na * $nb/$rb->{$pat});
	}

	return ($cnt, $sum);
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
	$outFile = shift(@ARGV) if !defined $outFile && @ARGV > 0;
	die "The output file name is the same as the input file name\n" if defined $inFile && defined $outFile && $inFile eq $outFile;

	if ($helpFlag || (!$inFile1 && !$inFile2))
	{
		die("Arguments: -i in_file_1 -j in_file_2 [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}

	$verbose = 0 if !$verbose;
}
