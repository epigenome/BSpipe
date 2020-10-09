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

my ($helpFlag, $refFile, $faiFile, $inFile, $outFile, $verbose, $quiet, $baseQual, $mapQual);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"ref=s"			=> \$refFile,
	"fai=s"			=> \$faiFile,
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"bq=i"			=> \$baseQual,
	"mq=i"			=> \$mapQual,
) || die "\n";

checkOptions();

my (%refs, @seq, @cns, @line, $dir, $strand);

load($refFile) if !$faiFile;
process($inFile);
output($outFile);


#-------------------------------------------------------------------------------

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput("fasta2seq.pl -head $fileName |");

	while (<$in>)
	{
		my @a = split;
		$refs{$a[0]} =  uc $a[1];
	}

	close($in) if defined $fileName;
}

sub process
{
	my ($fileName) = @_;
	my $cmd = "samtools view $fileName";

	print "Processing ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput("$cmd |");

	while (<$in>)
	{
		#print;
		next if /^#/ || /^\s*$/;
		die "Wrong format in $_" if !/\bZT:z:(\w)/i;
		$dir = $1 eq 'C' ? +1 : -1;
		@line = split(/\t/);
		$line[1] = int($line[1]);
		next if $line[4] < $mapQual;
		$strand = strand($line[1]);
#		map {$cns[$_]{baseCode($a[2], $a[3]+$_-1, $dir)}[0]{$strand}++ if ord(substr($a[10], $_-1, 1)) >= $baseQual} split(/,/, $1) if /\bYM:z:(\S+)/i;
		map {add($_, 0) if ord(substr($line[10], $_-1, 1)) >= $baseQual} split(/,/, $1) if /\bYM:z:(\S+)/i;
		map {add($_, 1) if ord(substr($line[10], $_-1, 1)) >= $baseQual} split(/,/, $1) if /\bYU:z:(\S+)/i;
	}

	close($in) if defined $fileName;
}

sub add
{
	my ($pos, $flag) = @_;
	my $base = baseCode($line[2], $line[3]+$pos-1, $dir);
	$cns[$pos]{$base}[$flag]{$strand}++;
}

sub strand
{
	return $_[0] & 0x10 ? -1 : 1;
}

sub getSeq
{
	my ($id, $pos, $len) = @_;
	my $rseq;

	if ($faiFile)
	{
		if (!$seq[0] || $seq[0] ne $id)
		{
			$seq[1] = '';
			my $in = openInput("samtools faidx $refFile $id |");
			while (<$in>) { next if /^>/; chop; $seq[1] .= uc $_; }
			close($in);
			$seq[0] = $id;
		}
		$rseq = \$seq[1];
	}
	else
	{
		$rseq = \$refs{$id};
	}

	die "Unknown sequence for $id\n" if !$$rseq;
	return substr($$rseq, $pos, $len);
}

sub baseCode
{
	my ($rid, $pos, $dir) = @_;
	my $base = getSeq($rid, $pos+$dir-1, 1);
	$base =~ tr/ACGT/TGCA/ if $dir == -1;
	#print "($rid, $pos, $dir) = $base\n";
	#$base = 'H' if $base ne 'G';
	$nuc{$base}++;
	return $base;
}

sub output
{
	my ($fileName) = @_;

	print "Processing ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $out = openOutput($outFile);
	my @nuc = sort keys %nuc;
	my @sum;

	foreach my $b (@nuc)
	{
		for (my $i = 1; $i < @cns; $i++)
		{
			my ($mf, $mr) = ($cns[$i]{$b}[0]{1}||0, $cns[$i]{$b}[0]{-1}||0);
			my ($uf, $ur) = ($cns[$i]{$b}[1]{1}||0, $cns[$i]{$b}[1]{-1}||0);
			onPat($out, $b, $i, $mf, $mr, $uf, $ur);

			$sum[$i][0] += $mf;
			$sum[$i][1] += $mr;
			$sum[$i][2] += $uf;
			$sum[$i][3] += $ur;
		}
	}
		for (my $i = 1; $i < @cns; $i++)
		{
			onPat($out, '', $i, $sum[$i][0]||0, $sum[$i][1]||0, $sum[$i][2]||0, $sum[$i][3]||0);
		}
	close($out) if defined $fileName;
}

sub onPat
{
	my ($out, $b, $i, $mf, $mr, $uf, $ur) = @_;
	print $out join("\t", "C$b", 'forward', $i, $mf, $uf, meth($mf, $uf)), "\n";
	print $out join("\t", "C$b", 'reverse', $i, $mr, $ur, meth($mr, $ur)), "\n";
	print $out join("\t", "C$b", 'both'   , $i, $mf+$mr, $uf+$ur, meth($mf+$mr, $uf+$ur)), "\n";
}

sub meth
{
	my ($m, $u) = @_;
	return $m+$u ? $m/($m+$u) : 'NaN';
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

	if ($helpFlag || (!$refFile && !$faiFile) || !$inFile)
	{
		die("Arguments: <-ref ref_file | -fai fai_fiel> [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}

	if (!$refFile)
	{
		$refFile = $faiFile;
		$refFile =~ s/\.fai$//;
		die "Please specify reference sequence file if it is not the prefix of fai file\n" if $refFile eq $faiFile;
	}
}
