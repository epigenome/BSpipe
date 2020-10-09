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
			"a"=>"t",
			"t"=>"a",
			"c"=>"g",
			"g"=>"c",
			"m"=>"k",
			"r"=>"y",
			"y"=>"r",
			"w"=>"w",
			"s"=>"s",
			"k"=>"m",
			"b"=>"v",
			"d"=>"h",
			"h"=>"d",
			"v"=>"b",
			"n"=>"n"
);
			
use Getopt::Long qw(:config no_ignore_case);

my ($helpFlag, $inFile, $outPrefix, $verbose, $quiet, $minRead, $track, $name, $gzip);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outPrefix,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"min=i"			=> \$minRead,
	"name=s"			=> \$name,
	"track=s"		=> \$track,
	"gzip"			=> \$gzip,
) || die "\n";

checkOptions();

my (@old, @cur, %outs, %tar);
my $trackLine = '';

load($inFile);


#-------------------------------------------------------------------------------

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);

	if ($name)
	{
		$trackLine = "track name=$name itemRgb=On";
		$trackLine .= $track if $track;
		$trackLine .= "\n";
	}

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my @a = split(/[\t\r\n]/);

		die "Wrong format in $_" if $a[8] !~ /mread=(\d+);uread=(\d+)/;
		@cur = ($a[0], $a[3]-1, $a[4], $1+$2, $a[5], $a[6], $1, $2, lc($a[2]));
		onLineForPalindrome();

		if (@old)
		{
			if ($old[0] eq $cur[0] && $old[1] == $cur[1]) #CG
			{
				$cur[6] += $old[6];
				$cur[7] += $old[7];
				$cur[5] = '.';
				$cur[3] = $cur[6]+$cur[7];
				$cur[4] = sprintf("%.3f", $cur[3] == 0 ? -1 : $cur[6]/$cur[3]);
			}
			else
			{
				onLine('', @old);
			}
		}

		@old = @cur;
	}

	onLine('', @old) if @old;
	close($in) if defined $fileName;
	map {close($outs{$_}) || die "Error in closing $_\n"} keys %outs;
}

sub onLineForPalindrome
{
	if (!exists $tar{$cur[8]})
	{
		die "Line $.: $cur[8] doesn't have C\n" if $cur[8] !~ /^([^c]*)c/;
		$tar{$cur[8]} = [join('', map {$RC{$_}} reverse(split(//, $cur[8]))), length($1)];
	}

	if ($cur[8] eq $tar{$cur[8]}[0])
	{
		my $off = $tar{$cur[8]}[1];
		#onLine("s$cur[8]", $cur[0], $cur[5] eq '-' ? ($cur[2]-$off-1, $cur[2]-$off) : ($cur[1]+$off, $cur[1]+$off+1), @cur[3..$#cur]);
		onLine("s$cur[8]", @cur);
	}
}

sub onLine
{
	my ($type, @a) = @_;
	$type = $a[8] if !$type;
	my $out;

	if (exists $outs{$type})
	{
		$out = $outs{$type};
	}
	else
	{
		$out = $outs{$type} = openOutput("$outPrefix.$type.bed");
		print $out $trackLine;
	}

	if ($a[3] >= $minRead)
	{
		print $out join("\t", @a[0..$#a-2]), "\n";
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
	open($fd, $fileName =~ /.gz$/ ? "| gzip -c > $fileName" : $gzip ? "| gzip -c > $fileName.gz" : $fileName =~ /.bz(ip)?2$/ ? "| bzip2 -z -c > $fileName" : ">$fileName") || die("Open error: $fileName");
	return $fd;
}

sub checkOptions
{
	$inFile  = shift(@ARGV) if !defined $inFile  && @ARGV > 0;
	$outPrefix = shift(@ARGV) if !defined $outPrefix && @ARGV > 0;

	if ($helpFlag || !$outPrefix)
	{
		die("Arguments: [-n name] [-t track] [[-i] in_file] [-o] out_prefix [-v] [-q]\n"
		  . "\t\n"
		  );
	}

	$minRead = 0 if !$minRead;
}
