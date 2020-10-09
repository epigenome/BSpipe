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

my ($helpFlag, $siteFile, $inFile, $outFile, $verbose, $quiet, $exclude, $header, $minCov);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"filter=s"		=> \$siteFile,		## input file
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"exclude!"		=> \$exclude,
	"header!"		=> \$header,
	"cov=i"			=> \$minCov,
) || die "\n";

checkOptions();

my %sites;
load($siteFile) if $siteFile;
process($inFile);


#-------------------------------------------------------------------------------

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my @a = split(/[\t\r\n]/);
		$sites{$a[0]}{$a[1]}++;
	}

	close($in) if defined $fileName;
}

sub process
{
	my ($fileName) = @_;

	print "Processing ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	my $out = openOutput($outFile);

	while (<$in>)
	{
		next if /^\s*$/;
		if (/track|browse/ || ($header && /^#/))
		{
			print $out $_;
		}
		else
		{
			my @a = split(/[\t\r\n]/);
			print $out $_ if (!$siteFile || ($exclude ^ exists $sites{$a[0]}{$a[1]})) && (!defined $minCov || $a[3] >= $minCov);
		}
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

	if ($helpFlag || (!$siteFile && !$inFile))
	{
		die("Arguments: [-ex] -f filter_bed [-c min_cov] [-i] bed_file [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}

	$exclude = 0 if !$exclude;
}
