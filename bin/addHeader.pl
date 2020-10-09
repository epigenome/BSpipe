#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Author: Jeong-Hyeon Choi
## Date: 2009

## Description

## Usage:
## CHECK usage with option  -h

## $Id: perl_template,v 1.2 2007-09-21 22:51:23 jeochoi Exp $


use Getopt::Long qw(:config no_ignore_case);

my ($helpFlag, $inFile, $outFile, $verbose, $quiet, $header);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"header=s"		=> \$header,
) || die "\n";

checkOptions();

load($inFile);


#-------------------------------------------------------------------------------

sub load
{
	my ($fileName) = @_;

	my $in = openInput($fileName);
	my $out = openOutput($outFile);

	if ($header)
	{
		$header =~ s/\\t/\t/g;
		print $out $header, "\n";
	}

	while (<$in>)
	{
		print $out $_;
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

	if ($helpFlag)
	{
		die("Arguments: [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}
}
