#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Author: Jeong-Hyeon Choi
## Date: May 28 2008

## bisulphite a sequence under a given percent

## Usage:
## CHECK usage with option  -h

## $Id: perl_template,v 1.2 2007-09-21 22:51:23 jeochoi Exp $


use Getopt::Long qw(:config no_ignore_case);

my ($helpFlag, $inFile, $outPrefix, $verbose, $quiet);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outPrefix,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
) || die "\n";

checkOptions();


load($inFile);


#-------------------------------------------------------------------------------

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	my $cout = openOutput($outPrefix ? "$outPrefix.c.fa" : $outPrefix);
	my $gout = openOutput($outPrefix ? "$outPrefix.g.fa" : $outPrefix);

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;

		if (/^\s*>/)
		{
			print $cout $_;
			print $gout $_;
		}
		else
		{
			my $temp = $_; $temp =~ tr/Tt/Cc/;
			print $cout $temp;
			$temp = $_; $temp =~ tr/Aa/Gg/;
			print $gout $temp;
		}
	}
}

#-------------------------------------------------------------------------------

sub openInput
{
	my ($fileName) = @_;

	return *STDIN unless defined $fileName;

	my ($fd);
	open($fd, $fileName =~ /.gz$/ ? "zcat $fileName |" : $fileName) || die("Open error: $fileName");
	return $fd;
}

sub openOutput
{
	my ($fileName) = @_;

	return *STDOUT unless defined $fileName;

	my ($fd);
	open($fd, $fileName =~ /.gz$/ ? "| zcat -c > $fileName" : ">$fileName") || die("Open error: $fileName");
	return $fd;
}

sub checkOptions
{
	$inFile  = shift(@ARGV) if !defined $inFile  && @ARGV > 0;
	$outPrefix = shift(@ARGV) if !defined $outPrefix && @ARGV > 0;

	if ($helpFlag)
	{
		die("Arguments: [[-i] in_file] [[-o] out_prefix] [-v] [-q]\n"
		  . "\t\n"
		  );
	}
}
