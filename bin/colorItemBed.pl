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

my ($helpFlag, $inFile, $outFile, $verbose, $quiet, $track, $score);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"track=s"		=> \$track,
	"score=i"		=> \$score,
) || die "\n";

checkOptions();

load($inFile);


#-------------------------------------------------------------------------------

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	my $out = openOutput($outFile);
	my $trackLine = '';

	print $out "track $track itemRgb=On\n" if ($track);

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my @a = split(/[\t\r\n]/);
		print $out join("\t", @a[0,1,2], sprintf("%.2f", $a[$score]), $score != 4 && $a[3] eq '.' ? $a[$score] : 0, $a[5]||'.', @a[1,2], color($a[$score])), "\n";
	}

	close($in) if defined $fileName;
}

sub color
{
	my ($m) = @_;
	$m = 0 if $m < 0;
	$m = 1 if $m > 1;
	my $red = int($m*255);
	my $green = 255-$red;
	return "$red,$green,0";
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
		die("Arguments: [-track string][[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}

	$score = 5 if !defined $score;
	$score--;
}
