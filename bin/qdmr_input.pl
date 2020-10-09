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

my ($helpFlag, $inFile, $outFile, $verbose, $quiet);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
) || die "\n";

checkOptions();

my (@head, @data);

load($inFile);
process($outFile);


#-------------------------------------------------------------------------------
=format
#Key1 Key2  Key3  MU5   MU9   MU13  MU15  MU45  MU48  MCG18 MCG14 MCG8  MCG16 MCG10 MCGN
chr19 19539 19541 NA NA NA NA NA NA NA NA NA NA 1.000 NA

#1.2
777762   11
ID Chromosome*Start*End NPC_1 NBR   BC1063T  BC1133T  BC1142T  BC1063S  BC1133S  BC1142S  mayo22   mayo39   mayo59
chr1:496 chr1*496*498   94 95 96 55 97 100   50 94 100   97 89
=cut

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	@head = split(/[\t\r\n]/, <$in>);

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		push(@data, $_) if !/\bNA\b/;
	}

	close($in) if defined $fileName;
}

sub process
{
	my ($fileName) = @_;

	print "Processing ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $out = openOutput($outFile);
	print $out "#1.2\n";
	print $out @data+0, "\t", @head-3, "\n";
	print $out join("\t", 'ID', 'Chromosome*Start*End', @head[3..$#head]), "\n";

	foreach (@data)
	{
		my @a = split(/[\t\r\n]/);
		print $out join("\t", "$a[0]:$a[1]", "$a[0]*$a[1]*$a[2]", @a[3..$#head]), "\n";
	}

	close($out) if defined $fileName;
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
