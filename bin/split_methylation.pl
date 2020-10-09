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

my ($helpFlag, $inFile, $outDir, $verbose, $quiet, $context);

GetOptions(
	"h|?|help"    => \$helpFlag,	
	"input=s"     => \$inFile,		## input file
	"output=s"    => \$outDir,		## output file
	"verbose+"    => \$verbose,		## verbose output
	"quiet"       => \$quiet,
	"context=s"   => \$context,
) || die "\n";

checkOptions();

if (!$context)
{
	die "Can't get cytosine context from non-specified input file\n" if !$inFile;
	($context) = $inFile =~ /\.([^\.\\]+)\.methy\./;
	die "Can't get cytosine context from $inFile\n" if !$context;
}

load($inFile);


#-------------------------------------------------------------------------------

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	my @header = split(/[\t\r\n]/, <$in>);
	my %sams;
	map {$sams{$_}++} @header;
	my $start = 3;
	$start++ if $header[$start] =~ /^Key|ID/;

	foreach (keys %sams)
	{
		die "Duplicate sample name: $_\n" if $sams{$_} > 1;
	}

	my %outs;
	map {$outs{$_} = openOutput("$outDir/$header[$_].$context.bed.gz")} $start..$#header;

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my @a = split(/[\t\r\n]/);
		map {my $o = $outs{$_}; print $o join("\t", @a[0..2], '.', $a[$_]), "\n" if $a[$_] ne 'NA'} $start..$#header;
	}

	close($in) if defined $fileName;
	map {close($outs{$_})} keys %outs;
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
	open($fd, $fileName =~ /.gz(ip)?$/ ? "| gzip -c > $fileName" : $fileName =~ /.bz(ip)?2$/ ? "| bzip2 -z -c > $fileName" : ">$fileName") || die("Open error: $fileName");
	return $fd;
}

sub checkOptions
{
	$inFile  = shift(@ARGV) if !defined $inFile  && @ARGV > 0;
	$outDir = shift(@ARGV) if !defined $outDir && @ARGV > 0;
	die "The output file name is the same as the input file name\n" if defined $inFile && defined $outDir && $inFile eq $outDir;

	if ($helpFlag || !$outDir)
	{
		die("Arguments: [-c context] [[-i] in_file] [-o] out_dir [-v] [-q]\n"
		  . "\t\n"
		  );
	}
}
