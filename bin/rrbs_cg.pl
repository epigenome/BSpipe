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
my ($enzyme);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"enzyme=s"		=> \$enzyme,
) || die "\n";

checkOptions();

my ($out, $digestLen);
my ($oid, $oseq, $oref, $obeg, $oend);
my ($cid, $cseq, $cref, $cbeg, $cend);

preprocess();
load($inFile);


#-------------------------------------------------------------------------------

sub preprocess
{
	print STDERR "$enzyme ";
	my @a = split(/-/, $enzyme);
	die "No cutting site: $enzyme\n" if @a < 2;
	my $cut = length($a[0]);
	$common = length($a[1])-$cut;
	my $site = substr($a[1], 0, $common);
	$digestLen = rindex($site, 'CG');
	print STDERR "=> common=$common, lastCG=$digestLen\n";
}

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	$out = openOutput($outFile);

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		($cid, $cseq) = split;
		($cref, $cbeg, $cend) = $cid =~ /(.+):(\d+)-(\d+)/;
		die "Wrong format in $cid\n" if !$cend;
		onFragment() if $oid;
		($oid, $oseq, $oref, $obeg, $oend) = ($cid, $cseq, $cref, $cbeg, $cend);
	}
		onFragment() if $oid;

	close($in) if defined $fileName;
}

sub onFragment
{
	my $frag = $oend-$obeg+1+$common;

	if ($oseq =~ /^CG/i)
	{
		print $out join("\t", $oref, $obeg-1, $obeg+1), "\n";
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

	if ($helpFlag || !$enzyme)
	{
		die("Arguments: -e enzyme [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}
}
