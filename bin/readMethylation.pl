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

my ($helpFlag, $inFile, $outFile, $verbose, $quiet, $target_bases);

GetOptions(
	"h|?|help"    => \$helpFlag,	
	"input=s"     => \$inFile,		## input file
	"output=s"    => \$outFile,		## output file
	"verbose+"    => \$verbose,		## verbose output
	"quiet"       => \$quiet,
	'n=s'         => \$target_bases,
) || die "\n";

checkOptions();

my (@bases, $rid, $chr, $start, $end, $dir, %types, $out);
@bases = sort map {uc $_} split(/:/, $target_bases);

load($inFile);


#-------------------------------------------------------------------------------

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	$out = openOutput($outFile);

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my @a = split(/[\t\r\n]/);;

		onRead(@a) if (!$rid || $rid ne $a[1]);
		$start = $a[3] if $start > $a[3];
		$end = $a[4] if $end < $a[4];
		if ($a[7]) { $types{$a[2]}[0]++; }
		else       { $types{$a[2]}[1]++; }
	}

		onRead() if ($rid);
	close($in) if defined $fileName;
}

sub onRead
{
	print $out join("\t", $chr, $rid, $start, $end, $dir, map 
		{
			my ($m, $t) = exists $types{$_} ? ($types{$_}[0]||0, $types{$_}[1]||0) : (0, 0);
			$m = 0 if !$m;
			$t = 0 if !$t;
			$t += $m;
			($m, $t, $t ? sprintf("%.2f", $m/$t) : '-')
		} @bases), "\n" if $rid;

	($chr, $rid, $start, $end, $dir) = @_[0,1,3,4,6];
	%types = ();
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
	$outFile = shift(@ARGV) if !defined $outFile && @ARGV > 0;
	die "The output file name is the same as the input file name\n" if defined $inFile && defined $outFile && $inFile eq $outFile;

	if ($helpFlag)
	{
		die("Arguments: [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}
}
