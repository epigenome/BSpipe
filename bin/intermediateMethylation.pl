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

my ($helpFlag, $inFile, $outFile, $verbose, $quiet, $first, $second, $third);

GetOptions(
	"h|?|help"    => \$helpFlag,	
	"input=s"     => \$inFile,		## input file
	"output=s"    => \$outFile,		## output file
	"verbose+"    => \$verbose,		## verbose output
	"quiet"       => \$quiet,
	"1|first=s"   => \$first,
	"2|second=s"  => \$second,
	"3|third=s"   => \$third,
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

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my @a = split;
		$a[3] = '-' if $a[3] eq '0';
		$a[4] = '-' if $a[4] eq '0';
		my $state = $a[3] ne '-' ? $a[4] ne '-' ? $second : $first : $a[4] ne '-' ? $third : $second;
		print $out join("\t", @a[0..2], join(":", $state, @a[3,4])), "\n";
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
	open($fd, $fileName =~ /.gz(ip)?$/ ? "| gzip -c > $fileName" : $fileName =~ /.bz(ip)?2$/ ? "| bzip2 -z -c > $fileName" : ">$fileName") || die("Open error: $fileName");
	return $fd;
}

sub checkOptions
{
	$inFile  = shift(@ARGV) if !defined $inFile  && @ARGV > 0;
	$outFile = shift(@ARGV) if !defined $outFile && @ARGV > 0;
	die "The output file name is the same as the input file name\n" if defined $inFile && defined $outFile && $inFile eq $outFile;

	if ($helpFlag || !$first || !$second || !$third)
	{
		die("Arguments: -1 first_state -2 second_state -3 third_state [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}
}
