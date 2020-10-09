#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Author: Jeong-Hyeon Choi
## Date: 2008

## Description

## Usage:
## CHECK usage with option  -h

## $Id: perl_template,v 1.2 2007-09-21 22:51:23 jeochoi Exp $


use Getopt::Long qw(:config no_ignore_case);

my ($helpFlag, $inFile, $outFile, $verbose, $quiet);
my ($enzyme);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"enzyme=s"		=> \$enzyme,
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
) || die "\n";

checkOptions();

my ($cut, $out, $id);
my ($seq, $line) = ('', '');

preprocess();
load($inFile);


#-------------------------------------------------------------------------------

sub preprocess
{
#	print "$enzyme ";
	my @a = split(/-/, $enzyme);
	die "No cutting site: $enzyme\n" if @a < 2;
	$cut = length($a[0]);
	$enzyme = join('', @a);
#	print "=> $enzyme\n";
}

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	$out = openOutput($outFile);

	while ($line = <$in>)
	{
		next if $line =~ /^#/ || $line =~ /^\s*$/;
		chop($line);

		if ($line =~ /^>\s*(\S+)/)
		{
			my $temp = getId($1);
			check() if $id;
			($id, $seq) = ($temp, '');
		}
		else
		{
			$seq .= $line;
		}
	}

	if ($id) { check(); }
	close($in) if defined $fileName;
}

sub check
{
	while ($seq =~ /$enzyme/ig)
	{
		output(pos($seq)+$cut-length($enzyme));
	}
}

sub output
{
	my ($p) = @_;
	print $out join("\t", $id, $p, $p+2), "\n";
}

sub getId
{
	my ($str) = @_;
	my ($id) = $str =~ /ref\|([^\|]+)\|/;
	return $id || $str;
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
	open($fd, $fileName =~ /.gz$/ ? "| gzip -c > $fileName" : ">$fileName") || die("Open error: $fileName");
	return $fd;
}

sub checkOptions
{
	$inFile  = shift(@ARGV) if !defined $inFile  && @ARGV > 0;
	$outFile = shift(@ARGV) if !defined $outFile && @ARGV > 0;
	die "The output file name is the same as the input file name\n" if defined $inFile && defined $outFile && $inFile eq $outFile;

	if ($helpFlag || !$enzyme)
	{
		die("Arguments: -e restriction [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}
}
