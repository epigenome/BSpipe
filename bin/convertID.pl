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

my ($helpFlag, $tabFile, $inFile, $outFile, $verbose, $quiet, @fields, @values, $column);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"table=s"		=> \$tabFile,		## input file
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"f|field=i"		=> \@fields,
	"v|value=i"		=> \@values,
	"column=i"		=> \$column,
) || die "\n";

checkOptions();

my %hash;

load($tabFile);
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
		for (my $i = 0; $i < @fields; $i++)
		{
			$a[$fields[$i]-1] =~ s/\([^\)]+\)//g;
			$a[$values[$i]-1] =~ s/\([^\)]+\)//g;
			my @b = split(/[\|,]/, $a[$fields[$i]-1]);
			my @c = split(/[\|,]/, $a[$values[$i]-1]);
			foreach my $k (@b)
			{
				map {s/\.\d+$//; $hash{$k}{$_}++} @c;
			}
		}
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
		next if /^#/ || /^\s*$/;
		my @a = split(/[\t\r\n]/);
		next if !exists $hash{$a[$column-1]};
		print $out (keys %{$hash{$a[$column-1]}})[0], "\n";
	}

	close($in) if defined $fileName;
	close($out) if defined $outFile;
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

	if ($helpFlag || (!$tabFile && !$inFile) || !@fields || !@values || @fields != @values || !$column)
	{
		die("Arguments: -t map_file [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}
}
