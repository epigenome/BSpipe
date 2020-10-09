#!/usr/bin/perl -w
#!/bin/perl -w

## Author: JeongHyeon Choi
## Date: Jul 2 2004

## in silico conversion

## Usage:
## CHECK usage with option  -h

## $Id: fasta2seq.pl,v 1.6 2007-06-05 13:30:54 jeochoi Exp $


use Getopt::Long;

my ($helpFlag, $inFile, $outFile, $includeFile, $excludeFile, @types, $mode);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"i|infile:s"	=> \$inFile,		## input file
	"o|outfile:s"	=> \$outFile,		## output file
	"t|type:s"		=> \@types,			## type of name
	"mode=s"			=> \$mode,
);

checkOptions();

my $in  = openInput($inFile);
my $out = openOutput($outFile);

for (my $i = 0; <$in>; $i++)
{
	if ($i%4 == 0)
	{
		die "Wrong format: $i $_" if !/^@/;
	}
	elsif ($i%4 == 1)
	{
		if (uc $mode eq 'C')
		{
			y/Tt/Cc/;
		}
		elsif (uc $mode eq 'G')
		{
			y/Aa/Gg/;
		}
	}
	elsif ($i%4 == 2)
	{
		die "Wrong format: $i $_" if !/^\+/;
	}

	print $out $_;
}

close($in) if !$inFile;
close($out) if !$outFile;

#-------------------------------------------------------------------------------
sub openInput
{
	my ($fileName) = @_;

	return STDIN unless defined $fileName;

	my ($fd);
	open($fd, $fileName=~/\.gz$/ ? "zcat $fileName|" : $fileName) || die "Open error:$fileName";
	return $fd;
}

sub openOutput
{
	my ($fileName) = @_;

	return STDOUT unless defined $fileName;

	my ($fd);
	open($fd, ">$fileName") || die "Open error:$fileName";
	return $fd;
}

sub checkOptions
{
	$inFile  = shift(@ARGV) if !defined $inFile  && @ARGV > 0;
	$outFile = shift(@ARGV) if !defined $outFile && @ARGV > 0;
	die "The output file name is the same as the input file name\n" if defined $inFile && defined $outFile && $inFile eq $outFile;

	if ($helpFlag || !$mode)
	{
		die "Arguments: -m mode [[-i] input_file] [[-o] output_file]\n"
		. "\tmode should be either c, g, or cg\n";
	}
}

