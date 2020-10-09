#!/usr/bin/perl -w
#!/bin/perl -w

## Author: JeongHyeon Choi
## Date: Jul 2 2004

## Will retrun length of sequence

## Usage:
## CHECK usage with option  -h

## $Id: fasta_length.pl,v 1.5 2007-06-05 13:32:39 jeochoi Exp $


use Getopt::Long;

my ($helpFlag, $inFile, $outFile, $headerFlag, $fullFlag, @types);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"i|infile:s"	=> \$inFile,		## input file
	"o|outfile:s"	=> \$outFile,		## output file
	"header!"		=> \$headerFlag,
	"type=s"			=> \@types,
	"full!"			=> \$fullFlag,
);

checkOptions();

my $in  = openInput($inFile);
my $out = openOutput($outFile);
my $len = 0;
my $header;

while (<$in>)
{
	chomp;

	if (!/^>/)
	{
		s/\s+//g;
		$len += length;
	}
	else
	{
		my $tmp = $';
		onSeq($header, $len) if defined $header;
		$header = $tmp;
		$len = 0;
	}
}

onSeq($header, $len) if defined $header;

close($in);


sub onSeq
{
	my ($header, $len) = @_;
	my ($id) = $header =~ /^\s*(\S+)/;

	print $out $fullFlag ? $header : $id, "\t" if ($headerFlag);
	print $out "$len\n";
}

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

	if ($helpFlag)
	{
		die "Arguments: [-i] input_file [-o] output_file [-header [-full]| -noheader]\n";
	}

	$headerFlag = 1 if !defined $headerFlag;
}

