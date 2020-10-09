#!/usr/bin/perl -w
#!/bin/perl -w

## Author: JeongHyeon Choi
## Date: Jul 2 2004

## Will make one sequence from fasta file

## Usage:
## CHECK usage with option  -h

## $Id: fasta2seq.pl,v 1.6 2007-06-05 13:30:54 jeochoi Exp $


use Getopt::Long;

my ($helpFlag, $inFile, $outFile, $comment, $headerFlag, $fullFlag, @types);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"i|infile:s"	=> \$inFile,		## input file
	"o|outfile:s"	=> \$outFile,		## output file
	"c|comment!"	=> \$comment,
	"header!"		=> \$headerFlag,
	"full!"			=> \$fullFlag,
	"t|type:s"		=> \@types,			## type of name
);

checkOptions();
$comment = 1 unless defined $comment;

my $in  = openInput($inFile);
my $out = openOutput($outFile);
my $first = 1;

while (<$in>)
{
	chop;

	if (/^>\s*(\S+)/)
	{
		print $out "\n" if !$first;
		$id = $fullFlag ? $' : $1;
		print $out $id, "\t" if $headerFlag;
		$first = 0;
	}
	elsif (!$comment || ! /^(#|\/\/)/)
	{
#		s/\s+//g;
		print $out $_;
	}
}

print $out "\n";
close($in);


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
		die "Arguments: [[-i] input_file] [[-o] output_file] [-header|-noheader] [-c|-noc]\n";
	}
}

