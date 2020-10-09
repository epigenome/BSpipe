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
use File::Basename;

my ($helpFlag, $seqFile, $outPre, $verbose, $quiet);
my ($class, $depth, $metScore, $metCount, $umtScore, $umtCount, $difCount, $vValue);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"sequence=s"	=> \$seqFile,
	"output=s"		=> \$outPre,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"class=s"		=> \$class,
	"d|depth=i"		=> \$depth,
	"ms=f"			=> \$metScore,
	"mc=i"			=> \$metCount,
	"us=f"			=> \$umtScore,
	"uc=i"			=> \$umtCount,
	"dc=i"			=> \$difCount,
	"v=f"				=> \$vValue,
) || die "\n";

checkOptions();
my @tfiles;

my $outDir = -d $outPre ? $outPre : dirname($outPre);
-e $outDir || mkdir($outDir) || die "Can't create directory: $outDir\n";

foreach my $file (@ARGV)
{
	die "Not found: $file\n" if !-e $file;
	my $bfile = basename($file);
	my $ofile = "$outDir/$bfile.$$";
	push(@tfiles, $ofile);
	run("split.pl -i $file -o $ofile");
	run("java -cp $class DataFormat3Normalization $depth $ofile");
	merge("$ofile/*.norm", "${outPre}$bfile.norm");
	run("java -cp $class Identify_URsandMRs $umtScore $metScore $umtCount $metCount $ofile");
	merge("$ofile/*.hot", "${outPre}$bfile.hot");
}

my $dmrDir = "$outDir/dmr.$$";
mkdir $dmrDir if !-e $dmrDir;

if (@ARGV == 2)
{
	run("java -cp $class TwoSamples_DMRs $difCount $dmrDir @tfiles");
}
elsif (@ARGV > 2)
{
	run("java -cp $class MutipleSamples_DMRs $difCount $vValue $dmrDir @tfiles");
	merge("$dmrDir/*.spe", "${outPre}spe");
}

merge("$dmrDir/*.dmr", "${outPre}dmr");
merge("$dmrDir/*.con", "${outPre}con");

run("fasta_split.pl -i $seqFile -o $outDir/seq.$$");
run("java -cp $class ConservedAnalysis $outDir/seq.$$ $dmrDir $dmrDir");
merge("$dmrDir/*.cseq", "${outPre}cseq");
merge("$dmrDir/*.cgcoe", "${outPre}cgcoe");

run("java -cp $class DMRsandSample_specificAnalysis $outDir/seq.$$ $dmrDir $dmrDir");
merge("$dmrDir/*.dseq", "${outPre}dseq");
merge("$dmrDir/*.dgcoe", "${outPre}dgcoe");

END {system("rm -r $outDir/*.$$");}

#-------------------------------------------------------------------------------

sub run
{
	print "@_\n" if $verbose;
	!system("@_") || die "Error in @_\n";
}

sub merge
{
	my ($src, $dst) = @_;
	run("cat $src | perl -ne 'print if \$i++ == 0 || !/Chromosome/' > $dst");
}

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
	if ($helpFlag || !$class || !$seqFile || !$outPre || @ARGV < 2)
	{
		die("Arguments: -c class_dir [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}

	$depth = 5 if !defined $depth;
	$metScore = 0.7 if !defined $metScore;
	$metCount = 3 if !defined $metCount;
	$umtScore = 0.3 if !defined $umtScore;
	$umtCount = 3 if !defined $umtCount;
	$difCount = 4 if !defined $difCount;
	$vValue = 0.4 if !defined $vValue;
	if (-d $outPre)
	{
		$outPre .= '/' if $outPre !~ /\/$/;
	}
	else
	{
		$outPre .= '.';
	}
}
