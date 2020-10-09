#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Author: Jeong-Hyeon Choi
## Date: 2015

## Description

## Usage:
## CHECK usage with option  -h

## $Id: perl_template,v 1.2 2007-09-21 22:51:23 jeochoi Exp $


use Getopt::Long qw(:config no_ignore_case);
use File::Basename;

my ($helpFlag, $inFile1, $inFile2, $outFile, $verbose, $quiet, $refFile, $path, $cov, $dist, $diff, $tmpDir, $leave, $name);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"i|input1=s"	=> \$inFile1,		## input file
	"j|input2=s"	=> \$inFile2,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"reference=s"	=> \$refFile,
	"name=s"	      => \$name,
	"path:s"			=> \$path,
	"tmp=s"			=> \$tmpDir,
	"coverage=i"	=> \$cov,
	"distance=i"	=> \$dist,
	"difference=i"	=> \$diff,
	"leave"	      => \$leave,
) || die "\n";

checkOptions();

my @rmFiles;
my $tmp = $tmpDir || (exists $ENV{TMPDIR} ? $ENV{TMPDIR} : "/tmp/methclone.$$");
(-d $tmp) || mkdir($tmp) || die "Can't create temp directory: $tmp\n";

$inFile1 && (check($inFile1) || convert(\$inFile1));
$inFile2 && (check($inFile2) || convert(\$inFile2));

process() if $inFile2 && $outFile;

unlink(@rmFiles) if !$leave && !$tmpDir;

#-------------------------------------------------------------------------------

sub check
{
	my ($inFile) = @_;
	my $in = openInput("samtools view $inFile |");
	my $flag;

	for (my $i = 0; $i < 1000 && ($_ = <$in>); $i++)
	{
		if (/\bXM:Z:/i)
		{
			$flag = 1;
			last;
		}
	}

	close($in);
	return $flag;
}

sub convert
{
	my ($rinFile) = @_;
	my $tmpFile = ($tmpDir || ($leave ? dirname($outFile) : $tmp)) . '/' . basename($$rinFile);
	die "Use -tmp or don't specify -leave\n" if $$rinFile eq $tmpFile;
	my $cmd = "java -Xmx256m -jar " . dirname($0) . "/methcloneEncode.jar $refFile $$rinFile $tmpFile";
	run($cmd);
	push(@rmFiles, $tmpFile, "$tmpFile.bai");
	$$rinFile = $tmpFile;
	$tmpFile =~ s/\.bam/.bai/;
	if (-e $tmpFile) { rename($tmpFile, "$$rinFile.bai"); }
	else             { run("samtools index $$rinFile"); }
}

sub process
{
	my $cmd = $path || 'methclone';
	$cmd .= " $inFile1 $inFile2 $outFile $name";
	$cov  && ($cmd .= " $cov");
	$dist && ($cmd .= " $dist");
	$diff && ($cmd .= " $diff");

	run($cmd);
}

sub run
{
	my ($cmd) = @_;
	print "$cmd\n" if $verbose;

	system($cmd) && die  "Error in $$: $cmd\n";
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
	$name = 'sample' if !$name;

	if ($helpFlag || !$inFile1 || ($inFile2 && !$outFile) || (!$inFile2 && $outFile))
	{
		die("Arguments: -i bam_file [-j bam_file -o out_file] [-n name] [-c coverage] [-dist distance] [-diff difference] [-v] [-q]\n"
		  . "\t\n"
		  );
	}
}
