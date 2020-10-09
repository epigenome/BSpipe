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

my ($helpFlag, $inDir, $outFile, $grpFile, $verbose, $quiet, $suffix, $fullFlag);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \$inDir,		## input file
	"group=s"		=> \$grpFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"suffix=s"		=> \$suffix,
	"full"			=> \$fullFlag,
) || die "\n";

checkOptions();

my %sam2grp;
load($grpFile) if $grpFile;

my $count = 0;
my $out = openOutput($outFile);

foreach my $file (<$inDir/*>)
{
	my $bfile = basename($file);
	if ($bfile !~ /^(.+?)(_R*([12]))*\.($suffix)$/)
	{
		warn "Unexpected format: $bfile\n" if $count++ < 5;
		next;
	}

	my ($type, $sam) = ($3||1, $1);
	print $out join("\t", $fullFlag ? $file : $bfile, $type, $sam, group($sam)), "\n";
}

close($out);
warn "$count files were ignored\n" if $count;

#-------------------------------------------------------------------------------

sub group
{
	my ($sam) = @_;

	foreach my $pat (keys %sam2grp)
	{
		return $sam2grp{$pat} if ($sam =~ /$pat/);
	}

	return $sam;
}

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my @a = split;
		die "Duplicate sample: [$a[0]] to [$sam2grp{$a[0]}] [$a[1]]\n" if exists $sam2grp{$a[0]};
		$sam2grp{$a[0]} = $a[1];
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
	$inDir  = shift(@ARGV) if !defined $inDir  && @ARGV > 0;
	$outFile = shift(@ARGV) if !defined $outFile && @ARGV > 0;
	die "The output file name is the same as the input file name\n" if defined $inDir && defined $outFile && $inDir eq $outFile;

	if ($helpFlag || !$inDir)
	{
		die("Arguments: [-full] [-g sam2grp_file] [-i] in_dir [-s suffix] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}

	$suffix = "fastq|fq|fastq.gz|fq.gz" if !$suffix;
}
