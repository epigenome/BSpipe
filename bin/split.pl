#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Author: Jeong-Hyeon Choi
## Date: JAN 10 2006

## split a file by the different value for given field

## Usage:
## CHECK usage with option  -h

## $Id: split.pl,v 1.1 2007-08-31 19:43:15 jeochoi Exp $


use Getopt::Long qw(:config no_ignore_case);

my ($helpFlag, $inFile, $outDir, $verbose, $field, $delimeter, $prefix, $suffix, $bufSize, $size);
my ($exFile, %exclusions);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"i|inFile:s"	=> \$inFile,		## input file
	"o|outDir:s"	=> \$outDir,		## output file
	"prefix:s"		=> \$prefix,
	"suffix:s"		=> \$suffix,
	"v|verbose+"	=> \$verbose,		## verbose output
	"f|field:s"		=> \$field,		## to calculate
	"d|delimeter:s"	=> \$delimeter,		## to calculate
	"buf:i"			=> \$bufSize,
	"exclude:s"		=> \%exclusions,
	"Exclude:s"		=> \$exFile,
	"size:s"			=> \$size,
);

checkOptions();

load($exFile) if $exFile;

mkdir($outDir) if !-e $outDir;
die "There are files in $outDir\n" if <$outDir/*>;

my %buf;
my $no = 0;
my %files;
process($inFile);

#-------------------------------------------------------------------------------

sub load
{
	my ($fileName) = @_;

	my $in = openInput($fileName);

	while (<$in>)
	{
		next if /^\s*(#|\/\/)/;
		my @a = split;
		$exclusions{$a[0]}++;
	}
}

sub process
{
	my ($fileName) = @_;

	my $in = openInput($fileName);

	while (<$in>)
	{
		next if /^(#|\/\/)/ || /^\s*$/;

		my @a = $delimeter ? split(/$delimeter/) : split;
		next if exists $exclusions{$a[$field]};

		push(@{$buf{$a[$field]}}, $_);
		die "Too many files: " . keys(%buf)+0 . "\n" if keys(%buf) > 200;

		output($a[$field]) if (@{$buf{$a[$field]}} > $bufSize)
	}

	close($in) if $fileName;

	foreach my $id (keys %buf)
	{
		output($id) if (@{$buf{$id}} > 0);
	}
}

sub output
{
	my ($id) = @_;

	my $flag;
	my $fileName;

	if (exists $files{$id})
	{
		$fileName = $files{$id};
	}
	else
	{
		$fileName = "$outDir/";
		$fileName = "$prefix" if $prefix;
		$fileName .= $size == 1 ? $id : int($no++/$size);
		$fileName .= ".$suffix" if $suffix;
		$files{$id} = $fileName;
	}

	my $out = openInput((-e $fileName ? '>>' : '>') . $fileName);
	foreach my $l (@{$buf{$id}})
	{
		print $out $l;
	}
	@{$buf{$id}} = ();
	close($out);
}

#-------------------------------------------------------------------------------

sub openInput
{
	my ($fileName) = @_;

	return STDIN unless defined $fileName;

	my ($fd);
	open($fd, $fileName =~ /.gz(ip)?$/ ? "zcat $fileName |" : $fileName =~ /.bz(ip)?2$/ ? "bzcat $fileName |" : $fileName) || die("Open error: $fileName");
	return $fd;
}

sub openOutput
{
	my ($fileName) = @_;

	return STDOUT unless defined $fileName;

	my ($fd);
	open($fd, $fileName =~ /.gz$/ ? "| gzip -c > $fileName" : $fileName =~ /.bz(ip)?2$/ ? "| bzip2 -z -c > $fileName" : ">$fileName") || die("Open error: $fileName");
	return $fd;
}

sub checkOptions
{
	$inFile = shift(@ARGV) if !defined $inFile  && @ARGV > 0;
	$outDir = shift(@ARGV) if !defined $outDir && @ARGV > 0;

	if ($helpFlag || !$outDir)
	{
		die("Arguments: [-i in_file] [-o] out_dir [-s suffix] [-f field] [-d delimeter] [-s size] [-e exclude_id ..] [-E exclude_file] [-b buf_line_number] [-v]\n"
		  );
	}

	$field = 1 unless defined $field;
	$field--;
	$bufSize  = 100000 if !$bufSize;
	$size = 1 if !$size;
}
