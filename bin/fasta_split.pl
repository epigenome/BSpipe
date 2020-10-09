#!/usr/bin/perl -w
#!/bin/perl -w

## Author: JeongHyeon Choi
## Date: Jul 2 2004

## Will split flat of Pfam to separate by each family

## Usage:
## CHECK usage with option  -h

## $Id: fasta_split.pl,v 1.6 2006-12-11 16:28:21 jeochoi Exp $


use Getopt::Long qw(:config no_ignore_case);

my ($helpFlag, $inFile, $outDir, $verbose, $type, $suffix);
my ($exFile, %exclusions);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"i|inFile:s"	=> \$inFile,		## input file
	"o|outDir=s"	=> \$outDir,		## input file
	"t|type:s" 		=> \$type,
	"s|suffix:s" 	=> \$suffix,
	"v|verbose" 	=> \$verbose,
	"exclude:s"		=> \%exclusions,
	"Exclude:s"		=> \$exFile,
);

checkOptions();

mkdir($outDir) if !-e $outDir;
die "Unexist directory:$outDir\n" if ! -e $outDir;

load($exFile) if $exFile;

my $no = 0;
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
	my $fd;

	while (<$in>)
	{
		# >P15711 104K_THEPA 104 kDa microneme-rhoptry antigen.
		if (/^>/)
		{
			my $id = extractId($_);
			my $name  = "$outDir/$id";
			   $name .= ".$suffix" if (defined $suffix);
			print "$name: $_" if $verbose;
			if ($fd) { close($fd); $fd = undef; }
			open($fd, ">$name") || die("Open error:$name") if !exists $exclusions{$id};

			$no++;
		}

		print $fd $_ if defined $fd;
	}

	close($in);
}

sub extractId
{
	my ($buf) = @_;

	return $no if ($type eq 'no');

	die "Unknown id:\n" . substr($buf, 0, 60) . "...\n" if ($buf !~ /^>\s*(\S+)/);
	return $1 if ($type eq 'id');

	if ($type eq 'last')
	{
		die "Unknown id:\n" . substr($buf, 0, 60) . "...\n" if ($buf !~ /\|([^\|\s]+)\s/);
		return $1;
	}

	die "Unknown id:\n" . substr($buf, 0, 60) . "...\n" if ($buf !~ /$type\|([^|]+)/);
	return $1;
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

sub checkOptions
{
	$inFile = shift(@ARGV) if !defined $inFile  && @ARGV > 0;
	$outDir = shift(@ARGV) if !defined $outDir && @ARGV > 0;

	if ($helpFlag || !defined $outDir)
	{
		die("Arguments: [-i fasta_file] -o out_dir [-t id|no|gi|ref|last] [-s suffix] [-v]\n"
		  );
	}

	$type = 'id' if !defined $type;
}

