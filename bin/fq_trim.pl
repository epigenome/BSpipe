#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Author: Jeong-Hyeon Choi
## Date: Jul 31 2008

## trim low quality reads and bases

## Usage:
## CHECK usage with option  -h

## $Id: perl_template,v 1.2 2007-09-21 22:51:23 jeochoi Exp $


use Getopt::Long qw(:config no_ignore_case);

my %types = ( 'PHRED33' => 33, 'PHRED64' => 64, 'SOLEXA64' => 59 );
my ($helpFlag, $inFile, $outFile, $logFile, $verbose, $quiet, $type, $trim5, $trim3, $n5, $n3, $qual5, $qual3, $minLength);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"type=s"			=> \$type,
	"5q=i"			=> \$qual5,
	"3q=i"			=> \$qual3,
	"5l=i"			=> \$trim5,
	"3l=i"			=> \$trim3,
	"5n!"				=> \$n5,
	"3n!"				=> \$n3,
	"length=f"		=> \$minLength,
	"log=s"			=> \$logFile,
) || die "\n";

checkOptions();

my ($out, $log, $id, @buf, %stat);
my ($total, $trim, $clean) = (0, 0, 0);

load($inFile);
summary();

#-------------------------------------------------------------------------------

sub summary
{
	my $log = openOutput($logFile);

	foreach my $len (sort {$a<=>$b} keys %stat)
	{
		print $log "#reads of length $len: $stat{$len}\n";
	}

	print $log "Total   reads: $total\n";
	print $log "trimmed reads: $trim\n";
	print $log "clean   reads: $clean\n";
}

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	$out = openOutput($outFile);

	while (<$in>)
	{
		chop;

		if (!@buf && /^@\s*(\S+)/)
		{
			$id = $1;
		}

		push(@buf, $_);
		onSeq() if $id && @buf == 4;
	}

	close($in) if defined $fileName;
	close($out) if defined $outFile;
}

sub onSeq
{
	my $len = my $org = length($buf[1]);
	my ($left, $right) = (0, $len);

	$trim5 && ($left   = $trim5);
	$trim3 && ($right -= $trim3);

	if ($n3)
	{
		$buf[1] =~ /(N*)$/i;
		$right > $len-length($1) && ($right = $len-length($1));
	}

	if ($n5)
	{
		$buf[1] =~ /^(N*)/i;
		$left < length($1) && ($left = length($1));
	}

	if ($qual5)
	{
		$left++  while (ord(substr($buf[3], $left -1, 1)) < $qual5 && $left < $len-1);
	}

	if ($qual3)
	{
		$right-- while (ord(substr($buf[3], $right-1, 1)) < $qual3 && $right > 0)
	}

	$len = $right-$left;
	$len = 0 if $len < 0;
	if ($minLength > 1 ? $len >= $minLength : $len/$org >= $minLength)
	{
		print $out join("\n", $buf[0], substr($buf[1], $left, $len), $buf[2], substr($buf[3], $left, $len)), "\n";
		$clean++;
	}

	$total++;
	$stat{$len}++;
	$trim++ if $len != $org;

	($id, @buf) = ();
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

	if ($helpFlag || (!$type && ($qual5 || $qual3)) || (!$qual5 && !$qual3 && !$trim5 && !$trim3 && !$minLength))
	{
		die("Arguments: [-type PHRED33|PHRED64|SOLEXA64] [-5q quality] [-3q quality] [-5l trim] [-3l trim] [-l min_length] [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}

	if ($type)
	{
		die "Unknown quality score type\n" if !exists $types{$type};
		$qual5 += $types{$type} if $qual5;
		$qual3 += $types{$type} if $qual3;
	}

	$minLength = 1 if !$minLength;
}
