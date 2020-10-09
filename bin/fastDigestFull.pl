#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Author: Jeong-Hyeon Choi
## Date: 2008

## Description

## Usage:
## CHECK usage with option  -h

## $Id: perl_template,v 1.2 2007-09-21 22:51:23 jeochoi Exp $


use Getopt::Long qw(:config no_ignore_case);

my ($helpFlag, $inFile, $outFile, $verbose, $quiet);
my ($enzyme);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"enzyme=s"		=> \$enzyme,
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
) || die "\n";

checkOptions();

my ($cut, $out, $id);
my ($seq, $line, $pos, $index, $len) = ('', '', 0, 0, 0);

preprocess();
load($inFile);


#-------------------------------------------------------------------------------

sub preprocess
{
#	print STDERR "$enzyme ";
	my @a = split(/-/, $enzyme);
	die "No cutting site: $enzyme\n" if @a < 2;
	$cut = length($a[0]);
	$enzyme = join('', @a);
#	print STDERR "=> $enzyme\n";
}

sub load
{
	my ($fileName) = @_;

	print STDERR "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	$out = openOutput($outFile);

	while ($line = <$in>)
	{
		next if $line =~ /^#/ || $line =~ /^\s*$/;
		chop($line);

		if ($line =~ /^>\s*(\S+)/)
		{
			my $temp = getId($1);
			output(0) if $id;
			($id, $seq, $pos, $index, $len) = ($temp, '', 0, 0, 0);
		}
		else
		{
			check();
		}
	}

	if ($id) { output(0); }
	close($in) if defined $fileName;
}

sub check
{
	$seq .= $line;
	$len += length($line);
	#print "$index $len $seq\n";
	my $new = $pos;
	print "$pos $len $index\n" if $verbose;

	for ( ; $index < $len; )
	{
		my $j = index(uc($seq), uc($enzyme), $cut ? $index : $index+1);
		last if $j == -1;
		my $break = $j+$cut;
		print substr($seq, $break-$cut, length($enzyme)), "\n" if $verbose;
		output($new-$pos, $break);
		$new = $pos+$break;
		$index = $break;
	}

	if ($new != $pos)
	{
		print "N: $new $pos $len => " if $verbose;
		$seq = substr($seq, $new-$pos);
		$len -= $new-$pos;
		$pos = $new;
		print " $pos $len ", substr($seq, 0, length($enzyme)), "\n" if $verbose;
	}

	$index = $len - length($enzyme);
	$index = 0 if $index < 0;
}

sub output
{
	my ($p, $e) = @_;
	my $s = $e ? substr($seq, $p, $e-$p) : substr($seq, $p);
	print "$id $pos $p ", $e||' ', ' ', substr($seq, $e-$cut, length($enzyme)), "\n" if $verbose;
	print $out "$id:", $pos+$p+1, "-", $pos+$p+length($s), "\t$s\n";
}

sub getId
{
	my ($str) = @_;
	my ($id) = $str =~ /ref\|([^\|]+)\|/;
	return $id || $str;
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

	if ($helpFlag || !$enzyme)
	{
		die("Arguments: -e restriction [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}
}
