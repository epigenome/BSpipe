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

my ($helpFlag, $inFile, $outPre, $verbose, $quiet, $header, $colStr, $field, $high, $low, $cutoff);

GetOptions(
	"h|?|help"    => \$helpFlag,	
	"input=s"     => \$inFile,		## input file
	"output=s"    => \$outPre,		## output file
	"verbose+"    => \$verbose,		## verbose output
	"quiet"	     => \$quiet,
	"header"	     => \$header,
	"c|column=s"	  => \$colStr,
	"field=i"	  => \$field,
	"abs"	        => \$absFlag,
	"top=f"	     => \$high,
	"bottom=f"	  => \$low,
	"cutoff=f"	  => \$cutoff,
) || die "\n";

checkOptions();

my (@buf, @head);
my @cols = map {--$_} split(/,/, $colStr);
my $abs = $absFlag ? 'a' : '';
my $undefined = 'NA';

load($inFile);
@buf = sort { value($b->[@cols+0]) <=> value($a->[@cols+0]) } @buf;
$outPre .= "$head[@cols+0]." if $outPre && @head;
process();


#-------------------------------------------------------------------------------

sub value
{
	return $absFlag ? abs($_[0]) : $_[0];
}

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);

	if ($header)
	{
		my @a = split(/[\t\r\n]/, <$in>);
		@head = @a[@cols, $field];
	}

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my @a = split(/[\t\r\n]/);
		push(@buf, [@a[@cols, $field]]) if $a[$field] ne $undefined && (!$cutoff || value($a[$field]) >= $cutoff);
	}

	close($in) if defined $fileName;
}

sub process
{
	print "Processing ", $outPre || 'STDIN', " ...\n" if !$quiet;
	output('top', $high, 1) if $high;
	output('bottom', $low, 0) if $low;
}

sub output
{
	my ($pre, $num, $flag) = @_;

	my $out = openOutput(!$outPre ? *STDOUT : "${outPre}$abs$pre$num");
	print $out join("\t", @head), "\n";
	$num = int($num*@buf + 0.5);

	for (my $i = 0; $i < $num && $i < @buf; $i++)
	{
		print $out "$pre\t" if !$outPre;
		print $out join("\t", @{$buf[$flag ? $i : $#buf-$i]}), "\n";
	}

	close($out) if $outPre;
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
	$inFile  = shift(@ARGV) if !defined $inFile  && @ARGV > 0;
	$outPre = shift(@ARGV) if !defined $outPre && @ARGV > 0;

	if ($helpFlag || !$field || (!$high && !$low))
	{
		die("Arguments: -f field [-t cutoff_top] [-b cutoff_bottom] [-abs] [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}

	$field--;
	if ($outPre)
	{ 
		if (-d $outPre) { $outPre .= '/' if $outPre !~ /\/$/; }
		else            { $outPre .= '.'; }
	}
}
