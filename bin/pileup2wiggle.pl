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

my ($helpFlag, $inFile, $outFile, $verbose, $quiet, $track);
my ($name, $desc, $color);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"n|name=s"		=> \$name,
	"desc=s"			=> \$desc,
	"color=s"		=> \$color,
	"track!"			=> \$track,
) || die "\n";

checkOptions();

my @buf = ('', 0, 0, 0);
my $out;

load($inFile);


#-------------------------------------------------------------------------------

# chrM  620  T  13   ^$,^$,^$,^$,^$,^$,^$,^$,^$,^$,^$,^$,^$,  BCB?CCAAC@D:D
# seq1  60  T  T  66  0  99  13  ...........^~.^~.   9<<55<;<<<<<<
sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	$out = openOutput($outFile);

	if ($track)
	{
		print $out "track type=wiggle_0 name=\"$name\" autoScale=on visibility=full";
		print $out " description=\"$desc\"" if $desc;
		print $out " color=$color" if $color;
		print $out "\n";
	}

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my @a = split(/[\t\r\n]/);
		my $cov = @a > 6 ? $a[7] : $a[3];

		if ($a[0] ne $buf[0] || $a[1] > $buf[2]+1 || $cov != $buf[3])
		{
			onContig() if $buf[1] != 0;
			@buf = ($a[0], $a[1]-1, $a[1], $cov);
		}

		$buf[2] = $a[1];
	}

	onContig() if $buf[1] != 0;
	close($in) if defined $fileName;
}

sub onContig
{
	print $out join("\t", @buf), "\n";
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
	$outFile = shift(@ARGV) if !defined $outFile && @ARGV > 0;
	die "The output file name is the same as the input file name\n" if defined $inFile && defined $outFile && $inFile eq $outFile;

	if ($helpFlag || ($track && !$name))
	{
		die("Arguments: -n name [-d description] [-c color] [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\tcolor: r,g,b ex) 51,0,255\n"
		  );
	}
}
