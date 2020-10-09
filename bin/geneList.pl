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

my ($helpFlag, $inFile, $outPrefix, $verbose, $quiet, $iScoreStr, $sign, @cols);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outPrefix,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"s|score=s"		=> \$iScoreStr,
	"sign=s"			=> \$sign,
	"col=i"			=> \@cols,
) || die "\n";

checkOptions();

my %hash;
my @iScore = $iScoreStr ? _parseRangeArg1($iScoreStr) : ();

load($inFile);
process();


#-------------------------------------------------------------------------------

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my @a = split(/[\t\r\n]/);
		my $score = 0;
		map {$score += $a[$_]} @iScore;
		next if $a[-4] =~ /intergen/ || ($sign && ($sign eq '+' ? $score < 0 : $score));

		foreach my $b (@a[@cols])
		{
			foreach my $c (split(/\|/, $b))
			{
				my ($g, $fs) = $c =~ /(.+)\((.+)\)/;
				next if !$fs || $fs =~ /far/;

				my %ts;
				map { s/\d+$//; $ts{lc $_}++ } split(/,/, $fs);

				foreach (keys %ts)
				{
					push(@{$hash{$_}{$g}}, $score);
					push(@{$hash{'exon'}{$g}}, $score) if /cds|utr/;
					push(@{$hash{'body'}{$g}}, $score) if /intron|exon|cds|utr/;
					push(@{$hash{'all'}{$g}}, $score);
				}
			}
		}
	}

	close($in) if defined $fileName;
}

sub process
{
	foreach my $key (keys %hash)
	{
		print "$key\n";
		my $out = openOutput("$outPrefix.$key");
		print $out join("\t", '#Gene', 'Count');
		print $out "\t", join("\t", 'Min', 'Avg', 'Max') if $iScoreStr;
		print $out "\n";

		foreach my $gene (sort keys %{$hash{$key}})
		{
			print $out join("\t", $gene, summary($hash{$key}{$gene})), "\n";
		}

		close($out);
	}
}

sub summary
{
	my ($ra) = @_;
	my $num = @$ra;
	return $num if !$iScoreStr;

	my ($sum, $min, $max) = (0);
	map {$sum += $_; $min = $_ if !defined $min || $min > $_; $max = $_ if !defined $max || $max < $_;} @$ra;
	return ($num, $min, sprintf("%.3f", $sum/$num), $max);
}

sub _decrease
{
	my ($array, $value) = @_;
	map {$_-=$value} @$array;
}

sub _parseRangeArg1
{
	my @array = _parseRangeArg(@_);
	_decrease(\@array, 1);
	return @array;
}

sub _parseRangeArg
{
	my @ret;

	foreach my $arg (@_)
	{
		foreach my $r (split(/,/, $arg))
		{
			push(@ret, ($r =~ /(\d+)-(\d+)/) ? $1..$2 : $r);
		}
	}

	return @ret;
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

	$fileName =~ s/'//g;
	my ($fd);
	open($fd, $fileName =~ /.gz$/ ? "| gzip -c > $fileName" : $fileName =~ /.bz(ip)?2$/ ? "| bzip2 -z -c > $fileName" : ">$fileName") || die("Open error: $fileName");
	return $fd;
}

sub checkOptions
{
	$inFile  = shift(@ARGV) if !defined $inFile  && @ARGV > 0;
	$outPrefix = shift(@ARGV) if !defined $outPrefix && @ARGV > 0;

	if ($helpFlag || !$outPrefix)
	{
		die("Arguments: [[-i] in_file] [-o] out_prefix [-v] [-q]\n"
		  . "\t\n"
		  );
	}

	if (!@cols) { @cols = (-2, -1); }
	else        { @cols = map {--$_} @cols; }
}
