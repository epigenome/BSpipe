#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Author: Jeong-Hyeon Choi
## Date: Oct 14 2008

## join two or more files

## Usage:
## CHECK usage with option  -h

## $Id: perl_template,v 1.2 2007-09-21 22:51:23 jeochoi Exp $


use Getopt::Long qw(:config no_ignore_case);
use File::Basename;

my ($helpFlag, $faiFile, @inFiles, $outFile, $verbose, $quiet, $fieldStr, $keyStr, $valueStr, $missing, $common, $number, $headerFlag, $line);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"fai=s"			=> \$faiFile,		## input file
	"input=s"		=> \@inFiles,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"f|field=s"		=> \$fieldStr,
	"key=s"			=> \$keyStr,
	"v|value=s"		=> \$valueStr,
	"missing=s"		=> \$missing,
	"common"			=> \$common,
	"number"			=> \$number,
	"header"			=> \$headerFlag,
	"line=i"			=> \$line,
) || die "\n";

checkOptions();

my @sizes;
my @headers;
my $del = '!#:';
my @fields = _parseRangeArg1($fieldStr);
my @keys = split(/,/, $keyStr) if $keyStr;
my @values = _parseRangeArg1($valueStr) if $valueStr;
my (@ins, $out);
my %order;
loadIndexFile($faiFile);

foreach my $i (0..$#inFiles)
{
	print "Loading ", $inFiles[$i] || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($inFiles[$i]);
	$ins[$i] = $in;
}

my @buf;
my $flag = 1;
my @min = 0..$#inFiles;


while (@min)
{
	foreach my $i (@min)
	{
		$buf[$i] = load($i);
	}

	if ($flag)
	{
		outputHeader();
		$flag = 0;
	}

	@min = ();
	outputData();
}

map {close($_)} @ins;
close($out) if $out;

#-------------------------------------------------------------------------------

sub load
{
	my ($i) = @_;
	my @a;
	my $in = $ins[$i];

	while (<$in>)
	{
		next if (!$headerFlag && /^#/) || /^\s*$/;
		chop; s/$/a/;
		@a = split(/[\t\r\n]/);
		$a[-1] =~ s/a$//;

		if (!defined $sizes[$i])
		{
			$sizes[$i] = @a-@fields;
			if ($headerFlag)
			{
				$headers[0] = join("\t", map {s/"//g; $_} @a[@fields]) if $i == 0;
				$headers[$i+1] = join("\t", map {s/"//g; $_} @a[@values]);
			}
		}
		elsif (@values && $sizes[$i] != @a-@fields)
		{
			die "Inconsistent #fields in $_";
		}

		last;
	}

	return \@a;
}

sub outputHeader
{
	$out = openOutput($outFile);

	if ($headerFlag)
	{
		print join("\t", @headers), "\n";
		print $out join("\t", @headers), "\n";
	}
	else
	{
		print $out '#', join("\t", $keyStr ? @keys : (map {"Key".$_} 1..@fields));
		my @files = _getDistinctName(\@inFiles, '[\. ]');
		for (my $i = 0; $i < @files; $i++) { map {print $out "\t$files[$i]"} 1..@values+0; }
		print $out "\n";
	}
}

sub outputData
{
	foreach my $i (0..$#inFiles)
	{
		next if !@{$buf[$i]};
		die "$buf[$i][0] is not in $faiFile\n" if !exists $order{$buf[$i][0]};

		if (@min == 0 || $order{$buf[$i][0]} < $order{$buf[$min[0]][0]})
		{
			@min = ($i);
		}
		elsif ($order{$buf[$i][0]} == $order{$buf[$min[0]][0]})
		{
			if    ($buf[$i][1] <  $buf[$min[0]][1]) { @min = ($i); }
			elsif ($buf[$i][1] == $buf[$min[0]][1]) { push(@min, $i); }
		}
	}

	return if !@min || ($common && @min != @inFiles);

	print $out join("\t", @{$buf[$min[0]]}[@fields]);
	my $j = 0;
	foreach my $i (0..$#inFiles)
	{
		print $out "\t", join("\t", $j < @min && $i == $min[$j] ? @{$buf[$min[$j++]]}[@values] : map {$missing} 0..$#values);
	}
	print $out "\n";
}

sub loadIndexFile
{
	my ($fileName) = @_;
	my $in = openInput($fileName);

	while (<$in>)
	{
		next if (/^#/) || /^\s*$/;
		my @a = split(/[\t\r\n]/);
		$order{$a[0]} = keys %order;
	}

	close($in);
}

#-------------------------------------------------------------------------------

sub _getDistinctName
{
	my ($rfiles, $del) = @_;
	$del = '[\._]' if !$del;

	my @arrays;
	map { my @a = split(/$del/, fileparse($_)); push(@arrays, \@a) } @$rfiles;

	my %strs;
	for (my $i = 0; $i < @arrays; $i++)
	{
		map { $strs{$_}{$i}++ } @{$arrays[$i]};
	}

	foreach my $s (keys %strs)
	{
		next if keys(%{$strs{$s}}) < @arrays;
		for (my $i = 0; $i < @arrays; $i++)
		{
			for (my $j = 0; $j < @{$arrays[$i]}; $j++)
			{
#				if ($arrays[$i][$j] eq $s) { $arrays[$i][$j] = ''; last; }
				if ($arrays[$i][$j] eq $s) { $arrays[$i][$j] = ''; }
			}
		}
	}

	return map {join('', @$_)} @arrays;
}

sub _parseRangeArg1
{
	my @array = _parseRangeArg(@_);
	_decrease(\@array, 1);
	return @array;
}

sub _decrease
{
	my ($array, $value) = @_;
	map {$_-=$value} @$array;
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
	push(@inFiles, @ARGV) if !@inFiles  && @ARGV > 0;

	if ($helpFlag || !$fieldStr || !$valueStr || !@inFiles || !$faiFile)
	{
		die("Arguments: -f key [-header] [-m missing_value] [-v value] -fai ref_fai [-i] in_file .. [-o out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}

	$missing = '' if !defined $missing;
}
