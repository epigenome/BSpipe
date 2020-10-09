#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Author: Jeong-Hyeon Choi
## Date: Nov 26 2008

## convert or add hugo IDs

## Usage:
## CHECK usage with option  -h

## $Id: perl_template,v 1.2 2007-09-21 22:51:23 jeochoi Exp $


use Getopt::Long qw(:config no_ignore_case);

my $Null = '.';
my ($helpFlag, $hugoFile, $inFile, $outFile, $verbose, $quiet, $colFrom, $colTo, $colInput, $addFlag, $remain);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"hugo=s"			=> \$hugoFile,
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"field=i"		=> \$colFrom,
	"value=s"		=> \$colTo,
	"column=s"		=> \$colInput,
	"add"				=> \$addFlag,
	"remain"			=> \$remain,
) || die "\n";

checkOptions();

my %hugos;
my @fields = _parseRangeArg1($colInput);
my @values = _parseRangeArg1($colTo);

load($hugoFile);
process($inFile);


#-------------------------------------------------------------------------------

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

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my @a = split(/[\t\r\n]/);
		if ($a[$colFrom])
		{
			foreach my $from (split(/[\s,]+/, $a[$colFrom]))
			{
				foreach (@values)
				{
					next if (!$a[$_]);
					$hugos{$from} = $a[$_];
					last;
				}
			}
		}
	}

	close($in) if defined $fileName;
}

sub process
{
	my ($fileName) = @_;

	print "Processing ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	my $out = openOutput($outFile);

	while (<$in>)
	{
		next if /^\s*$/;
		my @a = split(/[\t\r\n]/);

		if (/^#/)
		{
			push(@a, @a[@fields]) if $addFlag;
			print $out join("\t", @a), "\n";
		}
		else
		{
			foreach my $f (@fields)
			{
				my $id;

				if ($a[$f] && $a[$f] ne $Null)
				{
					my @b = split(/[ \|]/, $a[$f]);
					my %ids;
					foreach my $i (@b)
					{
						next if $i =~ /^-*\d+$/;
						my ($k, $fs) = $i =~ /^([^\(]+)\((.*)\)/;
						if ($k)
						{
							die "Wrong format in $i\n" if !$fs;
							map {$ids{$hugos{$k}}{$_}++} split(/,/, $fs) if exists $hugos{$k};
						}
						else
						{
							$ids{$hugos{$i}}++ if exists $hugos{$i};
						}
					}
					$id = join('|', map {$_ . (%{$ids{$_}} ? ("(" . join(',', sort keys %{$ids{$_}}) . ")") : '')} sort keys %ids) if keys %ids;
				}
		
				$id = $remain ? $a[$f] : $Null if !$id;

				if ($addFlag) { push(@a, $id); }
				else          { $a[$f] = $id; }
			}
			print $out join("\t", @a), "\n";
		}
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
	open($fd, $fileName =~ /.gz$/ ? "| gzip -c > $fileName" : $fileName =~ /.bz(ip)?2$/ ? "| bzip2 -z -c > $fileName" : ">$fileName") || die("Open error: $fileName");
	return $fd;
}

sub checkOptions
{
	$inFile  = shift(@ARGV) if !defined $inFile  && @ARGV > 0;
	$outFile = shift(@ARGV) if !defined $outFile && @ARGV > 0;
	die "The output file name is the same as the input file name\n" if defined $inFile && defined $outFile && $inFile eq $outFile;

	if ($helpFlag || (!$hugoFile && !$inFile) || !defined $colFrom || !defined $colInput || !defined $colTo)
	{
		die("Arguments: -hugo hugo_file -f col_in_hugo_file -c col_in_input_file [-a] [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}

	die "Wrong column index for hugo file: $colFrom\n" if $colFrom < 1; $colFrom--;
}
