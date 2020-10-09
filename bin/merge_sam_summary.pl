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

my ($helpFlag, $inFile, $outFile, $verbose, $quiet);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
) || die "\n";

checkOptions();

my (%sections, %subsect, %hash);
my %keywords = ('TotalGiven' => -120, 'TotalInFile' => -110, 'Paired' => -100, 'Single' => -90, 'Read1' => -80, 'Read2' => -60, 'TotalMapped' => -50, 'UniquelyMapped' => -40, 'MultiplyMapped' => -20, 'U   MapQual<=0' => -30, 'M   MapQual<=0' => -10, 'PairMapped' => 0);

load($inFile);
output();

#-------------------------------------------------------------------------------

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	my $flag = 0;
	my ($sec);

	while (<$in>)
	{
		next if /^\s*$/;
		if (/^#/)
		{
			if ($flag)
			{
 				die "Unexpected error:\n$_" if !$sec;
				$subsect{$sec} = $_;
			}
			else
			{
				if (!exists $sections{$_})
				{
					$sections{$_} = keys(%sections) + 0;
				}
				$sec = $_;
			}
			$flag = 1;
		}
		else
		{
			my @a = split(/[\t\r\n]/);
			$hash{$sec}{$a[0]}[1] += $a[1];
			$hash{$sec}{$a[0]}[2] += $a[4] if defined $a[4];
			$hash{$sec}{$a[0]}[3] += $a[5] if defined $a[5];
			$flag = 0;
		}
	}

	close($in) if defined $fileName;
}

sub output
{
	my $out = openOutput($outFile);
	my ($total, $mapped, $multi) = (0);

	foreach my $sec (sort {$sections{$a} <=> $sections{$b}} keys %sections)
	{
		print $out $sec;
		print $out $subsect{$sec} if exists $subsect{$sec};

		if ($mapped)
		{
			$mapped = 0;
			map {$mapped += $_->[1]} values %{$hash{$sec}};
		}

		foreach my $k (sort {
			my $aa = exists $keywords{$a} ? $keywords{$a} : $a =~ /(\d+)/ ? $1 : 10000;
			my $bb = exists $keywords{$b} ? $keywords{$b} : $b =~ /(\d+)/ ? $1 : 10000;
			$aa <=> $bb || $a cmp $b
			} keys %{$hash{$sec}})
		{
			my $a = $hash{$sec}{$k};
			$total = $a->[1] if !$total && $k =~ /^(TotalGiven|TotalInFile)$/;
			$mapped = $a->[1] if !$mapped && $k eq 'TotalMapped';
			$multi  = $a->[1] if !$multi && $k eq 'MultiplyMapped';
			$paired  = $a->[1] if !$paired && $k eq 'Paired';

			print $out join("\t", $k, $a->[1]);
			if ($sections{$sec})
			{
				print $out "\t", percent($a->[1], $sec =~ /^#Pair/ ? $paired : $total) if $sec =~ /^#Pair/ ? $paired : $total;
				print $out "\t", percent($a->[1], $sec =~ /Reference/ && $sec =~ /multiply/ ? $multi : $mapped) if $k ne 'TotalMapped' && $mapped && $k ne 'PairMapped';
				print $out "\t$a->[2]\t$a->[3]" if defined $a->[2];
			}
			die "Unexpected error: $sec, @$a\n" if ($a->[2] && !defined $a->[3]) || ($a->[3] && !defined $a->[2]);
			print $out "\n";
		}

		print $out "\n";
	}
}

sub percent
{
	return sprintf("%.2f", 100*$_[0]/$_[1]);
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

	if ($helpFlag)
	{
		die("Arguments: [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}
}
