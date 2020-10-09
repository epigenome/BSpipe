#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Author: Jeong-Hyeon Choi
## Date: Apr 16 2011

## pick the best gene for each feature from the output of geneClass.pl

## Usage:
## CHECK usage with option  -h

## $Id: perl_template,v 1.2 2007-09-21 22:51:23 jeochoi Exp $

my $Null = '.';
my ($out5, $end5, $body, $end3, $out3, $outside, $none) = ("farup", "upstream", "body", "downstream", "fardown", "intergenic", "none");
my ($exon, $intron, $utr5, $cds, $utr3) = ('exon', 'intron', "5UTR", 'CDS', "3UTR");
my @colNames = ('Ref', 'Start', 'End', 'Strand', 'Name', $out5, $end5, $utr5, $cds, $intron, $utr3, $end3, $out3, 'Alias');
my ($iRef, $iStart, $iEnd, $iStrand, $iName, $i5out, $i5end, $i5utr, $iCds, $iIntron, $i3utr, $i3end, $i3out, $iAlias) = 0..20;
my %encode = ('U' => $end5, '5' => $utr5, 'E' => $exon, 'I' => $intron, '3' => $utr3, 'D' => $end3);

use Getopt::Long qw(:config no_ignore_case);

my ($helpFlag, $inFile, $outFile, $sumFile, $verbose, $quiet, $cntFile, $priority);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"summary=s"		=> \$sumFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"count=s"		=> \$cntFile,
	"priority=s"   => \$priority,
) || die "\n";

checkOptions();

my %feats;
my %pr = parse($priority);

load($inFile);
process($outFile);


#-------------------------------------------------------------------------------

=format
#ref  start end   ID 5'end 5'UTR CDS   intron   3'UTR 3'end
#chr1  852952   854817   -  NR_026874   856558-856559:497.36:856150-856861:85:26.43:0.30   -  -  -  -  852213-852214:1632.35:851806-852582:190:54.33:0.31  18|852935-852936:509.04:852593-853403:91:26.43:0.31   FLJ39609
=cut 
sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	my $cnt = openOutput($cntFile) if $cntFile;
	my $sum = openOutput($sumFile) if $sumFile;

	while (<$in>)
	{
		if (/^#/ || /^\s*$/)
		{
			print $cnt $_ if $cntFile;
			s/($out5|$out3)\t//g;
			print $sum $_ if $sumFile;
			next;
		}

		my @a = split(/[\t\r\n]/);
		my @c = @a;

		for (my $i = $i5out; $i <= $i3out; $i++)
		{
			if ($a[$i] eq $Null)
			{
				$c[$i] = 0;
				next;
			}

			my @b = split(/\|/, $a[$i]);
			$c[$i] = @b;

			foreach my $c (@b)
			{
				my @d = split(/:/, $c, 4);
				my ($start, $end) = split(/\-/, $d[0]);
				die "Wrong format in $c\n$_" if !defined $end;
				my $dist = int(($start + $end) / 2) - ($a[$iStrand] eq '-' ? $a[$iEnd] : $a[$iStart]); # TSS

				if (  !exists $feats{$a[$iRef]} || !exists $feats{$a[$iRef]}{$start}
					||  $pr{$colNames[$feats{$a[$iRef]}{$start}[0]]} >  $pr{$colNames[$i]}
					|| ($pr{$colNames[$feats{$a[$iRef]}{$start}[0]]} == $pr{$colNames[$i]} && abs($feats{$a[$iRef]}{$start}[10]) > abs($dist)) )
				{
					$feats{$a[$iRef]}{$start} = [$i, $end, $d[1], $d[2], $d[3], @a[$iStart, $iEnd, $iStrand, $iName, $iAlias], $dist];
				}
			}
		}

		print $cnt join("\t", @c), "\n" if $cntFile;
		print $sum join("\t", @a[$iRef..$iName, $i5end..$i3end, $iAlias]), "\n" if $sumFile && ($a[$i5end] ne $Null || $a[$i5utr] ne $Null || $a[$iCds] ne $Null || $a[$iIntron] ne $Null || $a[$i3utr] ne $Null || $a[$i3end] ne $Null);
	}

	close($in) if defined $fileName;
	close($cnt) if $cntFile;
	close($sum) if $sumFile;
}

sub process
{
	my ($fileName) = @_;

	print "Processing ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $out = openOutput($outFile);
	print $out "#Ref\tStart\tEnd\tInfo\tScore\tStrand\tDistance\tCategory\tName\tSymbol\tStart\tEnd\tStrand\n";

	foreach my $ref (sort
		{
			my ($aa) = $a =~ /(\d+)/ ? $1 : 100000000;
			my ($bb) = $a =~ /(\d+)/ ? $1 : 100000000;
			$aa <=> $bb || $a cmp $b
		}
		keys %feats)
	{
		foreach my $start (sort {$a <=> $b} keys %{$feats{$ref}})
		{
			my $cat = @colNames[$feats{$ref}{$start}[0]];
			$cat = $outside if $cat eq $out5 || $cat eq $out3;
			print $out join("\t", $ref, $start,
				$feats{$ref}{$start}[1],
				$feats{$ref}{$start}[2] || $Null,
				$feats{$ref}{$start}[3] || $Null,
				$feats{$ref}{$start}[4] || $Null,
				$feats{$ref}{$start}[7] eq '-' ? -$feats{$ref}{$start}[10] : $feats{$ref}{$start}[10],
				$cat,
				$feats{$ref}{$start}[8],
				$feats{$ref}{$start}[9],
				$feats{$ref}{$start}[5],
				$feats{$ref}{$start}[6],
				$feats{$ref}{$start}[7]), "\n";
		}
	}

	close($out) if defined $fileName;
}

sub parse
{
	my ($n, %prio) = (0);

	foreach my $c (split(//, $_[0]))
	{
		die "Unsupported character for categories: $c\n" if !exists $encode{$c};
		$prio{$body      } = $n   if $c =~ /[EI53]/;
		$prio{$cds       } = $n++ if $c eq 'E';
		$prio{$encode{$c}} = $n++;
	}

	$prio{$out5} = $n++;
	$prio{$out3} = $n++;
	$prio{$none} = $n++;
	return %prio;
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

	if ($helpFlag || !$priority)
	{
		die("Arguments: -p priority [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\tprority is a combination of U(upstream), 5(5'UTR), E(exon), I(intron), 3(3'UTR), and D(downstream)\n"
		  );
	}
}
