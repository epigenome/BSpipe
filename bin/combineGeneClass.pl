#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Author: Jeong-Hyeon Choi
## Date: May 13 2008

## Description

## Usage:
## CHECK usage with option  -h

## $Id: perl_template,v 1.2 2007-09-21 22:51:23 jeochoi Exp $

my $Null = '.';
my ($out5, $end5, $body, $end3, $out3, $outside, $none) = ("farup", "upstream", "body", "downstream", "fardown", "intergenic", "none");
my ($exon, $intron, $utr5, $cds, $utr3) = ('exon', 'intron', "5UTR", 'CDS', "3UTR");
my ($biend5, $biend3) = ("bi_$end5", "bi_$end3");
my @combinedTypes = ($biend5, $end5, $biend3, $end3, $utr5, $utr3, $cds, $exon, $intron, $outside);
my %encode = ('U' => $end5, '5' => $utr5, 'E' => $exon, 'I' => $intron, '3' => $utr3, 'D' => $end3);

use Getopt::Long qw(:config no_ignore_case);

my ($helpFlag, $inFile, $outFile, $verbose, $quiet, $field1, $field2, $scoreFlag, $binSize, $priority);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"f|field1=i"	=> \$field1,
	"g|field2=i"	=> \$field2,
	"score!"			=> \$scoreFlag,
	"bin=f"			=> \$binSize,
	"priority=s"   => \$priority,
) || die "\n";

checkOptions();

my (%class, %score);
my (%allclass, %comclass);
my (%subclass, %subfinal);
my $total = 0;

my %pr = parse($priority);
process($inFile);


#-------------------------------------------------------------------------------

sub getOneType
{
	my ($str) = @_;
	$str = 'none' if !$str;
	my %types;

	foreach my $type (split(/,+/, $str))
	{
		$types{$type}++;
	}

#	return $oneKeyFt->(\%types);
	my $t = (sortPriority(keys %types))[0];
	return $t eq $out5 || $t eq $out3 ? $outside : $t
}

sub getFinalType
{
	my ($a, $b) = @_;

	if ($a eq $b)
	{
		return ($a eq $end5 ? $biend5 : $a eq $end3 ? $biend3 : $a, $a eq $outside ? '' : '*');
	}

	my $t = (sortPriority($a, $b))[0];
	return $t ne $a ? ($b, '-') : ($a, '+');
}

sub getAllType
{
	my ($a, $b, $rtypes, $rdirs) = @_;
	my (%types, %dirs);
	map {if ($_ ne $out5 && $_ ne $out3) {$types{$_}++; $dirs{"$_+"}++;}} split(/,+/, $a);
	map {if ($_ ne $out5 && $_ ne $out3) {$types{$_}++; $dirs{"$_-"}++;}} split(/,+/, $b);
	if (!keys %types) { $types{$outside}++; $dirs{$outside}++; }

	my $comtype = join(',', sortPriority(keys %types));
	$comtype =~ s/,*body//;

	if (exists $types{$end5} && $types{$end5} > 1) { delete $types{$end5}; $types{$biend5}++; }
	if (exists $types{$end3} && $types{$end3} > 1) { delete $types{$end5}; $types{$biend3}++; }

	@$rtypes = keys %types;
	@$rdirs  = keys %dirs ;
	return $comtype;
}

sub sortPriority
{
	map { die "Undefined category in the priority: $_\n" if !exists $pr{$_} } @_;
	return sort {$pr{$a} <=> $pr{$b}} @_;
}

sub process
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	my $out = openOutput($outFile);
	my $sum = openOutput("$outFile.stat");
	my $flag = 1;

	while (<$in>)
	{
		next if /^\s*$/;
		my @a = split(/[\t\r\n]/);
		my (%types1, %types2);

		if ($flag)
		{
			my $h = $a[$field2]; $h =~ s/\(.*\)$/(Both)/;
			print $out join("\t", @a, $h, $h, 'Dist(Both)'), "\n";
			$flag = 0;
		}
		elsif (!/^#/)
		{
			my $dist1 = $a[$field1-1] =~ /^(-*\d+)\|/ ? $1 : $Null;
			my $dist2 = $a[$field2-1] =~ /^(-*\d+)\|/ ? $1 : $Null;
			my $dist = $dist1 eq $Null ? $dist2 : $dist2 eq $Null || abs($dist1) < abs($dist2) ? $dist1 : $dist2;
			$a[$field1-1] = $Null if $a[$field1-1] !~ s/^(-*\d+)\|//;
			$a[$field2-1] = $Null if $a[$field2-1] !~ s/^(-*\d+)\|//;

			my $type1 = getOneType($a[$field1]);
			my $type2 = getOneType($a[$field2]);
			my ($ctype, $dir) = getFinalType($type1, $type2);
			print "$type1 $type2 $ctype $dir\n" if $verbose;
			$class{$ctype}++;
			$score{int($a[4]/$binSize)}{$ctype}++ if $scoreFlag;

			my (@keys, @dirs);
			my $comtype = getAllType($a[$field1], $a[$field2], \@keys, \@dirs);
			$comclass{$comtype}++;
			for (my $i = 0; $i < @keys; $i++) { $allclass{$keys[$i]}++; }

			print $out join("\t", @a, $comtype, "$ctype$dir", $dist), "\n";
			$total++;
		}
		else
		{
			s/^#//;
			print $sum $_;
		}
	}

	close($in) if defined $fileName;

	foreach my $type (@combinedTypes)
	{
		print $sum join("\t", "all", $type, $allclass{$type} || 0), "\n";
	}

	foreach my $type (@combinedTypes)
	{
		print $sum join("\t", "priority", $type, $class{$type} || 0), "\n";
	}

	foreach my $type (sort {$comclass{$b} <=> $comclass{$a}} keys %comclass)
	{
		print $sum join("\t", "combination", $type, $comclass{$type} || 0), "\n";
	}

	if ($scoreFlag)
	{
		my @ss = sort {$a <=> $b} keys %score;
			print $sum join("\t", 'score', 'type', map {$_*$binSize} @ss), "\n";

		foreach my $type (@combinedTypes)
		{
			print $sum join("\t", 'score', $type);
			foreach my $s (@ss)
			{
				print $sum "\t", $score{$s}{$type} || 0;
			}
			print $sum "\n";
		}
	}

	print $sum "total\t$total\n";
}

sub parse
{
	my ($n, %prio) = (0);

	foreach my $c (split(//, $_[0]))
	{
		die "Unsupported character for categories: $c\n" if !exists $encode{$c};
		$prio{$cds       } = $n++ if $c eq 'E';
		$prio{$encode{$c}} = $n++;
		$prio{$body      } = $n++ if $c eq 'I';
	}

	$prio{$out5   } = $n++;
	$prio{$out3   } = $n++;
	$prio{$outside} = $n++;
	$prio{$none   } = $n++;
	return %prio;
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
	open($fd, $fileName =~ /.gz$/ ? "| zcat -c > $fileName" : ">$fileName") || die("Open error: $fileName");
	return $fd;
}

sub checkOptions
{
	$inFile  = shift(@ARGV) if !defined $inFile  && @ARGV > 0;
	$outFile = shift(@ARGV) if !defined $outFile && @ARGV > 0;
	die "The output file name is the same as the input file name\n" if defined $inFile && defined $outFile && $inFile eq $outFile;

	if ($helpFlag || !$field1 || !$field2)
	{
		die("Arguments: -f field_1 -g field_2 [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}

	$field1--; $field2--;
	$binSize = 1 if !defined $binSize;
}
