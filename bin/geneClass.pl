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
#my ($up5, $down5, $up3, $down3) = (1500, 500, 1500, 0);
#my ($up5, $down5, $up3, $down3) = (2000, 500, 500, 500);
my ($up5, $down5, $up3, $down3) = (0, 0, 0, 0);
my ($out5, $end5, $body, $end3, $out3, $outside, $none) = ("farup", "upstream", "body", "downstream", "fardown", "intergenic", "none");
my ($exon, $intron, $utr5, $cds, $utr3) = ('exon', 'intron', "5UTR", 'CDS', "3UTR");
my %encode = ('U' => $end5, '5' => $utr5, 'E' => $exon, 'I' => $intron, '3' => $utr3, 'D' => $end3);

use Getopt::Long qw(:config no_ignore_case);

my ($helpFlag, $inFile, $outFile, $geneFile, $verbose, $quiet, $field, $column, $bedIndex, $priority);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"gene=s"			=> \$geneFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"field=i"		=> \$field,
	"column=i"		=> \$column,
	"bed=i"			=> \$bedIndex,
	"up5=i"        => \$up5,
	"down5=i"      => \$down5,
	"up3=i"        => \$up3,
	"down3=i"      => \$down3,
	"priority=s"   => \$priority,
) || die "\n";

checkOptions();

my (%class, %final);
my (%geneclass, %genefinal);
my (%genes);

my %pr = parse($priority);
load($inFile);
process() if $geneFile;

#-------------------------------------------------------------------------------

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	my $out = openOutput($outFile);
	my $flag = 1;

	while (<$in>)
	{
		next if /^\s*$/;
		my @a = split(/[\t\r\n]/);
		my (%types, %subfeats, %detail);

		if ($flag)
		{
			print $out join("\t", @a[0..$column], $a[$field], $a[$field]), "\n";
			$flag = 0;
		}
		elsif (!/^#/)
		{
			my @b = split(/\|/, $a[$field]);

			if ($a[$field] eq $Null)
			{
				$types{$none}++;
			}
			else
			{
				my ($cbeg, $cend) = @a[1,2];
				for (my $i = 1; $i < @b; $i++)
				{
					my (%gtypes);
					my ($id, $gbeg, $gend, $dir) = $b[$i] =~ /(.+):(\d+)-(\d+)([+\-])/;
					$a[$field] = "-$a[$field]" if $b[0] && $i == 1 && $a[2] < $gbeg;
					my $extra;
					if (!$b[0] && $' =~ /;([^#]+)/)
					{
						map {s/\d+$//; $subfeats{$_}++;} split(/,/, $1);
						$extra = $1
					}
					my $glen = $gend - $gbeg + 1;
					my ($_down5, $_down3) = $glen - $down5 - $down3 > ($down5 >= $down3 ? $down5 : $down3) ? ($down5, $down3) : ($glen/3, $glen/3);

					# always gbeg is less than gend
					if ($dir eq '+')
					{
						if    ($cbeg >=  $gend+$up3)                                                     { $types{$out3}++; $gtypes{$out3}++; }
						elsif ($cend <=  $gbeg-$up5)                                                     { $types{$out5}++; $gtypes{$out5}++; }
						else
						{
							if (isOverlap($cbeg, $cend, $gbeg+$_down5, $gend-$_down3))                    { $types{$body}++; $gtypes{$body}++; }
							if ($up5 != -$down5 && isOverlap($cbeg, $cend, $gbeg-$up5,    $gbeg+$_down5)) { $types{$end5}++; $gtypes{$end5}++; }
							if ($up3 != -$down3 && isOverlap($cbeg, $cend, $gend-$_down3, $gend+$up3))    { $types{$end3}++; $gtypes{$end3}++; }
						}
					}
					else
					{
						if    ($cend <=  $gbeg-$up3)                                                     { $types{$out3}++; $gtypes{$out3}++; }
						elsif ($cbeg >=  $gend+$up5)                                                     { $types{$out5}++; $gtypes{$out5}++; }
						else
						{
							if (isOverlap($cbeg, $cend, $gbeg+$_down3, $gend-$_down5))                    { $types{$body}++; $gtypes{$body}++; }
							if ($up5 != -$down5 && isOverlap($cbeg, $cend, $gend-$_down5, $gend+$up5))    { $types{$end5}++; $gtypes{$end5}++; }
							if ($up3 != -$down3 && isOverlap($cbeg, $cend, $gbeg-$up3,    $gbeg+$_down3)) { $types{$end3}++; $gtypes{$end3}++; }
						}
					}

					addGene(\@a, $b[$i], \%gtypes) if $geneFile;
					if (exists $gtypes{$body})
					{
						die "$b[$i]\n" if (exists $gtypes{$out5} || exists $gtypes{$out3});
						delete $gtypes{$body};
					}
					$detail{$id} = join(',', sort keys %gtypes);
					if ($extra)
					{
						$detail{$id} .= "," if $detail{$id};
						$detail{$id} .= $extra;
					}
				}
			}

			if (exists $types{$body} || exists $types{$end5} || exists $types{$end3})
			{
				map {delete $detail{$_} if $detail{$_} =~ /out/} keys %detail;
				map {delete $types{$_} if /out/} keys %types;
			}

			foreach my $type (keys %types)
			{
				$class{$type}++;
			}

			$final{classPriority(keys %types)}++; 

			#print $out join("\t", @a[0..$column], "@types");
			print $out join("\t", @a[0..$column], join('|', $b[0], map {"$_($detail{$_})"} sort(keys(%detail))), join(',', sortPriority(keys %types, keys %subfeats)));

			if (exists $types{$body} && keys %subfeats)
			{
				map {$geneclass{$_}++} sort keys %subfeats;
				$genefinal{bodyPriority(keys %subfeats)}++;
			}

			print $out "\n";
		}
		else
		{
			print $out $_;
		}
	}

	close($in) if defined $fileName;

	print $out "#", $field+1, 'all';
	foreach my $type ($out5, $end5, $body, $end3, $out3, $none)
	{
		print $out "\t", $class{$type} || 0;
	}
	print $out "\n";

	print $out "#", $field+1, 'priority';
	foreach my $type ($out5, $end5, $body, $end3, $out3, $none)
	{
		print $out "\t", $final{$type} || 0;
	}
	print $out "\n";

	print $out "#", $field+1, 'geneall';
	foreach my $type ($utr5, $cds, $utr3, $exon, $intron)
	{
		print $out "\t", $geneclass{$type} || 0;
	}
	print $out "\n";

	print $out "#", $field+1, 'genepriority';
	foreach my $type ($utr5, $cds, $utr3, $exon, $intron)
	{
		print $out "\t", $genefinal{$type} || 0;
	}
	print $out "\n";
}

sub isOverlap 
{
	my ($x1, $x2, $y1, $y2, $allow) = @_;
	$allow = 0 if !defined $allow;

	($x1, $x2) = ($x2, $x1) if ($x1 > $x2);
	($y1, $y2) = ($y2, $y1) if ($y1 > $y2);
	return ($x1 <= $y1 && $y1 < $x2-$allow) || ($x1+$allow < $y2 && $y2 <= $x2) || ($y1 < $x2 && $x2 <= $y2) ? 1 : 0;
	# inclusive-exclusive
}

sub classPriority
{
	return (sortPriority(@_))[0];
}

sub bodyPriority
{
	return (sortPriority(@_))[0];
}

sub sortPriority
{
	map { die "Undefined category in the priority: $_\n" if !exists $pr{$_} } @_;
	return sort {$pr{$a} <=> $pr{$b}} @_;
}

sub addGene
{
	my ($rfeat, $hit, $rtypes) = @_;
	my $ft = classPriority(keys %$rtypes);
	return if $ft eq $none;

	my ($gid, $gbeg, $gend, $dir) = $hit =~ /(.+):(\d+)-(\d+)([+\-])/;
	$dir = $Null if !$dir;
	my $key = "$rfeat->[0]:$gid:$gbeg:$gend:$dir";
	$genes{$key}{'info'} = [$rfeat->[0], $gbeg, $gend, $dir, $gid] if !exists $genes{$key};

	if ($ft eq $body && $' =~ /;([^#]+)/)
	{
		my %subfeats;
		map {s/\d+$//; $subfeats{$_}++} split(/,/, $1);
		$ft = bodyPriority(keys %subfeats);
	}

	push(@{$genes{$key}{$ft}}, join(',', join(':', "$rfeat->[1]-$rfeat->[2]", $bedIndex > 5 ? @$rfeat[3..5] : $bedIndex > 4 ? @$rfeat[3,4] : $bedIndex > 3 ? $rfeat->[3]: '')));
}

sub process
{
	my $out = openOutput($geneFile);
	print $out join("\t", '#ref', 'start', 'end', 'strand', 'ID', $out5, $end5, $utr5, $cds, $intron, $utr3, $end3, $out3), "\n";

	foreach my $key (sort {
			$genes{$a}{'info'}[0] cmp $genes{$b}{'info'}[0] ||
			$genes{$a}{'info'}[1] <=> $genes{$b}{'info'}[1] ||
			$genes{$a}{'info'}[2] <=> $genes{$b}{'info'}[2] ||
			$genes{$a}{'info'}[4] cmp $genes{$b}{'info'}[4]
		} keys %genes)
	{
		my @buf;
		print $out join("\t", @{$genes{$key}{'info'}});

		foreach my $t ($out5, $end5, $utr5, $cds, $intron, $utr3, $end3, $out3)
		{
			print $out "\t", join('|', $t ne $cds ? !exists $genes{$key}{$t   } ? $Null : @{$genes{$key}{$t}}
			                                      : !exists $genes{$key}{$exon} ? !exists $genes{$key}{$t} ? $Null                  : @{$genes{$key}{$t}}
			                                                                    : !exists $genes{$key}{$t} ? @{$genes{$key}{$exon}} : @{$genes{$key}{($t,$exon)}});
		}

		print $out "\n";
	}

	close($out) if $geneFile;
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

	if ($helpFlag || !$field || !$column || !$bedIndex || !$priority)
	{
		die("Arguments: -c last_column_in_output -f field -b bed_index -p priority [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\tprority is a combination of U(upstream), 5(5'UTR), E(exon), I(intron), 3(3'UTR), and D(downstream)\n"
		  );
	}

	$bedIndex--; $field--; $column--;
}
