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
my ($up5, $down5, $up3, $down3) = (0, 0, 0, 0);
my ($out5, $end5, $body, $end3, $out3, $outside, $none) = ("farup", "upstream", "body", "downstream", "fardown", "intergenic", "none");
my %indexes = ($out5 => 3, $end5 => 4, $body => 5, $end3 => 6, $out3 => 7, $none => 8);
my %rindexes; map {$rindexes{$indexes{$_}} = $_} keys %indexes;

use Getopt::Long qw(:config no_ignore_case);

my ($helpFlag, $inFile, $outFile, $sumFile, $verbose, $quiet, $field, $classFlag, $scoreFlag, $binSize);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"field=i"		=> \$field,
	"quiet"			=> \$quiet,
	"up5=i"			=> \$up5,
	"down5=i"		=> \$down5,
	"up3=i"			=> \$up3,
	"down3=i"		=> \$down3,
	"class!"			=> \$classFlag,
	"score!"			=> \$scoreFlag,
	"bin=f"			=> \$binSize,
	"category"		=> \$category,
	"summary=s"		=> \$sumFile,
) || die "\n";

checkOptions();

my (%class, %final, %fina2, %infos, %score, %infoScore, %classScore, %feats, %feat2, %comb);

load($inFile);
summary() if $sumFile;


#-------------------------------------------------------------------------------

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	my $out = openOutput($outFile);

	while (<$in>)
	{
		next if /^\s*$/;
		my @a = split(/[\t\r\n]/);

		if (/^#/)
		{
			print $out $_;
		}
		else
		{
			my (@dats, %ihash, %phash, $flag, $dist);
			my $sc = $scoreFlag ? int($a[4]/$binSize) : 0;

			if ($a[$field] eq $Null)
			{
				push(@dats, [$Null, $none]); # ID and feature
				$dist = $Null;
			}
			else
			{
				my ($cbeg, $cend) = @a[1,2];
				my @b = split(/\|/, $a[$field]);
				$dist = $b[0];

				for (my $i = 1; $i < @b; $i++)
				{
					my ($id, $fbeg, $fend, $dir) = $b[$i] =~ /(.+):(\d+)-(\d+)([+\-])*/;
					die "Unknown id: $a[$field]\n" if !$id;
#					$a[$field] = "-$a[$field]" if $a[2] < $fbeg;
					$id =~ s/\_\d+$// if $category;
					my $part = classify($dir, $cbeg, $cend, $fbeg, $fend);
					push(@dats, [$id, $part, $fbeg, $fend, $id, $dir||$Null]);
					$flag = 1 if $part eq $end5 || $part eq $body || $part eq $end3;
				}

				for (my $i = 0; $i < @dats; $i++)
				{
					next if ($flag && ($dats[$i][1] eq $out5 || $dats[$i][1] eq $out3));

					if ($b[$i+1] =~ /#(\S+)/)
					{
						$infos{$1}{$dats[$i][1]}++;
						$infoScore{$sc}{$1}{$dats[$i][1]}++ if $scoreFlag;
					}

					$ihash{$dats[$i][0]} = $i;
					$phash{$dats[$i][1]}++;
				}

				if ($flag)
				{
					my @temp;
					map { my $p = $dats[$ihash{$_}][1]; if ($p eq $end5 || $p eq $body || $p eq $end3) { push(@temp, $dats[$ihash{$_}]); } } keys %ihash;
					@dats = @temp;
				}

				if ($sumFile)
				{
					foreach my $d (@dats)
					{
						push(@refs, $a[0]) if !$feats{$a[0]};
						@{$feats{$a[0]}{$d->[2]}} = @$d[3..5] if !exists $feats{$a[0]}{$d->[2]};
						$feats{$a[0]}{$d->[2]}[$indexes{$d->[1]}] += 1/keys(%ihash);
						$feat2{$a[0]}{$d->[2]}[$indexes{$d->[1]}] ++;
					}
				}
			}

#			print $out join("\t", @a[0..$#a-1],
#				join(',', map { $_->[0] } @dats),
#				join(',', map { $_->[1] } @dats),
#				$dist), "\n";
			print $out join("\t", @a[0..$#a-1], join(',', map { $_->[1] ne $none ? "$_->[0]($_->[1])" : $_->[0] } @dats), keys %phash, $dist), "\n";

			if ($a[$field] eq $Null)
			{
				$final{$none}++;
				$fina2{$none}++;
				$comb{$none}++;
			}
			else
			{
				for my $d (@dats)
				{
					if ($classFlag)
					{
						$class{$d->[0]}{$d->[1]}++;
						$classScore{$sc}{$d->[0]}{$d->[1]}++ if $scoreFlag;
					}
					$score{$sc}{$d->[1]}++ if $scoreFlag;
					$final{$d->[1]} += 1/keys(%ihash);
				}
				map { $fina2{$_}++ } keys(%phash);
				$comb{join(',', sort {$indexes{$a} <=> $indexes{$b}} keys %phash)}++;
			}
		}
	}

	close($in) if defined $fileName;

	my $log = openOutput($outFile ? "$outFile.stat" : $outFile);
	print $log "## INFO\n" if %infos;
	foreach my $inf (sort {$a =~ /^\d+/ && $b =~ /^\d+/ ? $a <=> $b : $a cmp $b} keys %infos)
	{
		my $sum = 0;
		print $log "$inf";

		foreach my $part ($out5, $end5, $body, $end3, $out3)
		{
			print $log "\t", $infos{$inf}{$part} || 0;
			$sum += $infos{$inf}{$part} || 0;
		}

		print $log "\t$sum\n";
	}

	if ($scoreFlag)
	{
		print $log "## SCORES\n";
		print $log join ("\t", '#', $out5, $end5, $body, $end3, $out3), "\n";
		foreach my $s (sort {$a <=> $b} keys %score)
		{
			my $sum = 0;
			print $log $s*$binSize;

			foreach my $part ($out5, $end5, $body, $end3, $out3)
			{
				print $log "\t", $score{$s}{$part} || 0;
				$sum += $score{$s}{$part} || 0;
			}

			print $log "\t$sum\n";
		}

		print $log "## SCOREINFO\n" if %infoScore;
		print $log join ("\t", '#', $out5, $end5, $body, $end3, $out3), "\n";
		foreach my $inf (sort {$a =~ /^\d+/ && $b =~ /^\d+/ ? $a <=> $b : $a cmp $b} keys %infos)
		{
			foreach my $s (sort {$a <=> $b} keys %score)
			{
				my $sum = 0;
				print $log "$inf\t", $s*$binSize;

				foreach my $part ($out5, $end5, $body, $end3, $out3)
				{
					print $log "\t", $infoScore{$s}{$inf}{$part} || 0;
					$sum += $infoScore{$s}{$inf}{$part} || 0;
				}

				print $log "\t$sum\n";
			}
		}

		if (%classScore)
		{
			print $log "## SCOREFEATURES\n";
			print $log join ("\t", '#', $out5, $end5, $body, $end3, $out3), "\n";
			foreach my $type (sort keys %class)
			{
				foreach my $s (sort {$a <=> $b} keys %score)
				{
					my $sum = 0;
					print $log "$type\t", $s*$binSize;

					foreach my $part ($out5, $end5, $body, $end3, $out3)
					{
						print $log "\t", $classScore{$s}{$type}{$part} || 0;
						$sum += $classScore{$s}{$type}{$part} || 0;
					}

					print $log "\t$sum\n";
				}
			}
		}
	}

	if (%class)
	{
		print $log "## FEATURES\n";
		print $log join ("\t", '#', $out5, $end5, $body, $end3, $out3), "\tsum\n";
		foreach my $type (sort keys %class)
		{
			my $sum = 0;
			print $log "$type";

			foreach my $part ($out5, $end5, $body, $end3, $out3)
			{
				print $log "\t", $class{$type}{$part} || 0;
				$sum += $class{$type}{$part} || 0;
			}

			print $log "\t$sum\n";
		}

		print $log join ("\t", '#', $end5, $body, $end3, $outside, $none, 'sum'), "\n";

		my $sum = 0;
		print $log "total";
		foreach my $part ($end5, $body, $end3, $out5, $none)
		{
			my $c = $final{$part} || 0;
			$c += $final{$out3} || 0 if $part eq $out5;
			$sum += $c;
			print $log "\t$c";
		}
		print $log "\t$sum\n";

		$sum = 0;
		print $log "total2";
		foreach my $part ($end5, $body, $end3, $out5, $none)
		{
			my $c = $fina2{$part} || 0;
			$c += $fina2{$out3} || 0 if $part eq $out5;
			$sum += $c;
			print $log "\t$c";
		}
		print $log "\t$sum\n";
	}

	print $log "#combination\n";
	map {print $log "$_\t$comb{$_}\n"} sort keys %comb;
}

sub classify
{
	my ($dir, $cbeg, $cend, $gbeg, $gend) = @_;
	my $type;

	if (!$dir || $dir eq '+')
	{
		if    ($cbeg >= $gend+$up3                           ) { $type = $out3; }
		elsif ($cend <= $gbeg-$up5                           ) { $type = $out5; }
		elsif ($cbeg <  $gend-$down3 && $cend >  $gbeg+$down5) { $type = $body; }
		elsif ($cbeg <= $gbeg+$down5 && $cend >  $gbeg-$up5  ) { $type = $end5; }
		else                                                   { $type = $end3; }
	}
	else
	{
		if    ($cend <= $gbeg-$up3                           ) { $type = $out3; }
		elsif ($cbeg >= $gend+$up5                           ) { $type = $out5; }
		elsif ($cbeg <  $gend-$down5 && $cend >  $gbeg+$down3) { $type = $body; }
		elsif ($cbeg <  $gend+$up5   && $cend >= $gend-$down5) { $type = $end5; }
		else                                                   { $type = $end3; }
	}

	return $type;
}

sub summary
{
	my $out = openOutput($sumFile);
	print $out join("\t", '#Ref', 'Start', 'End', 'Strand', (sort {$indexes{$a} <=> $indexes{$b}} keys %indexes), 'Sum'), "\n";
	my (%count, @total, @tota2);

	foreach my $rid (@refs)
	{
		foreach my $p (sort {$a<=>$b} keys %{$feats{$rid}})
		{
			print $out join("\t", $rid, $p, @{$feats{$rid}{$p}}[0..2]);
			my ($sum, $su2) = (0, 0);
			my @ts;

			foreach my $i (3..8)
			{
				my $n = $feats{$rid}{$p}[$i] || 0; # 1
				my $m = $feat2{$rid}{$p}[$i] || 0; # 1/n
				print $out "\t$m";
				push(@ts, $rindexes{$i}) if $m;
				$sum += $m;
				$su2 += $n;
				$total[$i] += $n;
				$tota2[$i] += $m;
			}

			print $out "\t$sum\n";
			$total[9] += $su2;
			$tota2[9] += $sum;
			$count{join(",", sort {$indexes{$a} <=> $indexes{$b}} @ts)}++;
		}
	}

	print $out join ("\t", '#', $out5, $end5, $body, $end3, $out3, 'none', 'sum'), "\n";
	print $out join("\t", '#total' , map {$_||0} @total[3..9]), "\n";
	print $out join("\t", '#total2', map {$_||0} @tota2[3..9]), "\n";
	map {print $out "#count\t$_\t$count{$_}\n"} sort keys %count;
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

	if ($helpFlag || !$field)
	{
		die("Arguments: -f field [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}

	$field--;
	$classFlag = 1 if !defined $classFlag;
	$binSize = 1 if !defined $binSize;
}
