#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Author: Jeong-Hyeon Choi
## Date: Feb 11 2008

## convert a gff format to a wiggle format

## Usage:
## CHECK usage with option  -h

## $Id: perl_template,v 1.2 2007-09-21 22:51:23 jeochoi Exp $


use constant { WIN=>0, CPG=>1, MET=>2, COV=>3, SUM=>4, VAR=>5, EXT=>6 };

use Getopt::Long qw(:config no_ignore_case);

my ($helpFlag, $inFile, $outFile, $verbose, $quiet, $window, $step);
my ($probeFile, $iRef, $sName, @iNames, $iStart, $iEnd, $iStrand, $upstream, $dnstream);
my $bin = 10000;

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"w|window=i"	=> \$window,
	"s|step=i"		=> \$step,
	"probe=s"		=> \$probeFile,
	"ref=s"			=> \$iRef,
	"name=s"			=> \$sName,
	"start=i"		=> \$iStart,
	"end=i"			=> \$iEnd,
	"strand=i"		=> \$iStrand,
	"upstream=i"	=> \$upstream,
	"downstream=i"	=> \$dnstream,
) || die "\n";

checkOptions();

my (%probes, $id, %refs, $out);
load($probeFile) if $probeFile;
process($inFile);

#-------------------------------------------------------------------------------

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	my $no = 0;

	#chr3RHet	51206	51242	1
	while (<$in>)
	{
		next if /^#/ || /^\s*$/ || /^chrom/;
		my @a = split(/[\t\r\n]/);
		my $pname = ++$no;
		map {$pname .= ":$_"} @a[@iNames] if @iNames;
		$pname =~ s/ //g;

		my ($dir) = defined $iStrand && $a[$iStrand] eq '-' ? 0 : 1;
		my ($wn, $st) = $window   ? ($window, $step) :
		                $upstream ? ($upstream, $upstream) :
		                $dnstream ? ($dnstream, $dnstream) : ($a[$iEnd]-$a[$iStart], $a[$iEnd]-$a[$iStart]);
		my ($left, $right) = ($upstream && $dir == 1) ? getPos($a[$iStart], $upstream, $wn, $st) : 
		                     ($dnstream && $dir == 0) ? getPos($a[$iStart], $dnstream, $wn, $st) : 
		                     ($upstream && $dir == 0) ? ($a[$iEnd], $a[$iEnd]+$upstream) : 
									($dnstream && $dir == 1) ? ($a[$iEnd], $a[$iEnd]+$dnstream) : @a[$iStart, $iEnd];

		for (my $i = $left; ($i+$i+$wn)/2 < $right; $i += $st)
		{
			map {push(@{$probes{$a[$iRef]}{$_}}, [$i, $i+$wn-1, $pname])} int($i/$bin)..int(($i+$wn)/$bin);
			#print "$a[$iRef] $i ", $i+$wn-1, " $pname\n";
		}
	}

	close($in);
}

sub getPos
{
	my ($pos, $stream, $_window, $_step) = @_;
	return ($pos - $_window - $_step * int(($stream-$_window/2) / $_step), $pos);
}

sub process
{
	my ($fileName) = @_;

	print "Processing ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	$out = openOutput($outFile);
	print $out "#Ref\tStart\tEnd\tCpGs\tMethylated\tCoverage\tMeanMethyl\tStdevMethyl";
	print $out "\tID" if $probeFile;
	print $out "\n";

	#chr3RHet	51206	51242	1
	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my @a = split(/[\t\r\n]/);

		onRef() if ($id && $id ne $a[0]);

		$id = $a[0];
		my ($s, $e);
		print if $verbose;

		if ($probeFile)
		{
			my %done;
			foreach my $b (int($a[1]/$bin)..int($a[2]/$bin))
			{
				next if !exists $probes{$a[0]}{$b};
				foreach my $f (@{$probes{$a[0]}{$b}})
				{
					if (!exists $done{$f->[2]}{$f->[0]} && $f->[0] <= $a[1] && $a[1] <= $f->[1])
					{
						add($f->[0], $f->[1], $a[3], $a[4], $f->[2]);
						$done{$f->[2]}{$f->[0]}++;
					}
				}
			}
		}
		else
		{
			$s = ($a[1]+$step-$window) / $step; $s = 0 if $s < 0;
			#$e = ($a[2]-1            ) / $step; # looks like a bug
			$e = ($a[2]-2            ) / $step;

			foreach my $p (int($s)..int($e))
			{
				add($p*$step+1, $p*$step+$window, $a[3], $a[4]);
			}
		}
	}

		onRef() if $id;
	close($in) if defined $fileName;
}

sub add
{
#	my ($start, $end, $pct, $cov, $id) = @_;
	my ($start, $end, $cov, $pct, $id) = @_;
	my $key = $probeFile ? $id : $start;
	$refs{$key}[WIN]  = $end;
	$refs{$key}[CPG] ++;
	$refs{$key}[MET] += int($cov*$pct+0.5);
	$refs{$key}[COV] += $cov;
	$refs{$key}[SUM] += $pct;
	$refs{$key}[VAR] += $pct*$pct;
	$refs{$key}[EXT]  = $probeFile ? $start : $id;
	print "@_ : @{$refs{$key}}\n" if $verbose;
}

sub onRef
{
	foreach my $k (sort {$probeFile ? $refs{$a}[EXT] <=> $refs{$b}[EXT] : $a<=>$b} keys %refs)
	{
		print $out join("\t", $id, ($probeFile ? $refs{$k}[EXT] : $k)-1, @{$refs{$k}}[WIN,CPG,MET,COV]);
		my $m = $refs{$k}[SUM]/$refs{$k}[CPG];
		printf $out "\t%.3f\t%.3f", $m, ($refs{$k}[VAR]/$refs{$k}[CPG] - $m*$m <= 0) ? 0 : sqrt($refs{$k}[VAR]/$refs{$k}[CPG] - $m*$m);
		print $out "\t$k" if $probeFile;
		print $out "\n";
	}

	%refs = ();
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
	open($fd, $fileName =~ /.gz$/ ? "| gzip -c > $fileName" : ">$fileName") || die("Open error: $fileName");
	return $fd;
}

sub checkOptions
{
	$inFile  = shift(@ARGV) if !defined $inFile  && @ARGV > 0;
	$outFile = shift(@ARGV) if !defined $outFile && @ARGV > 0;
	die "The output file name is the same as the input file name\n" if defined $inFile && defined $outFile && $inFile eq $outFile;

	if ($helpFlag || !($window ? $step : $probeFile) || ($probeFile && (!$iRef || !$iStart || !$iEnd || ($upstream && $dnstream))))
	{
		die("Arguments: -s step -s window_size [-p probe_file] [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  );
	}

	if ($probeFile)
	{
		$iRef--; $iStart--; $iEnd--;
		$iStrand-- if $iStrand;
		if ($sName)
		{
			@iNames = map {$_-1} split(/,/, $sName);
		}
	}
}
