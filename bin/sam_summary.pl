#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Author: Jeong-Hyeon Choi
## Date: Feb 8 2010

## summarize sam alignments


## Usage:
## CHECK usage with option  -h

## $Id: perl_template,v 1.2 2007-09-21 22:51:23 jeochoi Exp $

my ($iRead, $iFlag, $iRef, $iPos, $iMapQual, $iInsert, $iMis, $iGap, $iSoft, $iHard, $iJunc, $iLen) = 0..11;	# should be modified if the array structure is changed in the line with ### ARRAY
my ($tSingle, $tRead1, $tRead2) = 0..2;


use Getopt::Long qw(:config no_ignore_case);
use File::Temp;

my ($helpFlag, $samFlag, @inFiles, $outFile, $verbose, $quiet, $mapQual, $seqFile, $sumFile, $ram, $sortFlag, $useMd, $useFlagForUnmapped);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \@inFiles,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"total=i"		=> \$givenTotal,
	"fq=s"			=> \$seqFile,
	"summary=s"		=> \$sumFile,
	"s|sam"			=> \$samFlag,		## input file
	"mq=i"			=> \$mapQual,
	"ram=i"			=> \$ram,
	"sort"			=> \$sortFlag,
	"md!"				=> \$useMd,
	"flag!"			=> \$useFlagForUnmapped,
) || die "\n";

checkOptions();

my %rc = ('A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', 'N' => 'N', 'a' => 't', 't' => 'a', 'c' => 'g', 'g' => 'c', 'n' => 'n');
my %cigarIndex = ('I' => 0, 'D' => 1, 'S' => 2, 'H' => 3, 'M' => 4, 'N' => 5);
my $out;
my %refs;
my %pairmap;
my ($total, $duplicate, $failure, $mapped) = (0, 0, 0);
my ($single, $paired, $read1, $read2) = (0, 0, 0, 0);
my ($funique, $runique, $unique, $multi, $lowUni, $lowMul) = (0, 0, 0, 0, 0, 0);
my @inserts = (0, 0, 0, 0); # diff_ref improper proper sum
my (%mismatches, %indels, %errors, %softs, %hards);
my (%locs, %lens, %ulens);
my ($nmFlag, $mdFlag, $nmPresent, $mdPresent);
my $tmp = File::Temp->new(OPEN => 0);
my $tmpFile = $tmp->filename;

if (!$givenTotal)
{
	if    ($seqFile) { $givenTotal = loadSeq($seqFile); }
	elsif ($sumFile) { $givenTotal = loadSum($sumFile); }
}

if (@inFiles) { map {load($_)} @inFiles; }
else { load(undef); }
output();

warn "NM tags are absent\n" if !$nmPresent;
warn "MD tags are absent\n" if !$mdPresent;

#-------------------------------------------------------------------------------

sub onRead
{
	my ($rmatches) = @_;
	my $qid = $rmatches->[0][0];
	my %ids;

	foreach my $ra (@{$rmatches})
	{
		$ra->[$iFlag] = int($ra->[$iFlag]);
		my $stype = type($ra->[$iFlag]);
		die "The read $qid is a paired read, but not set either read 1 or 2\n" .  join("\n", @{$rmatches}) . "\n" if $stype == 0 && paired($ra->[$iFlag]);
		push(@{$ids{$stype}}, \@$ra);
	}

	if (exists $ids{$tSingle})
	{
		die "The read $qid is both single and paired\n" . join("\n", @{$rmatches}) . "\n" if exists $ids{$tRead1} || exists $ids{$tRead2};
		$single++;
		my $n = mapStat($ids{$tSingle});
		refStat($n, $ids{$tSingle});
	}
	else
	{
		$paired++;
		$read1++ if exists $ids{$tRead1};
		$read2++ if exists $ids{$tRead2};
		my $n1 = mapStat($ids{$tRead1});
		my $n2 = mapStat($ids{$tRead2});
		refStat($n1, $ids{$tRead1});
		refStat($n2, $ids{$tRead2});
		$pairmap{"$n1=$n2"}++;
		if ($n1 == 1 && $n2 == 1)
		{
			if    ($ids{$tRead1}[0][$iRef] ne $ids{$tRead2}[0][$iRef]) { $inserts[0]++; }
			elsif ($ids{$tRead1}[0][$iInsert] == 0)                    { $inserts[1]++; }
			else                                                       { $inserts[2]++; $inserts[3] += $ids{$tRead1}[0][$iInsert]; }
		}
	}
}

sub mapped
{
#	print " $_[0] & 0x4 ? 0 : 1 = ", $_[0] & 0x4 ? 0 : 1, "\n" if $_[0] == 4;
	return $_[0] & 0x4 ? 0 : 1;
}

sub paired
{
	return $_[0] & 0x1 ? 1 : 0;
}

sub type
{
	return $_[0] & 0x40 ? $tRead1 : $_[0] & 0x80 ? $tRead2 : $tSingle;
}

sub strand
{
	return $_[0] & 0x10 ? 1 : 0;
}

sub mapStat
{
	my ($ra) = @_;

	if (!defined $ra)
	{
		$code = -1;
	}
	elsif (@$ra > 1)
	{
		$code = 3;
		$multi++;
	}
	elsif (!mapped($ra->[0][$iFlag]))
	{
		$code = 0;
	}
	elsif (@$ra > 1)
	{
		$code = 2;
		$multi++;
		$lowMul++ if ($ra->[0][$iMapQual] <= $mapQual);
	}
	else
	{
		$code = 1;
		$lowUni++ if ($ra->[0][$iMapQual] <= $mapQual);
		strand($ra->[0][$iFlag]) ? $runique++ : $funique++;
		$mismatches{$ra->[0][$iMis]}++;
		$indels{$ra->[0][$iGap]}++;
		$errors{$ra->[0][$iMis]+$ra->[0][$iGap]}++;
		$softs{$ra->[0][$iSoft]}++;
		$hards{$ra->[0][$iHard]}++;
		$ulens{$ra->[0][$iLen]}++;
	}

	$locs{@$ra}++            if ($code > 0);
	$lens{$ra->[0][$iLen]}++ if ($code > 0);
	$total++                 if ($code >= 0);

	return $code;
}

sub refStat
{
	my ($code, $ra) = @_;
	return if $code == 0;

	foreach my $rb (@$ra)
	{
		$refs{$rb->[$iRef]}[($code == 1 ? 0 : 1)*2+strand($rb->[$iFlag])]++;
	}
}

sub output
{
	$out = openOutput($outFile);

	print $out "#Read information\n";
	print $out "TotalGiven\t$givenTotal\n" if $givenTotal;
	print $out "TotalInFile\t$total\n";
	print $out "Single\t$single\n" if $single;
	if ($paired)
	{
		print $out "Paired\t$paired\n";
		print $out "Read1\t$read1\n";
		print $out "Read2\t$read2\n";
	}

	print $out "\n#Mapping information\n";
	print $out "#Type\tSum\t%Total\t%Mapped\tForward\tReverse\n";
	$unique = $funique + $runique;
	$mapped = $unique + $multi;
	print $out "TotalMapped\t$mapped\t", percent($mapped), "\n";
	print $out join("\t", "UniquelyMapped", $funique+$runique, percent($funique+$runique, $mapped), $funique, $runique), "\n";
	print $out join("\t", "U   MapQual<=$mapQual", $lowUni, percent($lowUni, $mapped)), "\n";
	print $out join("\t", "MultiplyMapped", $multi, percent($multi, $mapped)), "\n";
	print $out join("\t", "M   MapQual<=$mapQual", $lowMul, percent($lowMul, $mapped)), "\n";

	outputPair() if $paired;

	outputHash('Location by all mapped reads', \%locs, $mapped);
	outputHash('Length of all mapped reads', \%lens, $mapped);
	outputHash('Length of uniquely mapped reads', \%ulens, $mapped);
	outputHash('Mismatch by uniquely mapped reads', \%mismatches, $unique) if $unique;
	outputHash('Indel by uniquely mapped reads', \%indels, $unique) if $unique;
	outputHash('Error by uniquely mapped reads', \%errors, $unique) if $unique;
	outputHash('Soft clip by uniquely mapped reads', \%softs, $unique) if $unique && keys %softs;
	outputHash('Hard clip by uniquely mapped reads', \%hards, $unique) if $unique && keys %hards;

	outputRef('uniquely', 0, 1) if $unique;
	outputRef('multiply', 2, 3) if $multi;

	print $out "\n";
}

sub outputPair
{
	print $out "\n#Pair information\n";
	print $out "#number\tSum\t%Total\t%Mapped\n";
	my $pairMapped = 0;

	foreach my $key (sort keys %pairmap)
	{
		$pairMapped += $pairmap{$key} if $key =~ /[12]/;
	}

	foreach my $key (sort keys %pairmap)
	{
		my $v = $pairmap{$key} || 0;
		$key =~ s/-1/U/g;
		$key =~ s/0/N/g;
		$key =~ s/2/m/g;
		$key =~ s/3/M/g;
		print $out join("\t", $key, $v, pairPercent($v, $pairMapped)), "\n";
	}
}

sub outputHash
{
	my ($title, $hash, $basis) = @_;

	return if keys %$hash == 0;
	print $out "\n#$title\n";
	print $out "#number\tSum\t%Total\t%Mapped\n";

	my @ns = sort {$a <=> $b} keys %$hash;
	foreach my $n (0..$ns[-1])
	{
		next if $n == 0 && $title =~ /location/i;
		my $v = $hash->{$n} || 0;
		print $out join("\t", $n, $v, percent($v, $basis)), "\n"; 
	}
}

sub outputRef
{
	my ($title, @indices) = @_;

	print $out "\n#Reference information by $title mapped reads\n";
	print $out "#Ref\tSum\t%Total\t%Mapped\tForward\tReverse\n";
	foreach my $rid (sort {
			my $aa = $a =~ /(\d+)/ ? $1 : 100000;
			my $bb = $b =~ /(\d+)/ ? $1 : 100000;
			$aa <=> $bb || $a cmp $b
		} keys %refs)
	{
		my @n = map {$refs{$rid}[$_]||0} @indices;
		my $s = 0; map {$s += $_} @n;
		print $out join("\t", $rid, $s, percent($s, $indices[0]==0 ? $unique : $multi), @n), "\n";
	}
}

sub percent
{
	my $s = sprintf("%.2f", 100*$_[0]/($givenTotal||$total));
	$s .= "\t" . sprintf("%.2f", 100*$_[0]/$_[1]) if $_[1];
	return $s
}

sub pairPercent
{
	my $s = sprintf("%.2f", 100*$_[0]/$paired);
	$s .= "\t" . sprintf("%.2f", 100*$_[0]/$_[1]) if $_[1];
	return $s
}

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my ($in, $cmd, @matches);

	if ($sortFlag)
	{
		$cmd = "samtools view -bhS " . ($fileName || '-') . ' | ' if $samFlag;
		$cmd .= "samtools sort -n -m $ram " . (!$samFlag && $fileName ? $fileName : '-') . " $tmpFile";
		run($cmd);

		$cmd = "samtools view -h $tmpFile.bam |";
		$in = openInput($cmd);
	}
	else
	{
		$in = openInput($samFlag ? $fileName : "samtools view -h $fileName |");
	}

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my ($qid, $flag, $rid, $pos, $mapq, $cigar, $mid, $mpos, $isize) = split(/[\t\r\n]/);

		if (/^\@\s*(\S+)/)
		{
			if ($1 eq 'SQ')
			{
				if (/SN:(\S+)/ && !exists $refs{$1}) { $refs{$1} = [0, 0, 0, 0]; }
			}
		}
		else
		{
			if (@matches && $matches[0][0] ne $qid)
			{
				onRead(\@matches);
				@matches = ();
			}

			$flag = int($flag);
			if ($flag & 0x400) # d=0x400 (duplicate)
			{
				$duplicate++;
			}
			elsif ($flag & 0x200) # f=0x200 (failure)
			{
				$failure++;
			}
			else
			{
				#push(@{$reads{$qid}}, [$flag, $rid, $pos, $mapq, $mpos, $isize]);
				die "Wrong condition for mapping: $_\n" if !$useFlagForUnmapped && mapped($flag) != ($cigar ne '*' ? 1 : 0);
				my @cigars = parseCigar($cigar);
				my ($nm, $mis, $del);

				if (/\bNM:i:(\d+)/i) { $nm = $1; $nmFlag = 1; $nmPresent = 1; }
				else                 { $nm =  0; $nmFlag = 0; }

				if (/\bMD:Z:(\S+)/i)
				{
					($mis, $del) = misdel2($1); $mdFlag = 1;
					if ($del != $cigars[1])
					{
						warn "Inconsistent number of deletions: $del == $cigars[1]\n$_";
					}
					if (!$useMd && $nmFlag && $nm != $mis+$del+$cigars[0])
					{
						warn "Inconsistent number of errors: $nm == $mis+$del+$cigars[0]\n$_";
					}
					$mdPresent = 1;
				}
				else
				{
					($mis, $del) = (0, 0); $mdFlag = 0;
				}

				push(@matches, [$qid, $flag, $rid, $pos, $mapq, $isize, $mis, $cigars[0]+$cigars[1], $cigars[2], $cigars[3], $cigars[4], $cigars[0]+$cigars[4]]);	### ARRAY
			}
		}
	}

	onRead(\@matches) if @matches;
	die "Inconsistent 1:1 paired reads\n" if exists $pairmap{'1:1'} && $pairmap{'1:1'} != $inserts[0]+$inserts[1]+$inserts[2];
#	close($in) if $in != STDIN
}

# 8M1D46M
sub parseCigar
{
	my ($cigar) = @_;
	my @cigars = (0, 0, 0, 0, 0);

	while ($cigar =~ /(\d+)([IDSHMN])/ig) { $cigars[$cigarIndex{$2}] += $1; }

	return @cigars;
}

# [0-9]+(([ACGTN]|\^[ACGTN]+)[0-9]+)
sub misdel
{
	my ($md) = @_;
	my @a = split(/[\^A-Za-z]/, $md);
	my ($mis, $del) = (0, 0);

	for (my $i = 0; $i < $#a; $i++)
	{
		if    ($a[$i] eq '') { $del++; $mis-- if $a[$i-1] ne ''; } # deletion
		else                 { $mis++; }
	}

	#print "$md (@a) $mis $del\n";
	return ($mis, $del);
}

sub misdel2
{
	my ($md) = @_;
	my ($mis, $del) = (0, 0);

	while($md =~ /(\^[A-Za-z]+|[A-Za-z]+)/g){
		my $bases = $1;
		if($bases =~ /^\^/) {
			$del += length($');
		}
		else{
			$mis += length($bases);
		}
	}

	return ($mis, $del);
} 

sub loadSeq
{
	my ($path) = @_;
	my $no = 0;

	foreach my $fileName (split(/,/, $path))
	{
		print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
		my $in = openInput($fileName);

		while (<$in>) { $no++; }

		close($in) if defined $fileName;
	}

	return $no/4;
}

sub loadSum
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	my $no;

	while (<$in>)
	{
		if (/^TotalGiven\s+(\d+)/) { $no = $1; }
		elsif (/^TotalInFile\s+(\d+)/) { $no = $1 if !$no; last; }
	}

	close($in) if defined $fileName;
	return $no;
}

sub run
{
	my ($cmd) = @_;
	#print "$cmd\n";
	!system($cmd) || die "Error in $cmd\n";
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
	push(@inFiles, @ARGV) if !@inFiles  && @ARGV > 0;

	if ($helpFlag)
	{
		die("Arguments: [-t no_reads] [-sam] [-mq min_maq_qual] [-s] in_file in_file2 .. [-o out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}

	$mapQual = 0 if !$mapQual;
	$ram = 1000000000 if !$ram;
}
