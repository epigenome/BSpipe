#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Author: Jeong-Hyeon Choi
## Date: May 2 2008

## Given a locus, add the nearest feature in each strand

## Usage:
## CHECK usage with option  -h

## $Id: perl_template,v 1.2 2007-09-21 22:51:23 jeochoi Exp $


use Getopt::Long qw(:config no_ignore_case);

my $Null = '.';
my ($exon, $intron, $utr5, $cds, $utr3) = ('exon', 'intron', "5UTR", 'CDS', "3UTR");
my ($helpFlag, $inFile, $outFile, $verbose, $quiet);
my ($featureFile, $name, $iRef, $idStr, $iStart, $iEnd, $iDir, $unique, $iExonStart, $iExonEnd, $iCdsStart, $iCdsEnd, $iInfo, $allow, $headerFlag);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"i|input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"feature=s"		=> \$featureFile,
	"n|name=s"		=> \$name,
	"ref=i"			=> \$iRef,
	"id=s"			=> \$idStr,
	"start=i"		=> \$iStart,
	"end=i"			=> \$iEnd,
	"dir=i"			=> \$iDir,
	"unique!"		=> \$unique,
	"es=i"			=> \$iExonStart,
	"ee=i"			=> \$iExonEnd,
	"cs=i"			=> \$iCdsStart,
	"ce=i"			=> \$iCdsEnd,
	"info=i"			=> \$iInfo,
	"allow=i"		=> \$allow,
	"header"			=> \$headerFlag,
) || die "\n";

checkOptions();

my %features;
my %ids;
my $break = 10000;

my ($tId, $tStart, $tEnd) = 0..2;
my $tDir = defined $iDir ? $tEnd+1 : $tEnd;
my $tExonStart = defined $iExonStart ? $tDir+1 : $tDir;
my $tExonEnd = defined $iExonEnd ? $tExonStart+1 : $tExonStart;
my $tCdsStart = defined $iCdsStart ? $tExonEnd+1 : $tExonEnd;
my $tCdsEnd = defined $iCdsEnd ? $tCdsStart+1 : $tCdsStart;
my $tInfo = defined $iInfo ? $tCdsEnd+1 : $tCdsEnd;
my @iId = map {--$_} split(/,/, $idStr);

load($featureFile);
process($inFile);

#-------------------------------------------------------------------------------

sub getId
{
	my ($s) = @_;

	$s =~ s/ //g;
	return $unique ? $s : "${s}_".++$ids{$s};
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
		my $id = getId(defined $idStr ? join(":", @a[@iId]) : $name);

		my @v;
		push(@v, $id, @a[$iStart, $iEnd]);
		push(@v, $a[$iDir]) if defined $iDir;
		push(@v, $a[$iExonStart]) if defined $iExonStart;
		push(@v, $a[$iExonEnd]) if defined $iExonEnd;
		push(@v, @a[$iCdsStart, $iCdsEnd]) if (defined $iCdsStart && defined $iCdsEnd && $a[$iCdsStart] != $a[$iCdsEnd]);
		push(@v, $a[$iInfo]) if defined $iInfo;

		for (my $i = int(($a[$iStart]-$allow)/$break); $i <= int(($a[$iEnd]+$allow)/$break); $i++)
		{
			push(@{$features{$a[$iRef]}{$a[$iDir]}{$i}}, \@v);
			#print "$a[$iRef] $i\n";
		}

		$features{$a[$iRef]}{$a[$iDir]}{S} = int(($a[$iStart]-$allow)/$break) if !exists $features{$a[$iRef]}{$a[$iDir]}{S} || $features{$a[$iRef]}{$a[$iDir]}{S} > int(($a[$iStart]-$allow)/$break);
		$features{$a[$iRef]}{$a[$iDir]}{E} = int(($a[$iEnd  ]+$allow)/$break) if !exists $features{$a[$iRef]}{$a[$iDir]}{E} || $features{$a[$iRef]}{$a[$iDir]}{E} < int(($a[$iEnd  ]+$allow)/$break);
		}

	close($in) if defined $fileName;
}

sub process
{
	my ($fileName) = @_;

	print "Processing ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	my $out = openOutput($outFile);

	if ($headerFlag)
	{
		my $header = <$in>; $header =~ s/$/\t$name(Top)\t$name(Bot)/;
		print $out $header;
	}

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		chomp;
		my @a = split(/\t/, $_, -1);
		my @fs = search(@a[0,1,2], '+');
		my @ts = search(@a[0,1,2], '-');
		print $out join("\t", @a, join('|', @fs), join('|', @ts)), "\n";
	}

	close($in) if defined $fileName;
}

sub search
{
	my ($id, $start, $end, $strand) = @_;
	print "$id $start $end $strand\n" if $verbose;
	return ($Null) if !exists $features{$id}{$strand}{S};

	my (%done);
	my @ret1 = search2($features{$id}{$strand}, int($start/$break)  , $features{$id}{$strand}{E},  1, $start, $end, \%done);
	my @ret2 = search2($features{$id}{$strand}, int($start/$break)-1, $features{$id}{$strand}{S}, -1, $start, $end, \%done);
	my @ret = !@ret1 ? @ret2 : !@ret2 ? @ret1 : $ret1[0] < $ret2[0] ? @ret1 : $ret1[0] > $ret2[0] ? @ret2 : merge(\@ret1, \@ret2);

	print "result: @ret\n" if @ret && $verbose;
	return @ret ? @ret : ($Null);
}

sub merge
{
	my ($rb, $rf) = @_;
	my %hash;
	foreach my $f (@$rb[1..$#$rb], @$rf[1..$#$rf]) { $hash{$f}++; }
	return ($rb->[0], keys %hash);
}

sub search2
{
	my ($rfeat, $ij, $ik, $del, $start, $end, $rdone) = @_;
	my (@ret, $dist);

	for (my $i = $ij; $del > 0 ? $i <= $ik : $i >= $ik; $i += $del)
	{
		next if !exists $rfeat->{$i};
		print "    bin: $i\n" if $verbose;
		searchBin($rfeat->{$i}, $start, $end, $rdone, \@ret, \$dist);
		last if defined $dist && ($del > 0 ? $dist < 0 : $dist > 0); # this is last bin
	}

	print "    @ret\n" if $verbose;
	return @ret;
}

sub searchBin
{
	my ($rbins, $start, $end, $rdone, $rret, $rdist) = @_;

	foreach my $f (@$rbins)
	{
		next if exists $rdone->{$f};
		print "\t$f @$f\n" if $verbose > 1;
		$rdone->{$f}++;
		my ($d, @body);

		if (isOverlap($f->[$tStart], $f->[$tEnd], $start, $end))
		{
			$d = 0;
			searchExon($f, $start, $end, \@body) if defined $iExonStart && $iExonEnd;
		}
		elsif ($f->[$tEnd] <= $start)
		{
			$d = $start - $f->[$tEnd] + 1; # inclusive-exclusive
		}
		else # ($end <= $f->[$tStart])
		{
			$d = $end - $f->[$tStart] - 1; # inclusive-exclusive
		}

		next if (!defined $d);
		if (!@$rret || (abs($d) < $rret->[0] && abs($d) >= $allow))
		{
			@$rret = ();
			push(@$rret, abs($d), featureStr($f, \@body));
		}
		elsif (abs($d) == $rret->[0] || abs($d) <= $allow)
		{
			$rret->[0] = abs($d) if $rret->[0] > abs($d);
			push(@$rret, featureStr($f, \@body));
		}

		print "\tR: @$rret\n" if $verbose > 1;
		$$rdist = $d if !defined $$rdist || abs($$rdist) < abs($d);
	}
}

sub searchExon
{
	my ($f, $start, $end, $rbody) = @_;
	my @ss = split(/,/, $f->[$tExonStart]);
	my @es = split(/,/, $f->[$tExonEnd]);
	my ($lutr, $rutr, $cds) = $f->[$tDir] eq '-' ? ($utr3, $utr5, $cds) : ($utr5, $utr3, $cds);
	my ($cno, @cdss) = (0);
	print "CDS: $f->[$tCdsStart] $f->[$tCdsEnd]\n" if ($f->[$tCdsStart] && $f->[$tCdsEnd]) && $verbose;

	for (my $i = 0; $i < @ss; $i++)
	{
		my ($eno, $ino) = $f->[$tDir] eq '-' ? (@ss-$i, @ss-$i) : ($i+1, $i);
		my ($lno, $rno) = ($i+1, @ss-$i);
		print "if (isOverlap($start, $end, $ss[$i], $es[$i]))" if $verbose > 2;

		if ($f->[$tCdsStart] && $f->[$tCdsEnd] && isOverlap($f->[$tCdsStart], $f->[$tCdsEnd], $ss[$i], $es[$i]))
		{
			$cno = $cno ? $cno+1 : 1;
		}

		if (isOverlap($start, $end, $ss[$i], $es[$i]))
		{
			my $cflag;
			if (!$f->[$tCdsStart] || !$f->[$tCdsEnd])
			{
			}
			elsif ($es[$i] <= $f->[$tCdsStart])
			{
				push(@$rbody, "$lutr$lno");
			}
			elsif ($ss[$i] >= $f->[$tCdsEnd])
			{
				push(@$rbody, "$rutr$rno");
			}
			elsif ( ($f->[$tCdsStart] <= $ss[$i] && $es[$i] <= $f->[$tCdsEnd])
				  || ($f->[$tCdsStart] <= $start  && $end    <= $f->[$tCdsEnd]) )
			{
				$cflag = 1;
			}
			else
			{
				if ($ss[$i] < $f->[$tCdsStart])
				{
					if ($start < $f->[$tCdsStart]                           ) { push(@$rbody, "$lutr$lno"); }
					if ($end   > $f->[$tCdsStart] && $start < $f->[$tCdsEnd]) { $cflag = 1; }
				}
				if ($es[$i] > $f->[$tCdsEnd])
				{
					if ($end   > $f->[$tCdsEnd]                           ) { push(@$rbody, "$rutr$rno"); }
					if ($start < $f->[$tCdsEnd] && $end > $f->[$tCdsStart]) { $cflag = 1; }
				}
			}

			push(@cdss, $cno) if ($cflag);
			push(@$rbody, "exon$eno");
			print " => @$rbody" if $verbose > 2;
		}

		push(@$rbody, "intron$ino") if $i && isOverlap($start, $end, $es[$i-1]+1, $ss[$i]);
		print "\n" if $verbose > 2;
	}

	map {push(@$rbody, $cds . ($f->[$tDir] eq '-' ? $cno-$_+1 : $_))} @cdss;
}

sub featureStr
{
	my ($f, $g) = @_;

	map {print "\t\tG: $_\n"} sort @$g if @$g && $verbose;
	my $ret = "$f->[$tId]:$f->[$tStart]-$f->[$tEnd]";
	$ret .= $f->[$tDir] if defined $iDir;
	$ret .= ';' . join(',', sort @$g) if @$g;
	$ret .= '#' . $f->[$tInfo] if defined $iInfo && $f->[$tInfo];
	print "\tF: $ret\n" if $verbose;
	return $ret;
}

sub isOverlap 
{
	my ($l1, $r1, $l2, $r2) = @_;
	return ($l1 <= $l2 && $l2 < $r1) || ($l1 < $r2 && $r2 <= $r1) || ($l2 <= $l1 && $r1 <= $r2) ? 1 : 0;
	# inclusive-exclusive
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

	if ($helpFlag || (!$featureFile && !$inFile) || !$name || !$iRef || !$iStart || !$iEnd || !$idStr || !$iDir)
	{
		die("Arguments: -name feature_name -ref index_ref -id index_id -start index_start -end index_end -d index_dir -f feature_file [-i in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}

	die "-es and -ee should be used together\n" if !$iExonStart ^ !$iExonEnd;
	if (!$iCdsStart ^ !$iCdsEnd) { die "-cs and -ce should be used together\n"; }
	elsif ($iCdsStart && $iCdsEnd && !$iDir) { die "-cs and -ce need -d\n"; }

	$iRef--; $iStart--; $iEnd--;
	$iDir-- if $iDir;
	$iExonStart-- if $iExonStart;
	$iExonEnd-- if $iExonEnd;
	$iCdsStart-- if $iCdsStart;
	$iCdsEnd-- if $iCdsEnd;
	$iInfo-- if $iInfo;
	$allow = 0 if !$allow;
	$verbose = 0 if !$verbose;
}
