#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Author: Jeong-Hyeon Choi
## Date: Apr 28 2009

## find the best matches from bowtie outputs for bisulfite reads
## assume that file names includes \.[cg][cg]\.
## where the first represents the coversion for references
##       the second represents the coversion for reads
##       c means all Ts are coverted to Cs
##       g means all As are coverted to Gs
## so cc and cg allow c-t mismatch
##    gc and gg allow g-a mismatch
## also cc and gg have reads mapped to forward strand
##      cg and gc have reads mapped to reverse strand

my ($iQue, $iFlag, $iRef, $iPos, $iMapq, $iCigar, $iSeq, $iQual) = (0..5, 9, 10);

use Getopt::Long qw(:config no_ignore_case);

my ($helpFlag, @inFiles, $outFile, $verbose, $quiet, $refFile, $seqFile, $errorCutoff, $uniqueFlag);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \@inFiles,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"ref=s"			=> \$refFile,
	"seq=s"			=> \$seqFile,
	"error=i"		=> \$errorCutoff,
	"unique"			=> \$uniqueFlag,
) || die "\n";

checkOptions();

my (@ins, $out);
my ($line, @buf, @lookup);
my (%refs, %reads, @modes);

loadFasta($refFile);
my $rin = openFastq($seqFile);
my $rno = 0;

foreach my $fileName (@inFiles)
{
	my $in = openInput($fileName =~ /\.gz$/ ? "zcat $fileName | grep -v '^\@' | " : "grep -v '^\@' $fileName | ");
	die "Wrong file name: $fileName\n" if $fileName !~ /\.([cg][cg])\.[bs]am(\.gz)*$/;
	push(@modes, $1);
	push(@ins, $in);
}

$out = openOutput($outFile);
writeHeader();
process();


#-------------------------------------------------------------------------------

sub writeHeader
{
	my $in = openInput($inFiles[0]);

	while (<$in>)
	{
		last if !/^\@/;
		print $out $_;
	}

	close($in);
}

sub loadFasta
{
	my ($fileName) = @_;
	print "Loading $fileName ...\n" if !$quiet;
	my $in = openInput($fileName);
	my ($id, $seq);

	while (<$in>)
	{
		if (/^>\s*(\S+)/)
		{
			$refs{$id} = $seq if $id;
			($id, $seq) = ($1, '');
		}
		else
		{
			chop;
			$seq .= uc $_;
		}
	}

	$refs{$id} = $seq if $id;
	close($in);
}

sub openFastq
{
	my ($fileName) = @_;
	print "Loading $fileName ...\n" if !$quiet;
	return openInput($fileName);
}

sub loadFastq
{
	my ($match) = @_;
	my ($id, $seq);

	# Remove /[12] in read IDs in the FASTQ file because even paired reads are mapped in single read mode.
#	my $rid = "$match->[$iQue]/" . readType($match->[$iFlag]);
	my $rid = "$match->[$iQue]";
	my $flag = $rid =~ /\/[12]$/ ? 0 : 1; # read IDs in bam has /1 or /2?
	return @{$reads{$match->[-3] = $match->[$iQue]}} if exists $reads{$match->[$iQue]};
	return @{$reads{$match->[-3] = $rid           }} if exists $reads{$rid           };
	return @{$reads{$match->[-3]                  }} if exists $reads{$match->[-3]   };

	for (my $i = 0; <$rin>; $i++)
	{
		chop;
		if    ($i%4 == 0) { ($id) = /^\@(\S+)/; $id =~ s/\/[12]$// if $flag; }
		elsif ($i%4 == 1) { $seq = $_; }
		elsif ($i%4 == 2) { $reads{$id} = [++$rno, $seq]; }
		else              { return ($rno, $seq) if $id eq $rid || $id eq $match->[$iQue] || $id eq $match->[-3]; }
	}
}

sub process
{
	print "Processing ", join(' ', @inFiles), " ...\n" if !$quiet;

	while (1)
	{
		my $flag;

		for (my $i = 0; $i < @inFiles; $i++)
		{
			loadLine($i) if !$buf[$i];
			$flag = 1 if $buf[$i];
		}

		last if !$flag;
#		sameReads();
		onRead();
		print "==================\n" if $verbose;
	}
}

sub loadLine
{
	my ($index) = @_;
	my $in = $ins[$index];

	while (($_ = $lookup[$index] || <$in>))
	{
		next if /^\@/;
		$lookup[$index] = undef if $lookup[$index];
		chomp;

		my @a = split(/\t/);
		$a[$iFlag] = int($a[$iFlag]);
		next if isUnmapped($a[$iFlag]) || strand($a[$iFlag]) != bsStrand($modes[$index]);
		my @b = (@a[0..$iQual], join("\t", @a[$iQual+1..$#a]));
		push(@b, getReadId(@b[$iQue, $iFlag]), undef, undef); # add the elements for read id and #error

		last if ($buf[$index] && $buf[$index][0][-3] ne $b[-3]);

		push(@{$buf[$index]}, \@b);
	}

	$lookup[$index] = $_ if $_;
}

sub isUnmapped
{
	return $_[0] & 0x4 || $_[0] & 0x400 || $_[0] & 0x200 ? 1 : 0;
}

sub getReadId
{
	my ($id, $flag) = @_;
	$id =~ s/\s+.*$//;
	return $flag & 0x40 ? "$id/1" : $flag & 0x80 ? "$id/2" : $id;
}

sub strand
{
	return $_[0] & 0x10 ? 1 : 0;
}

sub sameReads
{
	my $id;

	foreach my $i (0..$#buf)
	{
		foreach my $a (@{$buf[$i]})
		{
			if (!$id) { $id = $a->[-2]; }
			elsif ($a->[-2] ne $id) { die "Different read ID $a->[-2] at file $i to $id\n"; }
		}
	}
}

sub onRead
{
	my ($rorder, $rid, @bests);

	foreach my $i (0..$#buf)
	{
		foreach my $a (@{$buf[$i]})
		{
			print join("\t", @$a[0..$#$a-3]), "\n" if $verbose;
			onMatch(bsMode($modes[$i]), $a) if !defined $a->[-1];
			@bests = () if @bests && ($rorder > $a->[-2] || ($rorder == $a->[-2] && $bests[0][-1] > $a->[-1])); # read order and #errors
			($rid, $rorder) = ($a->[-3], $a->[-2]) if !defined $rorder || $rorder > $a->[-2];
			push(@bests, $a) if !@bests || ($rorder == $a->[-2] && $bests[0][-1] == $a->[-1]);
		}
	}

	if (@bests)
	{
		map {print "M: @$_\n"} @bests if $verbose;
		print "---------------\n" if $verbose;

		if ((!defined $errorCutoff || $errorCutoff >= $bests[0][-1]) && (!defined $uniqueFlag || @bests == 1))
		{
			map {print $out join("\t", @$_[0..$#$_-3]), "\n"} @bests;
		}
	}

	if (defined $rorder)
	{
		print "$rorder\n" if $verbose;
		delete $reads{$rid};
		map {$buf[$_] = undef if $buf[$_] && (!defined $buf[$_][0][-1] || $buf[$_][0][-2] == $rorder)} (0..$#buf);
	}
	else
	{
		@buf = ();
	}
}

# $mode == 0 : CC
# $mode == 1 : CG
# $mode == 2 : GC
# $mode == 3 : GG
# HWI-EAS407_0023:2:108:1507:1216#0/1     0       chrY    27194928        255     51M     *       0       0       CGGGGTGGGGTGGGGTGGGGTGGGGTGGGGCGGGGTGGGGTGGGGTGGGGT     BBCCCCCCBCCCBACCCCBCCCABCCCCBCCCCBBCC?BBBBBABBBABBB XA:i:0  MD:Z:30T1T2G1C13        NM:i:4
# add the new tags and set read number and the number of error to the last two columns
sub onMatch
{
	my ($mode, $match) = @_;
#	my ($aliLen) = $match->[$iCigar] =~ /(\d+)[MDI]/;
#	die "Wrong format for alignment length: @$match\n" if !$aliLen;
	my ($rorder, $ss) = loadFastq($match);
	die "Unknown read sequence: $match->[$iQue]\n" if !$ss;
	$ss = rc($ss) if strand($match->[$iFlag]);
	$match->[-2] = $rorder;
	die "Unknown ref sequence: $match->[$iRef]\n" if !exists $refs{$match->[$iRef]};
	my ($error, $clip, $skip, $md, $rmethy, $runmethy) = align($mode, $match->[$iCigar], \$refs{$match->[$iRef]}, $match->[$iPos]-1, \$ss);
	$match->[$iSeq] = $ss;
	$match->[-4] =~ s/\b(MD:Z:)/$1$md\tZD:Z:/i;
	$match->[-4] =~ s/\b(NM:i:)/$1$error\tZE:i:/i;
	$match->[-4] .= "\tZT:Z:" . $mode;
	$match->[-4] .= "\tZM:i:" . scalar(@$rmethy  );
	$match->[-4] .= "\tZU:i:" . scalar(@$runmethy);
	$match->[-4] .= "\tYM:Z:" . join(",", @$rmethy  ) if @$rmethy;
	$match->[-4] .= "\tYU:Z:" . join(",", @$runmethy) if @$runmethy;
	$match->[-1] = $error+$clip;
}

sub splitCigar
{
	my ($s) = @_;
	my $chars = "MDISHN=X";
	my @t;
	die "Unsupported cigar character $& in $s\n" if $s =~ /[^$chars\d]/;
	while ($s =~ /(\d*)([$chars])/g) { push(@t, $1, $2); }
	for (my $i = 2; $i < @t; $i += 2) { $t[$i] += $t[$i-2]; }
	print "(@t)\n" if $verbose;
	return @t;
}

sub align
{
	my ($mode, $cigar, $rref, $pos, $rseq) = @_;
	my @cigars = splitCigar($cigar);

	print "$mode  R:", substr($$rref, $pos, $cigars[-2]), "\n   S:$$rseq\n" if $verbose;
	my ($mat, $cnt, $error, $clip, $skip, $md) = (0, 0, 0, 0, 0, '');
	my (@methy, @unmethy);
	my ($ir, $is, $il, $j) = (0, 0, 0, 0);

	for ( ; $il < $cigars[-2]; $il++)
	{
		print "   ($il $ir $is) $j ($mat $cnt $error) $md\n" if $verbose;

		if ($il >= $cigars[$j])
		{
			$j += 2;
			if ($cigars[$j+1] eq 'D') { $md .= ($mat+$cnt) . '^'; ($mat, $cnt) = (0, 0); }
		}

		if ($cigars[$j+1] eq 'I') { $is++; $error++; next; }
		if ($cigars[$j+1] eq 'D') { $md .= substr($$rref, $ir+$pos, 1); $ir++; $error++; next; }
		if ($cigars[$j+1] eq 'S') { $is++; $clip++; next; }
		if ($cigars[$j+1] eq 'H') { $is++; $clip++; next; }
		if ($cigars[$j+1] eq 'N') { $ir++; $skip++; next; }

		my $r = substr($$rref, $ir+$pos, 1);
		my $s = substr($$rseq, $is, 1);

		if ($r eq $s)
		{
			$mat++;
			push(@methy, $ir+1) if ($r eq $mode);
		}
		elsif (bisulfite($mode, $r) eq $s)
		{
			$cnt++;
			push(@unmethy, $ir+1);
		}
		else
		{
			$error++;
			$md .= ($mat+$cnt) . $r;
			($mat, $cnt) = (0, 0);
		}

		$is++; $ir++;
	}

	$md .= $mat+$cnt;
	if ($cigars[1] eq 'H' || $cigars[-1] eq 'H')
	{
		my ($lclip, $rclip) = ($cigars[1] eq 'H' ? $cigars[0] : 0, $cigars[-1] eq 'H' ? $cigars[-2]-$cigars[-4] : 0);
		$$rseq = substr($$rseq, $lclip, length($$rseq)-$lclip-$rclip);
	}
	print "$error $clip $skip $md (@methy) (@unmethy)\n" if $verbose;
	return ($error, $clip, $skip, $md, \@methy, \@unmethy);
}

sub bisulfite
{
	my ($mode, $base) = @_;
	return $mode eq 'C' ? $base eq 'C' ? 'T' : $base : $base eq 'G' ? 'A' : $base;
}

sub readType
{
	return $_[0] & 0x80 ? '2' : '1';
}

sub rc
{
	my ($seq) = @_;
	my $rc = reverse $seq;
	$rc =~ y/ACGTacgt/TGCAtgca/;
	return $rc;
}

sub bsMode
{
	return uc(substr($_[0], 0, 1));
}

sub bsStrand
{
	return uc(substr($_[0], 0, 1)) eq uc(substr($_[0], 1, 1)) ? 0 : 1;
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
	push(@inFiles, @ARGV) if !@inFiles  && @ARGV > 0;

	if ($helpFlag || !@inFiles || !$refFile || !$seqFile)
	{
		die("Arguments: -r ref_fasta -s read_fastq [-i] in_file [-o out_file] [-v] [-q]\n"
		  . "\t-error NNN : output only reads with N or less errors\n"
		  . "\t-unique : output only uniquely mapped reads\n"
		  );
	}
}
