#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Author: Jeong-Hyeon Choi
## Date: Apr 28 2010

## make a bam file with filtering criteria


## Usage:
## CHECK usage with option  -h

## $Id: perl_template,v 1.2 2007-09-21 22:51:23 jeochoi Exp $

my ($tSingle, $tRead1, $tRead2) = 0..2;

use Getopt::Long qw(:config no_ignore_case);

my ($helpFlag, $inFile, $outFile, $idxFile, $verbose, $quiet);
my ($samIn, $mismatch, $gap, $error, $mapQual, $samOut, $mem, $soft, $sortFlag, $multiFlag);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"si!"			=> \$samIn,		## input file
	"so!"			=> \$samOut,		## input file
	"faidx=s"		=> \$idxFile,
	"mismatch=f"		=> \$mismatch,
	"gap=f"			=> \$gap,
	"error=f"		=> \$error,
	"mq=i"			=> \$mapQual,
	"mem=i"			=> \$mem,
	"soft=f"		=> \$soft,
	"sort"			=> \$sortFlag,
	"multi!"			=> \$multiFlag,
) || die "\n";

checkOptions();

my %rc = ('A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', 'N' => 'N', 'a' => 't', 't' => 'a', 'c' => 'g', 'g' => 'c', 'n' => 'n');
my %cigarIndex = ('I' => 0, 'D' => 1, 'S' => 2, 'H' => 3);

my (@buf, $out, $hasUnique);
load($inFile);

#-------------------------------------------------------------------------------

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $tmpFile = "tmp.$$";
	my $cmd;

	if (!$sortFlag)
	{
		$cmd = $samIn ? $fileName : ("samtools view -h " . ($fileName || '-') . " | ");
	}
	elsif ($fileName ne '-')
	{
		$cmd = "samtools view -bhS " . ($fileName || '-') . ($idxFile ? " -t $idxFile" : '') . ' | ' if $samIn;
		$cmd .= "samtools sort -n -m $mem " . (!$samIn && $fileName ? $fileName : '-') . " $tmpFile";
		run($cmd);

		$cmd = "samtools view -h $tmpFile.bam |";
	}

	my $in = openInput($cmd);
	$out = $samOut ? openOutput($outFile) : openInput("| samtools view " . ($idxFile ? "-t $idxFile " : '') . "-ubS -o - - | samtools sort -m $mem -o $outFile -");

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;

		if (/^\@/)
		{
			print $out $_;
		}
		else
		{
			chop;
			my @b = split(/\t/);
			my @a = (@b[0..10], join("\t", @b[11..$#b]));
			if (@buf && $buf[0][0] ne $a[0])
			{
				onSeq();
				@buf = ();
			}
			push(@buf, \@a);
		}
	}

	onSeq() if @buf;
	close($in) if defined $cmd;
	close($out);

	if (!$hasUnique)
	{
		warn("No unique read\n");
	}
	elsif (!$samOut)
	{
		run("samtools index $outFile");
	}
	unlink("$tmpFile.bam");
}

sub run
{
	my ($cmd) = @_;
	#print "$cmd\n";
	!system($cmd) || die "Error in $cmd\n";
}

sub onSeq
{
	my (@read1, @read2);

	foreach my $a (@buf)
	{
		my ($qid, $flag, $rid, $pos, $mapq, $cigar, $mid, $mpos, $isize, $seq, $qual, $tag) = @$a;
		$flag = int($flag);
		# d=0x400 (duplicate)
		# f=0x200 (failure)
		next if ($flag & 0x400 || $flag & 0x200);

		$mapped = mapped($flag);
#		Wrong assumption
#		die "Wrong condition for mapping: @$a\n" if $mapped != ($cigar ne '*' ? 1 : 0);
		next if !$mapped;

		$flag & 0x80 ? push(@read2, $a) : push(@read1, $a);
	}

	mapStat(@read1);
	mapStat(@read2);
}

sub mapStat
{
	my @reads = @_;

#	print "LOG\t$reads[0][0]\tUNIQUE\t", @reads+0, "\n" if @reads > 1;
	return if $multiFlag ? @reads != 1 : @reads == 0;

	my ($qid, $flag, $rid, $pos, $mapq, $cigar, $mid, $mpos, $isize, $seq, $qual, $tag) = @{$reads[0]};
	$flag = int($flag);
	die "Unexpected flag: $flag\n@{$reads[0]}\n" if ($flag & 0x400 || $flag & 0x200);

#	print "LOG\t$qid\tMQ\t$mapq\n" if $mapq <= $mapQual;
	return if $mapq <= $mapQual;

	my @cigars = parseCigar($cigar);
#	print "LOG\t$qid\tSOFT\t$cigars[2]\n" if (defined $soft && $cigars[2] > ($soft =~ /\./ ? int($soft*length($seq)+0.5) : $soft));
	return if (defined $soft && $cigars[2] > ($soft =~ /\./ ? int($soft*length($seq)+0.5) : $soft));

	my ($nm, $mis, $del);
	die "No NM tag: @{$reads[0]}" if ($tag !~ /\bNM:i:(\d+)/i);
	$nm = $1;
#	print "LOG\t$qid\tERROR\t$nm\n" if (defined $error && $nm > ($error =~ /\./ ? int($error*length($seq)+0.5) : $error));
	return if (defined $error && $nm > ($error =~ /\./ ? int($error*length($seq)+0.5) : $error));

	die "No MD tag: @{$reads[0]}" if ($tag !~ /\bMD:Z:(\S+)/i);
	($mis, $del) = misdel2($1);
	die "Inconsistent number of deletions: $del == $cigars[1]\n@{$reads[0]}\n" if $del != $cigars[1];
	die "Inconsistent number of errors: $nm == $mis+$del+$cigars[0]\n@{$reads[0]}\n" if $nm != $mis+$del+$cigars[0];
#	print "LOG\t$qid\tMISMATCH\t$mismatch\n" if (defined $mismatch && $mis > $mismatch);
#	print "LOG\t$qid\tGAP\t$gap\n" if (defined $gap && $cigars[0]+$cigars[1] > $gap);
	return if (defined $mismatch && $mis > $mismatch) || (defined $gap && $cigars[0]+$cigars[1] > $gap);

	map {print $out join("\t", @$_), "\n"} @reads;
#	print "LOG\t$qid\tPASS\n";
	$hasUnique = 1 if @reads == 1;
}

# 8M1D46M
sub parseCigar
{
	my ($cigar) = @_;
	my @cigars = (0, 0, 0, 0);

	while ($cigar =~ /(\d+)([IDSH])/ig) { $cigars[$cigarIndex{$2}] += $1; }

	return @cigars;
}

sub mapped
{
#	print " $_[0] & 0x4 ? 0 : 1 = ", $_[0] & 0x4 ? 0 : 1, "\n" if $_[0] == 4;
	return $_[0] & 0x4 ? 0 : 1; # mapped => 1
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

# [0-9]+(([ACGTN]|\^[ACGTN]+)[0-9]+)
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

	if ($helpFlag || (!$samOut && !$outFile))
	{
		die("Arguments: [-e number_error] [-m number_mismatch] [-g number_gap] [-mq min_maq_qual] [-mem memory] [-si] [[-i] in_file] [-so] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}

	$error = (defined $mismatch ? $mismatch : 0) + (defined $gap ? $mismatch : 0) if !defined $error && (defined $mismatch || defined $gap);
#	$mismatch = $error if !defined $mismatch;
#	$gap = $error if !defined $gap;
#	die "error should be equal to the sum of mismatches and gaps: mismatch=$mismatch gap=$gap error=$error\n" if $error < $mismatch+$gap;
	$mapQual = 0 if !$mapQual;
	$mem = 2000000000 if !$mem;
	$multiFlag = 1 if !defined $multiFlag;
}
