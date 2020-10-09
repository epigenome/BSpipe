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

my ($helpFlag, $fixFile, $seqFile, $inFile, $logFile, $outFile, $verbose, $quiet);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"fix=s"			=> \$fixFile,		## input file
	"seq=s"			=> \$seqFile,		## input file
	"input=s"		=> \$inFile,		## input file
	"log=s"			=> \$logFile,
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
) || die "\n";

checkOptions();

my (%forFixs, %bakFixs);
my (%forPoss, %bakPoss);
my %done;
my $noym = 0;

loadFix($fixFile);
loadSeq($seqFile);
process($inFile);
summary();


#-------------------------------------------------------------------------------


#CCGG	+	2
#CCGG	-	3
sub loadFix
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my @a = split;
		my $i = $a[0] eq 'fix' ? 1 : 0;
		if    ($a[$i+1] eq '+') { push(@{$forFixs{$a[$i]}}, $a[$i+2]); }
		elsif ($a[$i+1] eq '-') { push(@{$bakFixs{$a[$i]}}, $a[$i+2]); }
#		else                 { die "Wrong format: $_\n"; }
	}

	close($in) if defined $fileName;
}

sub loadSeq
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput("fasta2seq.pl -head " . ($fileName || '/dev/stdin') . "|");

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my @a = split;

		$forPoss{$a[0]} = {} if !exists $forPoss{$a[0]};
		$bakPoss{$a[0]} = {} if !exists $bakPoss{$a[0]};
		map { findLocation(\$a[1], $_, $forFixs{$_}, $forPoss{$a[0]}) } keys %forFixs;
		map { findLocation(\$a[1], $_, $bakFixs{$_}, $bakPoss{$a[0]}) } keys %bakFixs;
	}

	close($in) if defined $fileName;
}

# ----CCGG----
#        Pos
#      Fix
sub findLocation
{
	my ($rref, $seq, $rpos, $rhash) = @_;

	while ($$rref =~ /$seq/ig)
	{
		#print join(" ", map {substr($$rref, pos($$rref)-length($seq)+$_-2, 4)} @$rpos), "\n";
		#map {my $t = pos($$rref)-length($seq)+$_; print "$t\n" if $t > 358300 && $t < 358400} @$rpos;
		#map {my $t = pos($$rref)-length($seq)+$_; print "$t\n" if $t > 622980 && $t < 630150} @$rpos;
		map {$rhash->{pos($$rref)-length($seq)+$_}++} @$rpos; # 1-based
	}
}

sub process
{
	my ($fileName) = @_;

	print "Processing ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	my $out = openOutput($outFile);

	while (<$in>)
	{
		if (/^\@/)
		{
			print $out $_;
			next;
		}

		my @a = split(/[\t\r\n]/);
		$a[1] = int($a[1]);
		my ($mode) = /ZT:Z:(\w)/;
		die "Missing bisulfite mode in $_" if !$mode;
		die "Wrong bisulfite mode ($mode)in $_" if $mode ne 'C' && $mode ne 'G';
		my $l = length($a[9]);
		my ($lclip, $rclip) = $a[5] =~ /^(?:(\d+)[SH])*.*?(?:(\d+)[SH])*$/;
		my (@poss, @targets);

		if ($mode eq 'C')
		{
			@targets = ($l, $l-1) if !$rclip || $rclip <= 2;
		}
		else
		{
			@targets = (1, 2) if !$lclip || $lclip <= 2;
		}

		if (@targets)
		{
			if (strand($a[1])) # backward
			{
				@poss = check($bakPoss{$a[2]}, $a[9], 'G', $a[3], @targets); # one-based ($l, 2)
			}
			else
			{
				@poss = check($forPoss{$a[2]}, $a[9], 'C', $a[3], @targets); # one-based (1, $l-1)
			}

			if (@poss)
			{
				print "@poss | @a" if $verbose;
				modify2(\@a, \@poss);
				map {$done{$a[2]}{$a[3]+$_-1}++} @poss;
			}
		}

		print $out join("\t", @a), "\n";
	}

	close($in) if defined $fileName;
	warn "$noym reads do not have YM tag\n" if $noym;
}

sub strand
{
	return $_[0] & 0x10 ? 1 : 0;
}

sub check
{
	my ($rhash, $seq, $base, $rpos, @qposs) = @_;
	my @ret;

	foreach my $p (@qposs)
	{
		print "$rpos+$p-1 || $base ne ", substr($seq, $p-1, 1), "\n" if $verbose;
		next if !exists $rhash->{$rpos+$p-1} || $base ne substr($seq, $p-1, 1);
		push(@ret, $p);
	}

	return @ret;
}

sub modify
{
	my($ra, $rp) = @_;
	my ($zmi, $zmv, $ymi);

	for (my $i = 10; $i < @$ra; $i++)
	{
		($zmi, $zmv) = ($i, $1) if $ra->[$i] =~ /^ZM:i:(\d+)/i;
		$ymi = $i if $ra->[$i] =~ /^YM:/i;
	}

	die "No ZM tags in @$ra\n" if !$zmi;
	die "No YM tags in @$ra\n" if !$ymi;

	if ($zmv < @$rp) { die "Unexpected error in @$ra\n"; }
	$ra->[$zmi] = 'ZM:i:' . ($zmv-@$rp);

	if ($zmv == @$rp)
	{
		splice(@$ra, $ymi, 1);
	}
	else
	{
		map {$ra->[$ymi] =~ s/\b$_\b,*//} @$rp;
		$ra->[$ymi] =~ s/,$//;
	}
}

sub modify2
{
	my($ra, $rp) = @_;
	my ($zmi, $zmv, $ymi);
	my ($zui, $zuv, $yui);

	for (my $i = 10; $i < @$ra; $i++)
	{
		if    ($ra->[$i] =~ /^ZM:i:(\d+)/i) { ($zmi, $zmv) = ($i, $1); }
		elsif ($ra->[$i] =~ /^YM:/i) { $ymi = $i; }
		elsif ($ra->[$i] =~ /^ZU:i:(\d+)/i) { ($zui, $zuv) = ($i, $1); }
		elsif ($ra->[$i] =~ /^YU:/i) { $yui = $i; }
	}

	my $total = $zmv + $zuv;
	return if !$total;

	die "No ZM tags in @$ra\n" if !$zmi;
	die "No ZU tags in @$ra\n" if !$zui;
	warn "No YM tags in @$ra\n" if $zmv && !$ymi && $noym++ < 5;
	die "No YM and YU tags in @$ra\n" if !$ymi && !$yui;
	die "Unexpected error for ZM and ZU in @$ra\n" if ($total < @$rp);

	foreach (@$rp)
	{
		$ymi && $ra->[$ymi] =~ s/\b$_\b,*// && $zmv--;
		$yui && $ra->[$yui] =~ s/\b$_\b,*// && $zuv--;
	}

	$ra->[$ymi] =~ s/,$// if $ymi;
	$ra->[$yui] =~ s/,$// if $yui;

	$ra->[$zmi] = "ZM:i:$zmv";
	$ra->[$zui] = "ZU:i:$zuv";

	$ymi && $zmv == 0 && splice(@$ra, $ymi, 1);
	$yui && $zuv == 0 && splice(@$ra, $yui, 1);
	return ($zmv == 0 && $zuv == 0);
}

sub summary
{
	my ($tpos, $tread) = (0, 0);

	my $out = openOutput($logFile ? $logFile : $outFile ? "$outFile.log" : undef);
	print $out "Chr\tPos\tNumRead\n";

	foreach my $id (sort keys %done)
	{
		foreach my $pos (sort {$a<=>$b} keys %{$done{$id}})
		{
			print $out "$id\t$pos\t$done{$id}{$pos}\n";
			$tpos ++;
			$tread += $done{$id}{$pos};
		}
	}

	print $out "Total\t$tpos\t$tread\n";
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

	if ($helpFlag || !$fixFile || (!$seqFile && !$inFile))
	{
		die("Arguments: -f fix_file -s seq_file [[-i] sam_file] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}
}
