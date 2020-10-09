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

my ($helpFlag, $pairFile, $name1, $name2, $inFile, $outFile, $verbose, $quiet);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"p|pair=s"		=> \$pairFile,		## input file
	"n1|name1=s"	=> \$name1,		## input file
	"n2|name2=s"	=> \$name2,		## input file
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
) || die "\n";

checkOptions();

my (@pairs, %sams, %sam2ind, @hs, $col, $format, @buf);

if ($pairFile) { load($pairFile); }
else           { put($name1, $name2); }
process($inFile);


#-------------------------------------------------------------------------------

sub put
{
	my @a = @_;
	$sams{$a[0]}++;
	$sams{$a[1]}++;
	push(@pairs, [$a[0], $a[1], $a[2]]);
}

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my @a = split;
		put(@a);
	}

	close($in) if defined $fileName;
}

sub _format
{
	my (%names, %freq);
	map {$names{$_}++} @_;
	map {$freq{$_}++} values %names;
	die "Wrong format in header: $_" if keys %freq != 1;
	return $names{$_[0]};
}

sub methylation
{
	return if $format == 1;

	my ($ra) = @_;

	if ($format == 3)
	{
		my $cnt = 0;
		map {$cnt++ if $ra->[$_] ne 'NA' && $ra->[$_+1] =~ /^\d+$/ && $ra->[$_+2] =~ /^\d+$/} values %sam2ind;
		$format = $cnt ? 5 : 4;
	}

	foreach my $i (values %sam2ind)
	{
		$ra->[$i] = $format == 2 ? ratio($ra->[$i], $ra->[$i+1]) : $format == 4 ? $ra->[$i+1] : ratio($ra->[$i+1], $ra->[$i+2]);
	}
}

sub ratio
{
	my ($a, $b) = @_;
	return !$b || $b eq 'NA' ? 'NA' : $a/$b;
}

sub process
{
	my ($fileName) = @_;

	print "Processing ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	my $out = openOutput($outFile);

	$_ = <$in>;
	die "No header: $_" if !/^#*Key/;
	@hs = split(/[\t\r\n]/);

	my $i;
	for ($i = 0; $i < @hs; $i++) { last if $hs[$i] !~ /Key\d/; }
	die "Unexpected header: $_" if $i == @hs;
	$col = $i;
	$format = _format(@hs[$col..$#hs]);

	for ( ; $i < @hs; $i += $format)
	{
		warn "$hs[$i] isn't defined\n" if !exists $sams{$hs[$i]};
		$sam2ind{$hs[$i]} = $i;
	}

	print $out join("\t", @hs[0..$col-1], outputPairs(1)), "\n";

	while (<$in>)
	{
		@buf = split;
		methylation(\@buf);
		print $out join("\t", @buf[0..$col-1], outputPairs(0)), "\n";
	}

	close($in) if defined $fileName;
	close($out) if defined $outFile;
}

sub outputPairs
{
	my $flag = shift;
	my @ret;

	foreach (@pairs)
	{
		if (!exists $sam2ind{$_->[0]}) { warn "$_->[0] in a pair is not in the input file\n" if $flag; next; }
		if (!exists $sam2ind{$_->[1]}) { warn "$_->[1] in a pair is not in the input file\n" if $flag; next; }

		push(@ret, $flag ? $_->[2] || "$hs[$sam2ind{$_->[0]}]-$hs[$sam2ind{$_->[1]}]"
		                 : $buf[$sam2ind{$_->[0]}] eq 'NA' ? 'NA' :
		                   $buf[$sam2ind{$_->[1]}] eq 'NA' ? 'NA' :
		                   sprintf("%.3f", $buf[$sam2ind{$_->[0]}]-$buf[$sam2ind{$_->[1]}]));
	}

	return @ret;
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

	if ($helpFlag || !(defined $pairFile ^ (defined $name1 && defined $name2)))
	{
		die("Arguments: <-p pair_file | -p1 name1 -p2 name2> [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}
}
