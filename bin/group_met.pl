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

my ($helpFlag, $grpFile, $inFile, $outFile, $verbose, $quiet, $confFlag);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"g|group=s"		=> \$grpFile,		## input file
	"input=s"		=> \$inFile,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"conf"			=> \$confFlag,
) || die "\n";

checkOptions();

my (@groups, %sam2ind, @hs, $col, $format);

load($grpFile);
process($inFile);


#-------------------------------------------------------------------------------

sub load
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	my %done;

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my @a = split;
		die "A group file should have at least 2 or 4 columns\n" if $confFlag ? @a != 4 : @a < 2;
		my $sam = $a[$confFlag ? 2 : 0];
		foreach my $g (@a[$confFlag ? (3) : (1..$#a)])
		{
			next if exists $done{$g}{$sam}; # ignore duplicates
			push(@{$groups{$g}}, $sam);
			$done{$g}{$sam}++;
		}
		$sam2ind{$sam} = -1;
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
	my ($rval, $rnam) = @_;
	my ($num, @sum) = (0);
	map {push(@sum, 0)} 1..($format >= 3 ? 3 : $format);

	if ($format == 3)
	{
		my $cnt = 0;
		map {$cnt++ if $rval->[$_+2] >= 1 && $rval->[$_+2] !~ /^\d+\./} values %sam2ind;
		$format = $cnt ? 5 : 4;
	}

	foreach my $n (@$rnam)
	{
		next if $format == 2 ? $rval->[$sam2ind{$n}+1] == 0 : $rval->[$sam2ind{$n}] eq 'NA';
		$num++;

		if    ($format < 5                                 ) { $sum[0] += $rval->[$sam2ind{$n}]; }
		elsif ($num == 1 || $sum[0] < $rval->[$sam2ind{$n}]) { $sum[0]  = $rval->[$sam2ind{$n}]; }
		$sum[1] += $rval->[$sam2ind{$n}+1] if $format > 1;
		$sum[2] += $rval->[$sam2ind{$n}+2] if $format > 3;
	}

	if    ($format == 1) {      $sum[0 ] = !$num ? 'NA' : sprintf("%.3f", $sum[0 ]/$num)      ; }
	elsif ($format == 4) { map {$sum[$_] = !$num ? 'NA' : sprintf("%.3f", $sum[$_]/$num)} 0..2; }

	return @sum;
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
		warn "$hs[$i] isn't defined\n" if !exists $sam2ind{$hs[$i]};
		$sam2ind{$hs[$i]} = $i;
	}

	my @groups = sort keys %groups;
	print $out join("\t", @hs[0..$col-1], map {my $v = $_; map {$v} 1..$format} @groups), "\n";

	while (<$in>)
	{
		my @buf = split;
		print $out join("\t", @buf[0..$col-1], map {methylation(\@buf, $groups{$_})} @groups), "\n";
	}

	close($in) if defined $fileName;
	close($out) if defined $outFile;
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

	if ($helpFlag || (!$grpFile && !$inFile))
	{
		die("Arguments: [-conf] -g group_file [[-i] in_file] [[-o] out_file] [-v] [-q]\n"
		  . "\t\n"
		  );
	}
}
