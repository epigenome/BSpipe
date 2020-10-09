#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Author: Jeong-Hyeon Choi
## Date: JAN 10 2006

## sort lines by the order of numbers in key field

## Usage:
## CHECK usage with option  -h

## $Id: sortn.pl,v 1.2 2007-07-18 20:51:57 jeochoi Exp $


use Getopt::Long;

my ($helpFlag, $inFile, $outFile, $verbose, $id, $start, $end, $delimeter, $commentFlag, $reverseFlag, $absFlag, $skip, $bedFlag);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"i|inFile:s"	=> \$inFile,		## input file
	"o|outFile:s"	=> \$outFile,		## output file
	"v|verbose+"	=> \$verbose,		## verbose output
	"id:i"			=> \$id,		## to calculate
	"start:i"		=> \$start,		## to calculate
	"end:i"			=> \$end,		## to calculate
	"d|delimeter:s"	=> \$delimeter,		## to calculate
	"c|comment!"	=> \$commentFlag,	## flag for removing comment
	"r|reverse"		=> \$reverseFlag,	## flag for reverse order
	"a|abs"			=> \$absFlag,	## flag for reverse order
	"skip=i"			=> \$skip,
	"bed"				=> \$bedFlag,
);

checkOptions();

my (@records, $out);

load($inFile);
@records = sort
{
	my $aa = $a->[$id] =~ /(\d+)/ ? $absFlag ? abs($1) : $1 : 1232314213;
	my $bb = $b->[$id] =~ /(\d+)/ ? $absFlag ? abs($1) : $1 : 1232314213;
	my $ret = ($aa<=>$bb || $a->[$id] cmp $b->[$id]);
	$ret = $absFlag ? abs($a->[$start]) <=> abs($b->[$start]) : $a->[$start] <=> $b->[$start] if !$ret;
	$ret = $absFlag ? abs($a->[$end  ]) <=> abs($b->[$end  ]) : $a->[$end  ] <=> $b->[$end  ] if defined $end && !$ret;
	$reverseFlag ? -$ret : $ret
}
@records;
output($outFile);


#-------------------------------------------------------------------------------

sub output
{
	my ($fileName) = @_;

	foreach (@records)
	{
		print $out join($delimeter, @$_), "\n";
	}
}

sub load
{
	my ($fileName) = @_;

	my $in = openInput($fileName);
	$out = openOutput($outFile);

	while (<$in>)
	{
		if ($skip-- > 0) { print $out $_; next; }
		next if $commentFlag && /^\s*(#|\/\/)/;

		chop;
		my @a = split(/$delimeter/);
		push(@records, \@a);
	}

	close($in) if $fileName;
}


#-------------------------------------------------------------------------------

sub openInput
{
	my ($fileName) = @_;

	return STDIN unless defined $fileName;

	my ($fd);
	open($fd, $fileName =~ /.gz(ip)?$/ ? "zcat $fileName |" : $fileName =~ /.bz(ip)?2$/ ? "bzcat $fileName |" : $fileName) || die("Open error: $fileName");
	return $fd;
}

sub openOutput
{
	my ($fileName) = @_;

	return STDOUT unless defined $fileName;

	my ($fd);
	open($fd, $fileName =~ /.gz$/ ? "| gzip -c > $fileName" : $fileName =~ /.bz(ip)?2$/ ? "| bzip2 -z -c > $fileName" : ">$fileName") || die("Open error: $fileName");
	return $fd;
}

sub checkOptions
{
	$inFile  = shift(@ARGV) if !defined $inFile  && @ARGV > 0;
	$outFile = shift(@ARGV) if !defined $outFile && @ARGV > 0;

	($id, $start, $end) = 1..3 if $bedFlag;

	if ($helpFlag || !$id || !$start)
	{
		die("Arguments: [-i in_file] [-o out_file] -id column -start column [-end column] [-skip line] [-d delimeter] [-comment | nocomment] [-v]\n"
		  );
	}

	$id--; $start--; $end-- if $end;
	$delimeter = "\t" if !$delimeter;
}
