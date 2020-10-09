#!/usr/bin/perl -w
#!/bin/perl -w

my ($end5, $end3, $body) = qw(promoter downstream body);

use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Storable qw(dclone);

# option
my ($helpFlag,$iInput,$iOut);

GetOptions(
        "h|?|help"              => \$helpFlag,
        "i|input=s"               => \$iInput,            ## input file 
        "o|out=s"               => \$iOut,            ## output file 
) || die "\n";

checkOptions();

my ($path, $name, $in, $out);

foreach (<$iInput.*.ann.stat>)
{
	$path = $_;
	$name = $path;
	$name =~ s/$iInput.//;
	$name =~ s/\.ann\.stat//;

	print ">$path\n";
	$in = openInput($path);
	my $s = <$in>;
	seek($in, 0, 0) || die "Rewind error\n";;

	$out = openOutput($iOut) if !$out;

	if ($s =~ /^\d+all/)
	{
		geneSummary()
	}
	else #if ($s =~ /^## (INFO|SCORES)/)
	{
		featureSummary();
	}

	close($in);
}

=refseq
priority    both_5'end      0
priority    5'end   489
priority    both_3'end      11
priority    3'end   218
priority    5'UTR   140
priority    3'UTR   107
priority    CDS     736
priority    exon    0
priority    intron  2040
priority    intergenic      2009
priority    none    0
=cut
sub geneSummary
{
	my %refseqTitle = ( "upstream" => 'promoter', "downstream" => 'downstream', 'CDS' => 'exon');
	my @refseqTitle = ('promoter', "5'UTR", "exon", "3'UTR", "intron", "downstream", "intergenic");
	my (%refseq, %strs);

	while (<$in>){
		chomp();
		my @values = split("\t",$_);
		next if (($values[0] ne 'all' && $values[0] ne "priority" && $values[0] ne 'combination') || $values[1] eq 'none');
		$strs{$values[0]} .= $values[1]."\t".$values[2]."\n";
		$values[1] =~ s/bi_//;
		$values[1] =~ s/([53])UTR/$1'UTR/g;
		$values[1] = join(',', map {exists $refseqTitle{$_} ? $refseqTitle{$_} : $_} split(/,/, $values[1]));
		$values[1] =~ s/exon,// if $values[1] =~ /exon.+exon/;
		$refseq{$values[0]}{$values[1]} += $values[2];
	}

	foreach my $id (sort keys %refseq) {
		my $_path = $path;
		print $out "######\t$name.$id\t$path\n";
		print $out "$strs{$id}\n";
		print $out "pie chart for gene summary ($id)\n";
		map {print $out "$_\t", $refseq{$id}{$_}||0, "\n"} $id ne 'combination' ?  @refseqTitle : sort {$refseq{$id}{$b} <=> $refseq{$id}{$a}} keys %{$refseq{$id}};
		print $out "\n";
	}
}

sub featureSummary
{
	my $last_line = "";
	my ($flag, @feats);

	while (<$in>){
		chomp();
		$last_line = $_ if /^total\s/;

		if ($flag)
		{
			if (/^#/) { $flag = 0; }
			else      { push(@feats, [split(/\t/, $_)]); }
		}
		elsif (/^## FEATURES/)
		{
			<$in>;
			$flag = 1;
		}
	}

	print $out "######\t$name\t$path\n";
	classSummary(\@feats) if @feats < 20;

	my @cgivalues= split("\t", $last_line);
	if ($cgivalues[0] eq "total"){
		print $out "upstream\t".$cgivalues[1]."\n";
		print $out "body\t".$cgivalues[2]."\n";
		print $out "downstream\t".$cgivalues[3]."\n";
		print $out "outside\t".$cgivalues[4]."\n";

		print $out "\npie chart\n";
		print $out "shore\t".($cgivalues[1]+$cgivalues[3])."\n";
		print $out "body\t".($cgivalues[2])."\n";
		print $out "outside\t".$cgivalues[4]."\n";
	}
	else{
		die("Error in $path\n");
	}
	print $out "\n";
}

sub classSummary
{
	my ($ra) = @_;

	print $out "feature plot\n";
	foreach my $rb (@$ra) {
		print $out join("\t", $rb->[0], $rb->[2]+$rb->[3]+$rb->[4]), "\n";
	}

	print $out "\n";
}

### FEATURES
##       5'out   5'end   body    3'end   3'out
#DNA     2       0       1       0       4       7
#LINE    2       0       3       0       9       14
#LTR     3       0       3       0       4       10
#Low_complexity  5       0       5       0       12      22
#SINE    27      0       17      0       10      54
#Satellite       0       0       0       0       1       1
#Simple_repeat   9       0       3       0       13      25
#scRNA   0       0       2       0       0       2
##       5'end   body    3'end   intergenic      none
#total   0       34      0       101     0

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
        open($fd, $fileName =~ /.gz$/ ? "| gzip -c > $fileName" : $fileName =~ /.bz(ip)?2$/ ? "| bzip2 -z -c > $fileName" : ">$fileName") ||
 die("Open error: $fileName");
        return $fd;
}


sub checkOptions
{
	$iInput = shift(@ARGV) if !defined $iInput && @ARGV > 0;
	$iOut = shift(@ARGV) if !defined $iOut && @ARGV > 0;

	if ($helpFlag || !$iInput) {
		die("Arguments:  [-i input prefix till before .refgene.ann.stat or .cgi.ann.stat ex) annot/res/mouse.win.fisher.cis3d-sham48h.cis3d-sham48h.fdr0.05diff0.25.xls.hypo ] [-o output files]"."\t\n");
	}
}

