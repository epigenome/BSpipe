#!/usr/bin/perl -w
#!/bin/perl -w
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Storable qw(dclone);

# option
my ($helpFlag,$iReadFile,$iTmpDir,$iMem,$iBamFile,$iOut);



#perl extreads.pl -r fastq -b bam -o samfile
GetOptions(
        "h|?|help"              => \$helpFlag,
        "r|read=s"               => \$iReadFile,            ## input file 
        "d|tmpdir=s"               => \$iTmpDir,            ## input file 
       # "s|split=i"               => \$iSplitNum,            ## number of split reads 
        "m|mem=f"               => \$iMem,            ## number of split reads 
        "b|bam=s"               => \$iBamFile,            ## input file 
	"o|out=s"		=> \$iOut,		## output file
        # "verbose+"              => \$verbose,           ## verbose output
        # "quiet"                 => \$quiet,
        # "stdev=f"               => \$stdevDiff,
        # "distance=f"    => \$distDiff,
        # "min=i"                 => \$min,
        # "max=i"                 => \$max,
        # "si"                            => \$samIn,
        # "so"                            => \$samOut,
) || die "\n";

checkOptions();
my ($idxID,$idxFlag,$idxRef,$idxPos,$idxMapq,$idxCigar,$idxFd6,$idxFd7,$idxFd8,$idxSeq,$idxQual) = (0..10);


#print STDERR "Loading id from $iReadFile .....\n";

# check split num for memeory 
# 5g -> 4m reads
# globals
my $splitNum = 1300000*$iMem;
#my $splitNum = 10*$iMem;
# for bam file splitting
makeDir($iTmpDir);
my $gfout = openBamOutput($iOut);
my $fbamin = openBamInput($iBamFile);
my $numofbamreads = 0;
my %bamreads;
my $partorder = 0;
while (<$fbamin>){
        ## key:id, value: number of duplication
	chomp();
	# writing header
	if (/^\@/){ 
		print $gfout $_."\n";
		next;
	}
	if (/^(\S+)\t(.*)/){
		my $bamid = trim($1);
		my $line = $bamid."\t".trim($2);
		if (exists($bamreads{$bamid})){
			push @{$bamreads{$bamid}},$line;
		}
		else{
			my @dupreads;
			push(@dupreads,$line);
			$bamreads{$bamid} = \@dupreads;
		}

		if ($numofbamreads != 0 && (($numofbamreads % $splitNum) == 0)|| eof($fbamin)){
			print STDERR "\n\t$numofbamreads reads in bam file have been processed...\n";
			# process bam
			#print "BamReads for part.".$partorder."\n";
			#while ((my $key,my $arr) = each (%bamreads)){
			#	print $key."\t".@{$arr}."\n";
			#}
			splitBamandSort($iReadFile,\%bamreads,$iTmpDir,$partorder);
			$partorder++;
			%bamreads = ();
			
		}	
		$numofbamreads++;
	}
}
close($fbamin);
# merge and sort part.i files;
my @fparts;
my @buflines;
for (my $i=0;$i<$partorder;$i++){
	#print "tmpdir:".$iTmpDir."/part.".$i."\n";
	my $partfilename = $iTmpDir."/part.".$i;
	if (defined ($partfilename)){
		$fparts[$i] = openInput($partfilename);
		$buflines[$i] = "";
	}
	else{
		die "There is no part file\n";
	}
}
while(1){
	my $numoffend = 0;
	for (my $i=0;$i<$partorder;$i++){
		if (eof($fparts[$i]) && $buflines[$i] eq ""){
			$numoffend++;
		}
	}
	last if ($numoffend == $partorder);
	# get the read info from buflines
	for (my $i=0;$i<$partorder;$i++){
		if ($buflines[$i] eq "" && !eof($fparts[$i])){
			my $fp = $fparts[$i];
			defined ($buflines[$i] = <$fp> ) or last;
			chomp($buflines[$i]);
			#print "buflines[$i]:".$buflines[$i]."\n";
		}
	}
	# find min index
	
	my $min ;
	for (my $i=0;$i<$partorder;$i++){
		if ($buflines[$i] ne ""){
			if ($buflines[$i] =~/^(\d+)\t/){
				$min = $1;
			}
		}
	}
	for (my $i=0;$i<$partorder;$i++){
		if ($buflines[$i] ne ""){
			if ($buflines[$i] =~/^(\d+)\t/){
				if ($min >= $1){
					$min = $1;
				}
			}
		}
	}
	for (my $i=0;$i<$partorder;$i++){
		if ($buflines[$i] ne ""){
			if ($buflines[$i] =~/^(\d+)\t(.*)/){
				if ($min == $1){
					print $gfout $2."\n";
					$buflines[$i] = "";
				}
			}
		}
	}
	
	
}

for (my $i=0;$i<$partorder;$i++){
	close($fparts[$i]);
}
close($gfout);
cleanTmpFiles($partorder);

sub splitBamandSort{
	my $ffq = openInput(shift());
	my $rbamreads = shift();
	my $tmpDir = shift();
	my $partorder = shift();
	my $fout = openOutput($tmpDir."/part.".$partorder);


	my @orderedbamreads;
	my @fqorders;
	my $fqidx =0;
	
	for (my $i = 0; my $line = <$ffq>; $i++){
		if ($i%4 == 0) { 
			chomp($line);
			if ($line =~ /^\@(\S+)/){ 
				my $id = $1;
				# global
				if (exists($rbamreads->{$id})){
					push(@orderedbamreads,$rbamreads->{$id});
					push(@fqorders,$fqidx);
				}
				elsif ($id =~ s/\/[12]$//){
					push(@orderedbamreads,$rbamreads->{$id});
					push(@fqorders,$fqidx);
				}
				else{ #process $id /1 or /2
					my $repid1 = $id."/1";
					my $repid2 = $id."/2";
					my $repid = "";
					if (exists($rbamreads->{$repid1})){
						$repid = $repid1;	
					}
					elsif (exists($rbamreads->{$repid2})){
						$repid = $repid2;	
					}
					if ($repid ne ""){
						push(@orderedbamreads,$rbamreads->{$repid});
						push(@fqorders,$fqidx);
					}
				}
				$fqidx++;
			}
		}
	}
	close($ffq);
	for(my $i=0;$i<@orderedbamreads;$i++){
		foreach(@{$orderedbamreads[$i]}){
			print $fout $fqorders[$i]."\t".$_."\n";
		}
	}
	close($fout);
}
sub cleanTmpFiles{
	my $partnum = shift();
	print STDERR "Deleting splitted files\n";
	#print "# of splitted files:".$numpart."\n";
	# delete tmp files and directory
	#unlink glob "$iTmpDir/* $iTmpDir/.*";
	for (my $i=0;$i<$partnum;$i++){
		unlink glob "$iTmpDir/fq.part.".$i;
		unlink glob "$iTmpDir/part.".$i;
	}
	#if ($iTmpDir ne "/tmp"){
		#rmdir $iTmpDir;
	#}
}
sub makeDir{
        my $directory = shift();
        if (!-d $directory) {
		if ($directory ne "/tmp"){
                	mkdir $directory, 0777 or die $directory." is not invalid\n";
		}
        }
}

sub trim {
        my @result = @_;

        foreach (@result) {
                s/^\s+//;
                s/\s+$//;
        }

        return wantarray ? @result : $result[0];
}

sub openBamInput
{
        my ($fileName) = @_;
        return *STDIN unless defined $fileName;

        my ($fd);
        open($fd, $fileName =~ /.bam?$/ ? "samtools view -h $fileName |": $fileName=~/.gz(ip)?$/?" zcat $fileName|": $fileName) || die("Open error: $fileName");
        return $fd;
}
sub openInput
{
        my ($fileName) = @_;

        return *STDIN unless defined $fileName;

        my ($fd);
        open($fd, $fileName =~ /.gz(ip)?$/ ? "zcat $fileName |" : $fileName =~ /.bz(ip)?2$/ ? "bzcat $fileName |" : $fileName) || die("Open error: $fileName");
        return $fd;
}
sub openBamOutput{
        my ($fileName) = @_;

        return *STDOUT unless defined $fileName;

        my ($fd);
        open($fd, $fileName =~ /.bam$/ ? "| samtools view -bS - > $fileName": $fileName =~ /.gz$/ ? "| gzip -c > $fileName": "> $fileName") || die("Open error: $fileName");
        #open($fd,"| samtools calmd -bS - $refFile > $fileName") || die("Open error: $fileName");
        return $fd;
}

sub openOutput
{
        my ($fileName) = @_;

        my ($fd);
        open($fd, $fileName =~ /.gz$/ ? "| gzip -c > $fileName" : $fileName =~ /.bz(ip)?2$/ ? "| bzip2 -z -c > $fileName" : ">$fileName") || die("Open error: $fileName");
        return $fd;
}



sub checkOptions
{
	$iReadFile  = trim(shift(@ARGV)) if !defined $iReadFile  && @ARGV > 0;
	$iMem  = trim(shift(@ARGV)) if !defined $iMem  && @ARGV > 0;
	$iTmpDir  = trim(shift(@ARGV)) if !defined $iTmpDir  && @ARGV > 0;
	$iBamFile  = trim(shift(@ARGV)) if !defined $iBamFile  && @ARGV > 0;
        $iOut = shift(@ARGV) if !defined $iOut && @ARGV > 0;

	#if ($helpFlag ||!$iMem || $iMem <= 0.0|| !$iTmpDir ||!$iReadFile|| !$iBamFile || !$iOut){
	if ($helpFlag ||!$iMem || $iMem <= 0.0|| !$iTmpDir ||!$iReadFile|| !$iBamFile ){
		die("Arguments: [[-m] memory(float giga bytes) > 0.0] [[-d] tmpDir]  [[-r] fastq ]  [[-b] bam or sam input file] [[-o] result sam or bam]"."\t\n");
	}
}

__END__
