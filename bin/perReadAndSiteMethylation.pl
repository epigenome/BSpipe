#!/usr/bin/perl

use strict;
use Getopt::Long;
use Data::Dumper;

my $qual;
my $minQual;
my $out;
my ($leftOff, $rightOff) = (100, -1);

my %IUPAC=(
			"A"=>"A",
			"T"=>"T",
			"C"=>"C",
			"G"=>"G",
			"M"=>"AC",
			"R"=>"AG",
			"Y"=>"CT",
			"W"=>"AT",
			"S"=>"CG",
			"K"=>"GT",
			"B"=>"CGT",
			"D"=>"AGT",
			"H"=>"ACT",
			"V"=>"ACG",
			"N"=>"ATCG"
);
			
my %trans=("A"=>"T",
    "T"=>"A",
    "G"=>"C",
    "C"=>"G",
    "N"=>"N"
);

my $samfile="";
my $output="";
my $seqfile="";

my ($incFile, %posIncluded);
my ($excFile, %posExcluded);
my $target_bases;
my %target_info;
my $inputProgram="";

my %result;
my %buf; # for final result
my %refseqMeth;
my $startSeq=0;
my $program="";
my $readId="";

my $preChrom="";
my $preMatchEnd="";

my $in_sam;
my $in_ref;
my %ref_index;
my $readseq_Fun;

sub checkBasFilter{
	my ($iupacStr)=@_;
	
	my @iupacS=split(//, $iupacStr);
	foreach (@iupacS){
		if(! exists $IUPAC{$_}){
			return $_;
		}
	}
	
	return -1;
}			

sub basFilter{
	my ($seq,$iupacStr)=@_;
	my @seqs=split(//,uc($seq));
	my @iupacS=split(//,$iupacStr);
	
	my $length=0;
	if( length($seq)<length($iupacStr) ){
		$length=length($seq);
	}else {
		$length=length($iupacStr);
	}
	
	for(my $i=0; $i<$length; ++$i){
		#print $seqs[$i]." ".$IUPAC{$iupacS[$i]}."\n";
		if(index($IUPAC{$iupacS[$i]},$seqs[$i])<0){
			return 0;
		}
	}
	#print '1'."\n";
	return 1;
	
}
sub DNAreverse {
	my $seq=shift;
 	my @char=split(//,$seq );
	my @reverse=reverse(@char);
	
	my $re="";

	foreach my $element (@reverse){
		if(exists $trans{$element}){
		 	$re=$re.$trans{$element};
		}else {
			#print $element."\n";
		}
	}
	return $re;
}

sub getTag{
	my($line_array,$start)=@_;
	
	my %tags;
	
	for(my $i=$start ; $i<scalar(@$line_array) ; $i++ ){
		my $tag=$$line_array[$i];
		my @splits=split(/:/,$tag);
		
		#print Dumper(\@splits);
		
		if(scalar(@splits) == 3){
			$tags{$splits[0]}=$splits[2];
		}
		
	}
	
	return %tags;
}

sub getStand{
	my $seqname="";
	my $sequence_ref="";
	
	foreach  (keys %refseqMeth){
			$seqname=$_;
			$sequence_ref=&$readseq_Fun($in_ref,\%ref_index,$seqname);
		
			my $position_hash=$refseqMeth{$seqname};
			my @positions=keys %$position_hash;
		
			my ($chrome,$start,$end)=parseSeqname($seqname);
			
			#print "$seqname .. $sequence \n";
		
			foreach my $position (@positions) {
				my ($methylatedN,$downstream)=getMdownstream($sequence_ref,$position-$start+1);
				if(exists $result{$chrome}{$position}){
					$result{$chrome}{$position}[3] = $methylatedN;
					$result{$chrome}{$position}[4] = $downstream;
				}
			}
	}
	
	#my $i=0;

	#while(my $line=<FASTA>){
	#	chomp $line;
	
		#print $line."\n";
	#	if(index($line,">")>-1){
	#		if(length($seqname)>1){
				
				#printn "$seqname";
				#printn "$sequence";
	#			++$i;
				
	#			fillStandFollow($seqname,$sequence,$result_hash,$refseqMeth_hash);
				
	#		}
	#		$seqname=$line;
	#		$seqname=~s/\>//;
		
	#		$sequence="";
	#		next;		
	#	}
	#	$line=~s/\s//;
	#	$sequence=$sequence.$line;
	#}
	#	fillStandFollow($seqname,$sequence,$result_hash,$refseqMeth_hash);
}

sub parseSeqname{
	my $seqname=$_[0];
	
	my ($chrom,$start,$end)=split(/[:,\.]+/,$seqname);
	$start=1  unless (defined  $start); 
	return ($chrom,$start,$end);
}

sub fillStandFollow{
	if(scalar(@_)<4){
		return ;
	}
	my ($seqname,$sequence,$result_hash,$refseqMeth_hash)=@_;
	
	if(exists $$refseqMeth_hash{$seqname}){
		my $position_hash=$$refseqMeth_hash{$seqname};
		my @positions=keys %$position_hash;
		
		my ($chrome,$start,$end)=parseSeqname($seqname);
		
		foreach my $position (@positions) {
			if(exists $$result_hash{$chrome}{$position}){
				my ($methylatedN,$downstream)=getMdownstream($sequence,$position-$start+1);
				$$result_hash{$chrome}{$position}[3]=$methylatedN;
				$$result_hash{$chrome}{$position}[4]=$downstream;
			}
		}
	}else{
		return ;
	}
}

sub getMdownstream{
	my ($refseq_ref,$position)=@_;

	my $methylatedN = uc(substr($$refseq_ref,$position-1,1));
	my $downstream  = $methylatedN eq 'C' ? uc(substr($$refseq_ref, $position-1+$leftOff,$rightOff-$leftOff+1)) :
	                             DNAreverse(uc(substr($$refseq_ref, $position-1-$rightOff,$rightOff-$leftOff+1)));

	while ($position++-1+$leftOff < 0)
	{
		$downstream = "-$downstream";
	}

	return ($methylatedN,$downstream);
}

sub getPgname{
	my ($line)=@_;
	my @temp=split(/\t/,$line);
	
	shift @temp;
	
	foreach (@temp){
		my @headtag=split(/:/,$_);
		if( $headtag[0] eq "ID"){
			return $headtag[1];
		}
	}
	
	return "";
}

sub readFilter{
	my ($file,$filter_hash)=@_;
	
	open(FILTER,$file) || die "can't open file $file.";
	
	while(<FILTER>){
		my ($chrome,$pos)=split(/\t/,$_);
		$$filter_hash{$chrome}{$pos}=0;
	}

	close(FILTER);
}

sub output{
	my ($chrome, $position, $offset, $length,$pattern,$updown,)=@_;
	my @ret;
	
	my	$methy=$result{$chrome}{$position}[0];
	my	$unmethy=$result{$chrome}{$position}[1];
	my	$reads=$result{$chrome}{$position}[2];
	my	$ncbase=$result{$chrome}{$position}[3];
	my $strand= $ncbase eq "C" ? "+" : "-" ;
	my $mQual=$result{$chrome}{$position}[5];
	my $uQual=$result{$chrome}{$position}[6];
	
	my $position2=$position;
	if($position =~  /^\d+$/){
		if($strand eq "+"){
			$position=$position+$offset;
			$position2=$position+$length-1;
		}else {
			$position2=$position-$offset;
			$position=$position2-$length+1;
		}
	}else {
		$position2=$position
	}
	
	$pattern= $pattern ne ""?$pattern:"$updown";
		
	my $score=0;
	if($qual){
		if($mQual+$uQual>0){
			$score=$mQual/($mQual+$uQual);
		}else {
			$score=-1;
		}
	}else{
		$score=$methy/($methy+$unmethy);
	}
		
	$score=sprintf("%.3f",$score);
	
	push(@ret,
		$chrome,
 		$readId,
	 	$pattern,
 		$position,
 		$position2,
 		$score,
 		$strand,
 		$methy,
		$unmethy);
 	if($qual){
 		push(@ret, $mQual, $uQual);
 	}
 	push(@ret, $updown);

	push(@{$buf{$chrome}{$position}{$strand}}, \@ret);
}

sub getQual{
	
	my ($quality, $pos)=@_;
	
	if(length($quality)==1){
		return 0;
	}
	my $pQual=substr $quality,$pos-1,1;
	return ord($pQual)-33;
}

sub printHelp{
	print " msam2gff.pl convert sam format to gff format.\n";
	print " msam2gff.pl -sam SamFile -o OutputFile -seq referenceFile.fasta ";
	print "[-qual] [-mq min_base_qual] [-inc includeFile] -f IUPACs[:IUPACs] [-p ProgramName]\n";
	#print "-n how many sites will be loaded in memory";
}


#sub outputTempFile{ 
#	my (result,$posfilter_hash,$tempfile)=@_;
#	
#	open(TEMP,">$tempfile") || die "can't open tempfile $tempfile";
#			outputFile(result,$posfilter_hash,*TEMP,1);
#			
#			 close TEMP;
#}

sub isOverlap{
	my ($chrom,$start,$preChrom,$preEnd)=@_;
	
	if( $chrom eq $preChrom && $preEnd>=$start){
		 return 1;
	}else {
		return 0;
	}
}

sub makeIndex{
	my  $file=shift;
	my $index_hash=shift;
	my $preoffset=0;
	my $seqname="";
	my $seqstart=0;
	my $seqlength=0;

	while(<$file>){
		chomp $_;
		if(/^>/){
			if(length $seqname > 1){
				$$index_hash{$seqname}=[$seqstart,$preoffset,$seqlength];
			}
			$seqname = substr($_,1);
			$seqstart=tell($file);
			$seqlength=0;
		}
		$seqlength+=length($_);
		$preoffset=tell($file);
	}
	if(! exists $$index_hash{$seqname}){
		$$index_hash{$seqname}=[$seqstart,$preoffset,$seqlength];
	}
}

sub getSequenceLength{
	my $index_hash=shift;
	my $seqname=shift;
	
	if(exists $$index_hash{$seqname}){
		return $$index_hash{$seqname}->[2];
	}
	return 0;
}

sub creadsequence{
	my $cSeq="";
	my $cSeqname="";
	my $rs=sub{
		my ($file, $index_hash,$seqname)=@_;
		if(!defined $seqname){
			return "";
		}
		if($cSeqname eq $seqname){
			return \$cSeq;
		}else {
			$cSeqname=$seqname;
			$cSeq=readsequence($file,$index_hash,$seqname);
			return \$cSeq;
		}
		
	};
	return $rs;
}

sub readsequence{
	my $file=shift;
	my $index_hash=shift;
	my $seqname=shift;
	
	if(!exists $$index_hash{$seqname}){
		return "";
	}
	
	my $start=$$index_hash{$seqname}->[0];
	my $end  =$$index_hash{$seqname}->[1];
	
	seek($file,$start,0);

	my $sequence="";

	read ($file,$sequence ,$end-$start);
	#print "readfile\n";
	$sequence=~ s/\n//g;
	$sequence=~ s/(.)/\U$1\E/g;
	return $sequence;
}

Getopt::Long::GetOptions(              
      'sam=s'    => \$samfile,            
       'o=s'    => \$output,           
       "seq=s"=>\$seqfile,          
       'include=s'    => \$incFile,           
       'exclude=s'    => \$excFile,           
       'n=s' => \$target_bases,
       'p=s' => \$inputProgram,
       'mq' => \$minQual,
       'qual' => \$qual); 

if (!$samfile && !$seqfile){
	printHelp();
	exit(1);
}

if($incFile){ readFilter($incFile,\%posIncluded); }
if($excFile){ readFilter($excFile,\%posExcluded); }

$target_bases = 'CG:CH' if !$target_bases;
if($target_bases ne ""){
	my @temp=split(/:/, uc($target_bases));
	my $check=checkBasFilter(join("",@temp));
	if($check!=-1){
		die "$target_bases: $check not in IUPAC.";
	}
	
	my @indexes;
	foreach (@temp){
		my $ob = $_;
		die "Bases $ob doesn't have C\n" if $ob !~ /^([^C]*)C(.*)$/;
		my $i = -length($1);
		my $j =  length($2);
		push(@indexes, $i);
		$leftOff = $i if $leftOff > $i;
		$rightOff = $j if $rightOff < $j;
	}
	print STDERR "Offset = ($leftOff $rightOff)\n";

	for (my $i = 0; $i < @temp; $i++){
		my $ob = $temp[$i];
		my @bs;
		map {push(@bs, $IUPAC{$_})} split(//, $ob);
		#print "$ob: (@bs)\n";

		my %comp;
		foreach my $c (@bs)
		{
			my @ks = keys %comp;
			map {$comp{$_} = $ob} split(//, $c) if !@ks;
			#print "\t$c: (@ks)\n";

			foreach my $k (@ks)
			{
				delete $comp{$k};
				map {$comp{"$k$_"} = $ob} split(//, $c);
			}
			#map {print "\t\t$_ = $comp{$_}\n"} sort keys %comp;
		}

		map {$target_info{$indexes[$i]}{$_} = $comp{$_}} keys %comp;
	}
	map {my $i = $_; map {print STDERR "$i:$_ = $target_info{$i}{$_}\n"} sort keys %{$target_info{$i}}} sort keys %target_info;
}

$in_sam = openInput($samfile && $samfile =~ /\.bam$/ ? "samtools view -h $samfile |" : $samfile);
$out = openOutput($output);
$in_ref = openInput($seqfile);
makeIndex($in_ref,\%ref_index);
$readseq_Fun= creadsequence();

while(<$in_sam>){
	
	chomp $_;
	my $fistchar=substr($_,0,1);
	if(index($_,"\@PG")>-1){
			$program=getPgname($_);
			#print "$program\n";
			#if user input programme name, just use it.
		if($inputProgram ne ""){
			$program=$inputProgram;
		}
			
		next;
	}
		
	if (!($fistchar eq "\@")) {
		$startSeq=1;
		#print "$_  ".index($_,"\@PG")."\n";
	}

	if($startSeq>0){
		#print $_."\n";
		my @line=split(/\t/,$_);
		$readId=$line[0];
		my $refseq=$line[2];
		my ($chrome,$start,$end)=parseSeqname($line[2]);
		my $offset=$line[3];
		my $pos=$line[3];
		my $seq=$line[9];
		my $quality=$line[10];
		my $matchStart=$start+$offset;
		my $matchEnd=$matchStart + length($seq);
		
		my %tag=getTag(\@line,11);
		#print Dumper(\%tag);
		my @meths  = exists $tag{"YM"} ? split(/,/,$tag{"YM"}) : ();
		my @umeths = exists $tag{"YU"} ? split(/,/,$tag{"YU"}) : ();
		
# %results{chromosome}{position} = [ methyRead, unmethyRead, position, methyBase, adjBases, methyQual, unmethyQual ]
		for(my $i=0;$i<scalar(@meths);++$i){
			my $meth=$meths[$i];
			my $meth_Qual=getQual($quality,$meth);
			next if $meth_Qual < $minQual;
			my $meth_position=$start+$offset+$meth-2;
			
			if(exists $result{$chrome}{$meth_position}){
				 $result{$chrome}{$meth_position}[0] = $result{$chrome}{$meth_position}[0]+1;
				# print "hit\n"
				 $result{$chrome}{$meth_position}[5] = $result{$chrome}{$meth_position}[5]+$meth_Qual;
			}else {
				$result{$chrome}{$meth_position}=[(1, 0, $meth_position, "", "",$meth_Qual,0)];
			}

			if(exists $refseqMeth{$refseq}){
				my $position_hash=$refseqMeth{$refseq};
				if (!(exists $$position_hash{$meth_position})){
					$$position_hash{$meth_position}=0;
				}
			 }else {
			 	my %position=($meth_position=>0);
			 	$refseqMeth{$refseq}=\%position;
			 }
		}
		
		for(my $i=0;$i<scalar(@umeths);++$i){
			my $umeth=$umeths[$i];
			my $umeth_Qual=getQual($quality,$umeth);
			next if $umeth_Qual < $minQual;
			my $umeth_position=$start+$offset+$umeth-2;

			if(exists $result{$chrome}{$umeth_position}){
				 $result{$chrome}{$umeth_position}[1] = $result{$chrome}{$umeth_position}[1]+1;
				 $result{$chrome}{$umeth_position}[6] = $result{$chrome}{$umeth_position}[6]+$umeth_Qual;
				
			}else {
				$result{$chrome}{$umeth_position}=[(0, 1, $umeth_position, "", "", 0, $umeth_Qual)];
			}
			
			if(exists $refseqMeth{$refseq}){
				my $position_hash=$refseqMeth{$refseq};
				if (!(exists $$position_hash{$umeth_position})){
					$$position_hash{$umeth_position}=0;
				}
			 }else {
			 	my %position=($umeth_position=>0);
			 	$refseqMeth{$refseq}=\%position;
			 }
		}
		
		#if(scalar(keys %result) > $limit){
		#	 print scalar(%result)."\n";
		#	 getStand(\%refseqMeth,$seqfile);
		#	 outputTempFile(\%result,\%posIncluded,$output.".".$tempFileNumber);
		#	 $tempFileNumber++;
		#	 undef %result;
		#	 }
		
		$preMatchEnd=$matchEnd if ($preChrom && $preChrom ne $chrome) || !$preMatchEnd || $preMatchEnd < $matchEnd;
		$preChrom=$chrome;

		getStand();
		outputFile();
		undef %refseqMeth;
		undef %result;
	}	
}

$samfile || close($in_sam) || die "Can't close $samfile\n";

$seqfile || close($in_sam) || die "Can't close $samfile\n";

$output || close($out) || die "Can't close $samfile\n";

#print "\n";
#print scalar(%result);
#print "\n", scalar(keys %result), "\n";



sub outputFile{
#	foreach my $chr (sort keys %result){
#		foreach my $pos (sort {$a<=>$b} keys %{$result{$chr}}){
	foreach my $chr (keys %result){
		foreach my $pos (keys %{$result{$chr}}){
			if($incFile && !exists $posIncluded{$chr}{$pos}){ next; }
			if($excFile &&  exists $posExcluded{$chr}{$pos}){ next; }
	
			my	$updown= $result{$chr}{$pos}[4]; # bases in [$leftOff, $rightOff]
		
			if($target_bases ne "" ){
				#print "$pos, $updown\n";
#				foreach my $i (sort {$result{$chr}{$pos}[3] eq 'C' ? $a<=>$b : $b<=>$a} keys %target_info){ # $i = the position of C in $updown
				foreach my $i (keys %target_info){ # $i = the position of C in $updown
					foreach my $s (keys %{$target_info{$i}}){
						my $_updown = substr($updown, $i-$leftOff, length($s));
						next if $s ne $_updown;
						output($chr, $pos, $i, length($s), $target_info{$i}{$s}, $_updown);
					}
				}
			}else {
	 			output($chr, $pos, 0, length($updown),"","");
			}
		}
	}

	map {
		my $chr = $_;
		map {
			my $pos = $_;
			map {
				my $dir = $_;
				map {print $out join("\t", @$_), "\n"} @{$buf{$chr}{$pos}{$dir}}
			} sort keys %{$buf{$chr}{$pos}}
		} sort {$a <=> $b} keys %{$buf{$_}}
	} sort{$a cmp $b} keys %buf;
	undef %buf;
}

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

	return *STDOUT unless $fileName;

	my ($fd);
	open($fd, $fileName =~ /.gz$/ ? "| gzip -c > $fileName" : $fileName =~ /.bz(ip)?2$/ ? "| bzip2 -z -c > $fileName" : ">$fileName") || die("Open error: $fileName");
	return $fd;
}
