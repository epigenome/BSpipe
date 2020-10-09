#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

## Author: Jeong-Hyeon Choi
## Author: Dong-Sung Ryu
## Date: 2012

## Description

## Usage:
## CHECK usage with option  -h



use Getopt::Long qw(:config no_ignore_case);

my ($helpFlag, @inFiles, $inFai, $outFile, $verbose, $quiet);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"input=s"		=> \@inFiles,		## input file
	"fai=s"		=> \$inFai,		## input file
	"output=s"		=> \$outFile,		## output file
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
) || die "\n";

checkOptions();

my (%chrs, %refs, %sites);

loadChr($inFai);

my @fins;
my @bufs;

# check position sorting added by dsryu at Sep. 15th
my @befposbuf;

my $posUnit = 1000000; # 1M

for (my $i=0;$i<@inFiles;$i++){

	print STDERR "Opening $inFiles[$i] files.\n" if (!$quiet);
	$fins[$i] = openInput($inFiles[$i]);
	$bufs[$i] = "";

	# check position sorting added by dsryu at Sep. 15th
	$befposbuf[$i] = -1;
}

my $fout = openOutput($outFile);
foreach $key (sort {$chrs{$a} <=> $chrs{$b} } keys %chrs) {
	#print STDERR "\tProcessing ".$key.":".$chrs{$key}."...\n" if (!$quiet) ;
	my $chrorder = $chrs{$key};

	
	#file pointer
	my $posRange = 0;
	my $condition = 1;

	while ($condition ){
		
		# for all files ending check
		my $feofcond = 0;
			
		#For chromosome checking in each files
		my $chrcond = 0;


		for (my $i=0;$i<@fins;$i++){
	#		print STDERR "\t\tLoad line  ".$i." -th buf or file \n" if (!$quiet);
			my $line = "";
			my $fin = $fins[$i];
			if (eof($fin)){
	#			print STDERR "\t\tIncrease feofcond from $feofcond to ".($feofcond+1)." by ".$i."-th file eof\n" if (!$quiet);
				$feofcond++;
				
				# check end condition  sorting added by dsryu at Oct. 8th
				#------------
	#			print STDERR "\t\tIncrease feofcond from $feofcond to ".($feofcond+1)." by ".$i."-th file eof\n" if (!$quiet);
				$chrcond++;
				#------------
			}
			if ($bufs[$i] eq ""){
				if (!eof($fin)){	
					#$line = <$fin>;
					#chomp($line);
					while (!eof($fin)){
						$line = <$fin>;
						chomp($line);
						next if ($line =~ /^#/ || $line =~ /^\s*$/);
						my @a = split("\t",$line);
						if ($a[1] >= 0){
							$bufs[$i] = $line;
							last;
						}			
						else{
							print STDERR $line." has minus position.\n";
						}
					}
					# check position
					#$bufs[$i] = $line;
					
				}
				# commented by dsryu at Oct. 8th
				#------------
				#else{ # all processed in this file
				#	print STDERR "\t\tIncrease chrcond from $chrcond to ".($chrcond+1)." by ".$i."-th by file eof\n" if (!$quiet);
					#$chrcond++;
				#}
				#------------
			}
			else{
				$line = $bufs[$i];
			}
	#		print STDERR "Printing line:".$line." from $i-th\n" if (!$quiet);
			if ($line ne ""){ # this is already processed at if eof($fin) else $chrcond++
				my @a = split("\t",$line);
				if ($a[0] ne $key && !eof($fin)){
	#				print STDERR "\t\tIncrease chrcond from $chrcond to ".($chrcond+1)." by ".$i."-th by line key:$key"."-a[0]".$a[0]."\n" if (!$quiet);
					$chrcond++;
				}
			}
			

		}
		#--------------
		#print STDERR "\tChecking chrcond: $chrcond and fins".scalar @fins." , feofcond:$feofcond \n" if (!$quiet);
		#print STDERR "Cur chr(key):".$key."\n" if (!$quiet);
		#for (my $i=0;$i< @bufs;$i++){
			#print STDERR "\t\tbufs[$i]=".$bufs[$i]."\n" if (!$quiet);
		#}
		if ($chrcond == @fins|| $feofcond == @fins ){
			# check position sorting added by dsryu at Sep. 15th
			for (my $i=0;$i<@fins;$i++){
				$befposbuf[$i] = 0;
			}
			$condition = 0;
			#commentd by dsryu at Oct. 8th
			#---------------------
			#last;
			#---------------------
		}
		
		for (my $i=0;$i<@fins;$i++){

			my $fin = $fins[$i];
			my $line = "";
		
			#while (!eof($fin) ){
			# edited by dsryu to fix file eof bug at Sep. 16th 
			while (!eof($fin) || $bufs[$i] ne "" ){
				# get Line from buf or file
				if ($bufs[$i] ne ""){
					$line = $bufs[$i];
	#				print STDERR "\t\tLoad line from $i-th Buf:".$line."\n" if (!$quiet);
					$bufs[$i]="";
				}
				else{
					defined ($line = <$fin>) or last;
					chomp($line);
				}
				my @a = split("\t",$line);

				if ($a[0] ne $key){
					$bufs[$i] = $line;
					last;
				}
				# check position sorting added by dsryu at Sep. 15th
				if ($befposbuf[$i] >  $a[1]){
					print STDERR "Error:Befposbuf[i]:".$befposbuf[$i]."\t".$a[1]."\n";
					die ("Error: This ".$inFiles[$i]." bed file is not sorted in position order:".$line."\n")  ;
				}
			
				# edited by dsryu at Oct. 8th to process remained buf
				#--------------------------
				# check position
				#if ( $a[1] < ($posRange + $posUnit) && $a[1] >= $posRange && $a[0] eq $key){
				#edited by dsryu at Oct. 8th
				if ( $a[1] < ($posRange + $posUnit) && $a[1] >= $posRange && $a[0] eq $key || ($condition == 0)){
				#---------------------
					#then process
	#				print STDERR "\t\tProcess:$line\n" if (!$quiet);
					$sites{$a[0]}{$a[1]}[0] = $a[2];
					$sites{$a[0]}{$a[1]}[1] += $a[3]; # total
					$sites{$a[0]}{$a[1]}[2] += $a[6]; # methylated
					$sites{$a[0]}{$a[1]}[3]  = $sites{$a[0]}{$a[1]}[3] && $sites{$a[0]}{$a[1]}[3] ne $a[5] ? '.' : $a[5];
					$bufs[$i]="";
					
					# check position sorting added by dsryu at Sep. 15th
					$befposbuf[$i] = $a[1];
				
				}
				else{
					#print STDERR "\t\tNext range......at $posRange in $i th file\n" if (!$quiet);
					# this line should be processed at next time.
					$bufs[$i] = $line;


					#print STDERR "\t\tCur buf[$i] =".$bufs[$i]."\n" if (!$quiet);
					last;
				}

			
			}
		}
		$posRange += $posUnit;
		#file writing
		# print STDERR "\t\tWriting at curPos :".$posRange. " of ".$key."\n" if (!$quiet);
		foreach my $pos (sort {$a <=> $b} keys %{$sites{$key}})
		{
			print $fout join("\t", $key, $pos, $sites{$key}{$pos}[0], $sites{$key}{$pos}[1], sprintf("%.3f", $sites{$key}{$pos}[2]/$sites{$key}{$pos}[1]), $sites{$key}{$pos}[3], $sites{$key}{$pos}[2]), "\n";
		}
		$sites{$key}= ();
		undef $sites{$key};

		#added  by dsryu at Oct. 8th
		#--------------------------
		last if ($condition == 0);
		#if ($posRange > 547199719){
			#die "There is some problem in chr. range\n";
		#}
		#--------------------------
	}

}

for (my $i=0;$i<@inFiles;$i++){
	close($fins[$i]);
}
close($fout);



sub loadChr{
	my $fainame = shift();
	my $ffai =openInput($fainame); 
	my $order = 0;
	while (<$ffai>){
		chomp();
		if (/^(\S+)/){
			$chrs{$1} = $order;
			$order++;
		}
	}
	close($ffai);
}


#-------------------------------------------------------------------------------
sub trim {
        my @result = @_;

        foreach (@result) {
                s/^\s+//;
                s/\s+$//;
        }

        return wantarray ? @result : $result[0];
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

	return *STDOUT unless defined $fileName;

	my ($fd);
	open($fd, $fileName =~ /.gz$/ ? "| gzip -c > $fileName" : $fileName =~ /.bz(ip)?2$/ ? "| bzip2 -z -c > $fileName" : ">$fileName") || die("Open error: $fileName");
	return $fd;
}

sub checkOptions
{
	push(@inFiles, @ARGV) if !@inFiles  && @ARGV > 0;
 	$inFai  = trim(shift(@ARGV)) if !defined $inFai  && @ARGV > 0;

	#print $inFai."\n";
	if ($helpFlag || !@inFiles ||!$inFai )
	{
		die("Arguments: [-i] in_file [-f fai file] [-o out_file] [-q]\n"
		  . "\t\n"
		  );
	}
}
