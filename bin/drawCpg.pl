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
use GD::Image;

my ($helpFlag, $refFile, $tarFile, $outPrefix, $htmlFile, $verbose, $quiet, $nocrdFlag, $field, $strand, $posFlag, $metFlag, $covFlag, $numFlag, $imgWidth, $percent, $idFirst, $naFlag, $bedFlag);

GetOptions(
	"h|?|help"		=> \$helpFlag,	
	"target=s"		=> \$tarFile,		## input file
	"reference=s"	=> \$refFile,		## input file
	"s|sample=s"	=> \$samFile,		## input file
	"output=s"		=> \$outPrefix,	## output file
	"html=s"			=> \$htmlFile,
	"verbose+"		=> \$verbose,		## verbose output
	"quiet"			=> \$quiet,
	"nocoord!"		=> \$nocrdFlag,
	"field=i"		=> \$field,
	"strand=s"		=> \$strand,
	"pos!"			=> \$posFlag,
	"num!"			=> \$numFlag,
	"na!"				=> \$naFlag,
	"width=i"		=> \$imgWidth,
	"percent!"		=> \$percent,
	"id"				=> \$idFirst,
	"met!"			=> \$metFlag,
	"cov!"			=> \$covFlag,
	"bed!"			=> \$bedFlag,
) || die "\n";

my ($Fw, $Rv, $Both) = (0, 1, 3);

checkOptions();

my $binSize = 10000;
my ($topBotMargin, $sideMargin, $chHeight, $chWidth) = (20, 20, 14, 6);
my $cellWidth = 10;
my $cellHeight = 10;
my %colors;
my $maxColor = 150;
my $legendSpace = $metFlag ? $covFlag ? 60 : 40 : $covFlag ? 40 : 0;
my $longestName = 0;
my (%bins, %tars, %refs, @sams, @order, $image, $title);

loadTarget($tarFile);
loadReference($refFile) if $refFile;
loadSample($samFile);
output($outPrefix);


#-------------------------------------------------------------------------------

sub loadTarget
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);
	my $no = 0;

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my @a = split;
		my ($id, @loc) = @a[$idFirst ? (0..$#a) : ($#a,0..$#a-1)];
		@loc = split(/[:-]/, $loc[0]) if (@loc == 1);
		$loc[1] =~ s/,//g;
		$loc[2] =~ s/,//g;
		push(@order, $loc[0]) if !exists $tars{$loc[0]};
		$tars{$loc[0]}{$loc[1]} = [@loc[1,2], $id];
		map {push(@{$bins{$loc[0]}{$_}}, $tars{$loc[0]}{$loc[1]})} int($loc[1]/$binSize)..int($loc[2]/$binSize);
		$no++;
	}

	close($in) if defined $fileName;
	die "There should be one target with -nocoord\n" if $nocrdFlag && $no != 1;
}

sub loadReference
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);

	while (<$in>)
	{
		next if /^#/ || /^\s*$/;
		my @a = split(/[\t\r\n]/);
		my $b = int($a[1]/$binSize);
		next if !exists $bins{$a[0]}{$b};

		foreach my $v (@{$bins{$a[0]}{$b}})
		{
			if ($v->[0] <= $a[1] && $a[1] <= $v->[1])
			{
				$refs{$a[0]}{$a[1]} = 1;
			}
		}
	}

	close($in) if defined $fileName;
}

sub loadSample
{
	my ($fileName) = @_;

	print "Loading ", $fileName || 'STDIN', " ...\n" if !$quiet;
	my $in = openInput($fileName);

	while (<$in>)
	{
		next if /^\s*$/;
		my @a = split(/[\t\r\n]/);

		if (/^#/)
		{
			@sams = @a[$field..$#a];
			map {$longestName = length($_) if(length($_) > $longestName)} @a[$field..$#a];
		}
		elsif ($nocrdFlag) # assume one target
		{
			foreach my $u (values %tars)
			{
				foreach my $v (values %$u)
				{
					push(@$v, [$a[1], map {$percent ? $_/100 : $_} @a[$field..$#a]]);
				}
			}
		}
		elsif ($naFlag || !/\bNA\b/)
		{
			my $b = int($a[1]/$binSize);
			next if !exists $bins{$a[0]}{$b} || ($refFile && !exists $refs{$a[0]}{$a[1]});

			foreach my $v (@{$bins{$a[0]}{$b}})
			{
				if ($v->[0] <= $a[1] && $a[1] <= $v->[1])
				{
					push(@$v, [$a[1], map {$percent ? $_/100 : $_} @a[$field..$#a]]);
				}
			}
		}
	}

	die "Zero samples\n" if !@sams;
	close($in) if defined $fileName;
}

sub output
{
	my ($fileName) = @_;

	foreach my $chr (@order)
	{
		foreach my $beg (sort {$a<=>$b} keys %{$tars{$chr}})
		{
			my $rgene = $tars{$chr}{$beg};
			{
				next if @$rgene <= 3;
				$title = $rgene->[2] || "${chr}_$rgene->[0]_$rgene->[1]";
				print "Drawing $title ...\n" if !$quiet;

				if ($bedFlag)
				{
					my $out = openOutput("$outPrefix$title.bed");
					print $out join("\t", '#chr', 'position', @sams), "\n";
					foreach my $rsam (@$rgene[3..$#$rgene])
					{
						print $out join("\t", $chr, @$rsam), "\n";
					}
					close($out);
				}

				draw($chr, $rgene);
			}
		}
	}
}

sub draw
{
	my ($chr, $rgene) = @_; # (start, end, ID, CpG, CpG, ...)
	my $imgFile = "$outPrefix$title.png";
	my $cpgNum = @$rgene - 3;
	my ($width, $height);

	if ($imgWidth)
	{
		$width = $imgWidth;
		$cellWidth = ($width - 2*$sideMargin - ($longestName * $chWidth)) / $cpgNum;
	}
	else
	{
		$width = 2*$sideMargin + ($longestName * $chWidth) + ($cpgNum*$cellWidth);
		$width = 350 if($legendSpace && $width < 350);
	}

	my $posHeight = 0;
	if ($posFlag)
	{
		my $len = 0;
		map {my $n = length($rgene->[$_][0]-$rgene->[0]+1); $len = $n if $len < $n} 3..$#$rgene;
		$posHeight += $chWidth * $len;
	}

	$image = new GD::Image($width, 2*$topBotMargin + (@sams*($cellHeight*($strand==$Both?2:1)+3))+ $chHeight + $legendSpace + $posHeight);
	
	#print " 2*$sideMargin + ($longestName * $chWidth) + ($cpgNum*$cellWidth), 
	#	$topBotMargin + ($lineNum*$cellHeight)+ $chHeight + $topBotMargin + 70 \n";
	my ($stStart, $stEnd) = (0, 1);
	$stStart = 1 if($strand == $Rv);
	$stEnd = 0 if($strand == $Fw);

	$image->interlaced('true');
	
	initImage();

	my ($rowOffset, $columOffset, $str, $drawTop);
	my ($top, $left) = ($topBotMargin+$legendSpace, $sideMargin); 

	#print "$i, $arrGeneClust{geneName} - \n";
	
	my $label = ($rgene->[2] || '') . " $chr:$rgene->[0]-$rgene->[1]";
	$label .= " ($cpgNum CpGs)" if $numFlag;
	$image->string(GD::gdMediumBoldFont, $left, $top, $label, $colors{black});

	$top += $chHeight+3;
	for($j = 0; $j < @sams; $j++)
	{
		$left = $sideMargin;
		$image->string(GD::gdSmallFont, $left, $top, $sams[$j] ,$colors{black});
		if(defined $htmlFile && $j< @sams){
			printf($out " 
				<AREA SHAPE=RECT COORDS='%d,%d,%d,%d'
					HREF=\"%s.html\"
					ALT='%s'
					title='%s'>\n", 
						$left, $top, $left+length($sams[$j])*($chWidth), $top+$chHeight,
						lc($sams[$j]),$sams[$j],$sams[$j]);
		}

		$left = $sideMargin + ($longestName * $chWidth) + 3;
		for($cid = 0; $cid < $cpgNum; $cid++)
		{
			###### set coverage for undetected genes in cell lines.
#			if($#{$arrGeneClust[$cid]->{$sams[$j]}} == -1){ 
#				for($k = 0; $k < $arrGeneClust[$cid]->{cpgNum}; $k++){
#					$arrGeneClust[$cid]->{$sams[$j]}->[$k]->[$Fw]->{cover} = 0;
#					$arrGeneClust[$cid]->{$sams[$j]}->[$k]->[$Rv]->{cover} = 0;
#				}
#			}
		
			for(my $st = $stStart; $st <= $stEnd; $st++)
			{
				$left = $sideMargin + ($longestName * $chWidth) + 3;
				for($k = 0; $k < $cpgNum; $k++)
				{
					if($strand == $Both && $st == $Rv) {	$drawTop = $top + $cellHeight;	}
					else{	$drawTop = $top;	}

					$image->rectangle($left, $drawTop, $left+$cellWidth, $drawTop+$cellHeight,$colors{black});
				
					#my $rate = $refs{$chr}{$rgene->[$k+3]}[$j];
					my $rate = $rgene->[$k+3][$j+1]; # (start, end, id, ((pos, value, ...) ..))
					$rate = undef if $rate && $rate eq 'NA';
					my $cover = 1000;
#					$cover = $arrGeneClust[$cid]->{$sams[$j]}->[$k]->[$st]->{cover};
					$cover = $maxColor if($cover > $maxColor);
					$cover = $maxColor-$cover;

					if(defined $rate)
					{
						$image->filledRectangle($left+1, $drawTop+1, $left+$cellWidth-1, $drawTop+$cellHeight-1, $colors{cover}[$cover/10]);
#						$str = $sams[$j] . ":" . $arrGeneClust[$cid]->{cpgPos}->[$k]->[$st]->{cpgid};			

						if($j < @sams  && $sams[$j] ne "all")
						{
#							$fileNames[$j] = basename($fileNames[$j]);
#							$fileNames[$j] =~ s/\.bsmapper\.gff$/\.bsmapper\.best\.msa/g;
#							printf($out " 
#							<AREA SHAPE=RECT COORDS='%d,%d,%d,%d'
#								HREF=\"https://projects-dev.cgb.indiana.edu/jeochoi/cgi-bin/methylation/2009_10_09/bisulfite.pl?id=%s&file=%s\"
#								ALT='%s'
#								title='%s'>\n", 
#									$left+1, $drawTop+1, $left+$cellWidth-1, $drawTop+$cellHeight-1,
#									$str, $fileNames[$j],$str.":".$sams[$j],$str.":".$sams[$j]) if(defined $htmlFile);
						}

						$height = $cellHeight - $rate*$cellHeight+1;
						$height = $cellHeight-1 if($height > $cellHeight-1);
						$image->filledRectangle($left+1, $drawTop+$height, $left+$cellWidth-1, $drawTop+$cellHeight-1, $colors{cpg}[$cover/10]);				
					}

					$left += $cellWidth;
				}
			}
		}
#			print "\n";
		$top += $cellHeight*($strand==$Both?2:1)+3;
	}

	if ($posFlag)
	{
		$top += $posHeight;
#		for($cid = 0; $cid < $cpgNum; $cid++)
		{
#			for(my $st = $stStart; $st <= $stEnd; $st++)
			{
				$left = $sideMargin + ($longestName * $chWidth) + 3;
				for($k = 0; $k < $cpgNum; $k++)
				{
					#$image->stringUp(GD::gdSmallFont, $left, $top, $rgene->[$k+3]-$rgene->[0]+1, $colors{black});
					$image->stringUp(GD::gdSmallFont, $left, $top, $rgene->[$k+3][0]-$rgene->[0]+1, $colors{black});
					$left += $cellWidth;
				}
			}
		}
	}

	open(ARAFH, ">$imgFile") || die print "Could not create an image file in temp directory : $!\n";
	binmode ARAFH;
	print ARAFH $image->png;
	close ARAFH;

	if(defined $htmlFile){	
		print $out "</map>\n<img usemap='#genemap' border=1 src='$imgFile'>\n</body></html>\n";
		close($out);
	}
}


sub initImage{
	$colors{white} = $image->colorAllocate(255,255,255);
	$colors{black} = $image->colorAllocate(0,0,0);
	$colors{navy} = $image->colorAllocate(0,0,150);	
	$colors{gray} = $image->colorAllocate(100,100,100);

	my ($top, $left, $i, $k);
	for($i = 0; $i <= $maxColor/10; $i++){
		$colors{cpg}[$i] = $image->colorAllocate($i*12,$i*12,255);
	}
	for($i = 0; $i <= $maxColor/10; $i++){
		$colors{cover}[$i] = $image->colorAllocate(255,190+($i*4),$i*11);
	}
	
	$top = $topBotMargin;
=head
	$top = $topBotMargin - 10;
	
	$image->filledRectangle($sideMargin, $top, $sideMargin+ 180, $top+16, $colors{navy});
	$image->rectangle($sideMargin, $top, $sideMargin+ 180, $top+17, $colors{black});
	$image->string(GD::gdMediumBoldFont, $sideMargin+3, $top+2 , "Interactive visualization" ,$colors{white});
	
	$top += 20;	
	$image->string(GD::gdSmallFont, $sideMargin, $top, "Click cell line names or tables to see detail information"  ,$colors{navy});
	
	$top += 20;
=cut
	
	$k = 150; ### for the left space for the string 'Coverage' ..
	$left = $sideMargin+$k+1;
	#$top = $topBotMargin;
	##########################3###### color legend.
#	$image->string(GD::gdMediumBoldFont, $sideMargin+$k, $top , "low" ,$colors{black});
#	$image->string(GD::gdMediumBoldFont, $sideMargin+$k+$maxColor-30, $top, "high" ,$colors{black});
#	$top += $chHeight;

	$j = $maxColor/10;

	if ($metFlag)
	{
		$image->string(GD::gdMediumBoldFont, $sideMargin, $top + 1, "CpG Methylation" ,$colors{black});

		#($top, $left) = ($topBotMargin+$chHeight+21, $sideMargin+$k+1);
		for($i = 0; $i <= $j; $i++){
			$image->filledRectangle($left, $top+1, $left+10, $top+16, $colors{cover}[0]);
			$image->filledRectangle($left, $top+($j-$i)+1, $left+10, $top+16,$colors{cpg}[0]);
			$left += 10;
		}
		#$top = $topBotMargin+$chHeight+21;
		$image->rectangle($sideMargin+$k, $top, $sideMargin+($maxColor+10)+2+$k, $top+17,$colors{black});
		$top += 21;
	}

	if ($covFlag)
	{
		$image->string(GD::gdMediumBoldFont, $sideMargin, $top + 1, "Read coverage" ,$colors{black});
		$left = $sideMargin+$k+1;
		#($top, $left) = ($topBotMargin+$chHeight, $sideMargin+$k+1);
		for($i = 0; $i <= $j; $i++){
			$image->filledRectangle($left, $top+1, $left+10, $top+8,$colors{cover}[$j-$i]);
			$image->filledRectangle($left, $top+9, $left+10, $top+16,$colors{cpg}[$j-$i]);
			$left += 10;
		}
		$image->rectangle($sideMargin+$k, $top, $sideMargin+($maxColor+10)+2+$k, $top+17,$colors{black});
	}
	################################ until here, color legend.
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
	if ($helpFlag || !$tarFile || !$samFile || (!$field && !$nocrdFlag))
	{
		die("Arguments: [-na] [-num] [-met] [-cov] [-percent] [-w width] [-pos] [-strand 0|1|3] -t target_file [-r ref_file(bed)] -s sample_file <-f field | -nocoord> [-bed] -o out_file [-v] [-q]\n"
		  );
	}

	if ($nocrdFlag)
	{
		$field = 0;
		undef $strand;
		undef $posFlag;
	}
	else
	{
		$field--;
	}

	if (!defined $strand)
	{
		$strand = $Fw;
	}
	elsif ($strand != $Fw && $strand != $Rv && $strand != $Both)
	{
		die "Please enter either $Fw, $Rv, or $Both for -strand\n";
	}

	if (! -d $outPrefix) { $outPrefix = "$outPrefix." if $outPrefix !~ /\.$/; }
	else                 { $outPrefix = "$outPrefix/" if $outPrefix !~ /\/$/; }
}
