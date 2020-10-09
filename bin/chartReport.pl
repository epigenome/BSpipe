#!/usr/bin/perl -w
#!/bin/perl -w
#!/usr/local/bin/perl -w

use strict;
use warnings;
use SOAP::Lite;
use HTTP::Cookies;
use Getopt::Long qw(:config no_ignore_case);

my @avail_idTypes = qw (AFFYMETRIX_3PRIME_IVT_ID AFFYMETRIX_EXON_GENE_ID AFFYMETRIX_SNP_ID AGILENT_CHIP_ID AGILENT_ID AGILENT_OLIGO_ID ENSEMBL_GENE_ID ENSEMBL_TRANSCRIPT_ID ENTREZ_GENE_ID FLYBASE_GENE_ID FLYBASE_TRANSCRIPT_ID GENBANK_ACCESSION GENOMIC_GI_ACCESSION GENPEPT_ACCESSION ILLUMINA_ID IPI_ID MGI_ID OFFICIAL_GENE_SYMBOL PFAM_ID PIR_ID PROTEIN_GI_ACCESSION REFSEQ_GENOMIC REFSEQ_MRNA REFSEQ_PROTEIN REFSEQ_RNA RGD_ID SGD_ID TAIR_ID UCSC_GENE_ID UNIGENE UNIPROT_ACCESSION UNIPROT_ID UNIREF100_ID WORMBASE_GENE_ID WORMPEP_ID ZFIN_ID);

my @avail_categories = qw (BBID BIND BIOCARTA BLOCKS CGAP_EST_QUARTILE CGAP_SAGE_QUARTILE CHROMOSOME COG_NAME COG_ONTOLOGY CYTOBAND DIP EC_NUMBER ENSEMBL_GENE_ID ENTREZ_GENE_ID ENTREZ_GENE_SUMMARY GENETIC_ASSOCIATION_DB_DISEASE GENERIF_SUMMARY GNF_U133A_QUARTILE GENETIC_ASSOCIATION_DB_DISEASE_CLASS GOTERM_BP_2 GOTERM_BP_1 GOTERM_BP_4 GOTERM_BP_3 GOTERM_BP_FAT GOTERM_BP_5 GOTERM_CC_1 GOTERM_BP_ALL GOTERM_CC_3 GOTERM_CC_2 GOTERM_CC_5 GOTERM_CC_4 GOTERM_MF_1 GOTERM_MF_2 GOTERM_CC_FAT GOTERM_CC_ALL GOTERM_MF_5 GOTERM_MF_FAT GOTERM_MF_3 GOTERM_MF_4 HIV_INTERACTION_CATEGORY HIV_INTERACTION_PUBMED_ID GOTERM_MF_ALL HIV_INTERACTION KEGG_PATHWAY HOMOLOGOUS_GENE INTERPRO OFFICIAL_GENE_SYMBOL NCICB_CAPATHWAY_INTERACTION MINT PANTHER_MF_ALL PANTHER_FAMILY PANTHER_BP_ALL OMIM_DISEASE PFAM PANTHER_SUBFAMILY PANTHER_PATHWAY PIR_SUPERFAMILY PIR_SUMMARY PIR_SEQ_FEATURE PROSITE PUBMED_ID REACTOME_INTERACTION REACTOME_PATHWAY PIR_TISSUE_SPECIFICITY PRINTS PRODOM PROFILE SMART SP_COMMENT SP_COMMENT_TYPE SP_PIR_KEYWORDS SCOP_CLASS SCOP_FAMILY SCOP_FOLD SCOP_SUPERFAMILY UP_SEQ_FEATURE UNIGENE_EST_QUARTILE ZFIN_ANATOMY UP_TISSUE TIGRFAMS SSF UCSC_TFBS);


my ($helpFlag, $geneFile, $idTypeStr, $userListName, $categoryLists, $output, $account, $idFlag, $catFlag);

GetOptions(
	"h|?|help"      => \$helpFlag,
	"i|input=s"     => \$geneFile,	    ## gene file  
	"t|type=s"      => \$idTypeStr,	    ## id type
	"l|listname=s"  => \$userListName,	    ## list file
	"c|category=s"  => \$categoryLists,	    ## category idx
	"o|ouput=s"     => \$output,	    ## output file
	"a|account=s"   => \$account,	    ## output file
	"si|showid"	    => \$idFlag,
	"sc|showcat"    => \$catFlag,
) || die "\n";

checkOptions();


my $GENEFILE = openInput($geneFile);
my $startLine = 0;
my $fLen = length($geneFile);
my $fExt =substr($geneFile,$fLen-3,$fLen) ;
$startLine = 3 if ($fExt eq "gct");
# print "StartLine:$startLine\n";
my $linenum = 0;


my $tmpinputIds="";
while(<$GENEFILE>){
	chomp;
	my $line = $_;
	$linenum++;
	next if ($linenum <= $startLine);
	my @cols = split("\t",$line);
	$tmpinputIds.=$cols[0].",";
}
close($GENEFILE);


#$ENV{PERL_LWP_SSL_VERIFY_HOSTNAME}=0;
my $soap = SOAP::Lite			     
     -> uri('http://service.session.sample')		
     -> proxy('http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService',
#     -> proxy('http://129.43.1.164/webservice/services/DAVIDWebService',
		ssl_opts => [ SSL_verify_mode => 'SSL_VERIFY_NONE' ],
		cookie_jar => HTTP::Cookies->new(ignore_discard=>1));

 #user authentication by email address
 #For new user registration, go to http://david.abcc.ncifcrf.gov/webservice/register.htm
my $check = $soap->authenticate($account)->result;
#  	print "\nUser authentication: $check\n";

if (lc($check) eq "true") { 
 
 #list conversion types
 my $conversionTypes = $soap ->getConversionTypes()->result;
 #print  "\nConversion Types: \n$conversionTypes\n"; 
	 
 #list all annotation category names
 my $allCategoryNames= $soap ->getAllAnnotationCategoryNames()->result;	 	  	
 #print  "\nAll available annotation category names: \n$allCategoryNames\n"; 
 
	my $inputIds = substr($tmpinputIds,0,length($tmpinputIds)-1);
	# print "\n\nTmpIds:".$tmpinputIds."\n";
	#print "\n\nInputIds:".$inputIds."\n";

	# idType
	# my $idType = 'AFFYMETRIX_3PRIME_IVT_ID';
	#my $idType = $avail_idTypes[$idTypeStr];
	my $idType = $idTypeStr;

	print "IdType:".$idType."\n";


# list name 
#  my $listName = 'make_up';
my $listName =  $userListName;
my $listType = 0;

#to add background list, set listType=1
my $list = $soap ->addList($inputIds, $idType, $listName, $listType)->result;
#print "\n$listName of list was mapped\n"; 
  	
#list all species  names
my $allSpecies= $soap ->getSpecies()->result;	 	  	
# print  "\nAll species: \n$allSpecies\n"; 
#list current species  names
my $currentSpecies= $soap ->getCurrentSpecies()->result;	 	  	
#print  "\nCurrent species: \n$currentSpecies\n"; 

#set user defined species 
#my $species = $soap ->setCurrentSpecies("1")->result;

#print "\nCurrent species: \n$species\n"; 
 

# category
#set user defined categories 
# my $categories = $soap ->setCategories("BBID,BIOCARTA,COG_ONTOLOGY,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,UP_SEQ_FEATURE")->result;
my $categories = $soap ->setCategories($categoryLists)->result;


#to user DAVID default categories, send empty string to setCategories():
# my $categories = $soap ->setCategories("")->result;



print "Valid categories: \n$categories\n\n";  

 
# open (chartReport, ">", "chartReport.txt");
open (chartReport, ">", $output);
print chartReport "Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n";
#close chartReport;

#open (chartReport, ">>", "chartReport.txt");
#getChartReport 	
my $thd=0.1;
my $ct = 2;
my $chartReport = $soap->getChartReport($thd,$ct);
	my @chartRecords = $chartReport->paramsout;
	#shift(@chartRecords,($chartReport->result));
	#print $chartReport->result."\n";
  	print "Total chart records: ".(@chartRecords+1)."\n";
# 	print "\n ";
	#my $retval = %{$chartReport->result};
	
	# added by dsryu
	if (defined($chartReport->result)){
	

	my @chartRecordKeys = keys %{$chartReport->result};
	
	#print "@chartRecordKeys\n";
	
	my @chartRecordValues = values %{$chartReport->result};
	
	my %chartRecord = %{$chartReport->result};
	my $categoryName = $chartRecord{"categoryName"};
	my $termName = $chartRecord{"termName"};
	my $listHits = $chartRecord{"listHits"};
	my $percent = $chartRecord{"percent"};
	my $ease = $chartRecord{"ease"};
	my $Genes = $chartRecord{"geneIds"};
	my $listTotals = $chartRecord{"listTotals"};
	my $popHits = $chartRecord{"popHits"};
	my $popTotals = $chartRecord{"popTotals"};
	my $foldEnrichment = $chartRecord{"foldEnrichment"};
	my $bonferroni = $chartRecord{"bonferroni"};
	my $benjamini = $chartRecord{"benjamini"};
	my $FDR = $chartRecord{"afdr"};
	
		print chartReport "$categoryName\t$termName\t$listHits\t$percent\t$ease\t$Genes\t$listTotals\t$popHits\t$popTotals\t$foldEnrichment\t$bonferroni\t$benjamini\t$FDR\n";
	
	
		for my $j (0 .. (@chartRecords-1))
		{			
			%chartRecord = %{$chartRecords[$j]};
			$categoryName = $chartRecord{"categoryName"};
			$termName = $chartRecord{"termName"};
			$listHits = $chartRecord{"listHits"};
			$percent = $chartRecord{"percent"};
			$ease = $chartRecord{"ease"};
			$Genes = $chartRecord{"geneIds"};
			$listTotals = $chartRecord{"listTotals"};
			$popHits = $chartRecord{"popHits"};
			$popTotals = $chartRecord{"popTotals"};
			$foldEnrichment = $chartRecord{"foldEnrichment"};
			$bonferroni = $chartRecord{"bonferroni"};
			$benjamini = $chartRecord{"benjamini"};
			$FDR = $chartRecord{"afdr"};			
			print chartReport "$categoryName\t$termName\t$listHits\t$percent\t$ease\t$Genes\t$listTotals\t$popHits\t$popTotals\t$foldEnrichment\t$bonferroni\t$benjamini\t$FDR\n";				 
		}		  	
	
	}
	close chartReport;
#	print "\n$output generated\n";
} 
# added by dsryu
#----------------------------------------------
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
	$geneFile  = shift(@ARGV) if !defined $geneFile  && @ARGV > 0;
	$idTypeStr  = shift(@ARGV) if !defined $idTypeStr  && @ARGV > 0;
	$userListName  = shift(@ARGV) if !defined $userListName  && @ARGV > 0;
	if (!defined $categoryLists|| $categoryLists eq ""){
		$categoryLists = "";
	}
	$output  = shift(@ARGV) if !defined $output  && @ARGV > 0;

	if ($idFlag ) { print join("\n", "Avaliable ID types: ", sort(@avail_idTypes), ""); exit; }
	if ($catFlag) { print join("\n", "Avaliable categories: ", sort(@avail_categories), ""); exit; }

	print "geneFile:".$geneFile."\n";
	print "idTypeStr:".$idTypeStr."\n";
	print "listName:".$userListName."\n";
	
	print "categoryLists:".$categoryLists."\n";
	print "output:".$output."\n";
	# if ($helpFlag || !$geneFile ||!$idTypeStr ||!$userListName || !$output){
	if ($helpFlag || !$account || !$geneFile ||!defined $idTypeStr ||!$userListName || !defined $categoryLists || !$output)
	{
		die("Arguments: [-show] -a account [[-i] geneFile] [[-t] idtype(string)] [[-l] list title] [[-c] categorylist('category1,category2...')] [[-o] out_file]\n"
		  . "\tFor new user registration, go to http://david.abcc.ncifcrf.gov/webservice/register.htm\n"
		  );
	}
}

__END__
		
