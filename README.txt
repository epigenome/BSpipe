BSpipe is an end-to-end anlysis pipeline for BS-seq


********************************************************************************
INSTALL
********************************************************************************
1. Uncompress a distribution file to an install directory.
shell> tar zxf bspipe-N.N.tgz

2. If necessary, set the path of bspipe to PATH.


********************************************************************************
PREREQUISITE
********************************************************************************
1. Install a mapping program. Tested version are
Bowtie   0.12.7 and 1.0.0
Bowtie2  2.1.0
BWA      0.5.9 and 0.7.8
SOAP2    2.21

2. Create a mapping configuration file or modify a file in conf/
--------------------------------------------------------------------------------
conf/bowtie2.conf
--------------------------------------------------------------------------------
#bowtie2 configuration

PROGRAM_DIR=/cluster/sw/encap/bowtie2-2.1.0/bin
PROGRAM=bowtie2
WATSON='QUALTYPE --norc --sensitive -k 20 -t -x REF SEQ'
CRICK='QUALTYPE --nofw --sensitive -k 20 -t -x REF SEQ'
INDEX='bowtie2-build -o 3 SEQ SEQ'
PHRED64=--phred64-quals
PHRED33=
--------------------------------------------------------------------------------

3. Install Bioconductor and the following libraries.
aod
cluster
corrgram
fpc
geneplotter
getopt
gplots
plotrix
plyr
reshape

4. For DMR annotation, modify or create a configuration file in conf/
and download necessary files in the first column, e.g., refGene.txt.gz,
from the UCSC genome browser or http://sourceforge.net/projects/bspipe/files/.
--------------------------------------------------------------------------------
conf/hg19ucsc.conf
--------------------------------------------------------------------------------
refGene.txt.gz	refseq	T	2	13	T	T	3	5	6	4	2000	2000	.	10	11	7	8	.
#ensGene.txt.gz	ensembl	T	2	13	T	F	3	5	6	4	2000	2000	.	10	11	7	8	.
#wgEncodeGencodeCompV17.txt.gz	gencode	T	2	13	T	T	3	5	6	4	2000	2000	.	10	11	7	8	.
wgRna.txt.gz	rna	T	5	.	T	F	2	3	4	7	0	0	.	.	.	.	.	.
#wgRna.txt.gz	rnaclass	F	10	.	F	F	2	3	4	7	0	0	.	.	.	.	.	-category
#lincRNAsTranscripts.txt.gz	linc	T	2	.	T	F	3	5	6	4	2000	2000	.	10	11	7	8	.
cpgIslandExt.txt.gz	cgi	F	5	.	F	F	2	3	4	.	2000	2000	11	.	.	.	.	.
#rmsk.txt.gz	repeat	F	11	.	F	F	6	7	8	10	0	0	.	.	.	.	.	.
rmsk.txt.gz	repcls	F	12	.	F	F	6	7	8	10	0	0	.	.	.	.	.	-category
#rmsk.txt.gz	repfam	F	13	.	F	F	6	7	8	10	0	0	.	.	.	.	.	.
#snp138.txt.gz	snp138	F	5	.	F	F	2	3	4	7	0	0	16	.	.	.	.	-category -noclass
--------------------------------------------------------------------------------

5. If a grid engine will be used, then create a grid engine configuration file or modify conf/sge.conf.
--------------------------------------------------------------------------------
conf/sge.conf
--------------------------------------------------------------------------------
...
hg19	index	'qsub -l h_vmem=20g -cwd -V -j y -o bsindex$$.log -b y -sync y'
hg19	bowtie	'qsub -l h_vmem=5g -cwd -V -j y -b y -sync y -o /dev/null'
hg19	bowtie2	'qsub -l h_vmem=5g -cwd -V -j y -b y -sync y -o /dev/null'
hg19	bwaw	'qsub -l h_vmem=9g -cwd -V -j y -b y -sync y -o /dev/null'
hg19	bwa	'qsub -l h_vmem=9g -cwd -V -j y -b y -sync y -o /dev/null'
hg19	soap2w	'qsub -l h_vmem=7g -cwd -V -j y -b y -sync y -o /dev/null'
common	common	'qsub -l h_vmem=2g -cwd -V -j y -b y -sync y -o /dev/null'
common	methylation	'qsub -l h_vmem=1g -cwd -V -j y -b y -sync y'
common	correlate	'qsub -l h_vmem=15g -cwd -V -j y -b y -sync y'
common	group	'qsub -l h_vmem=15g -cwd -V -j y -b y -sync y'
common	sample	'qsub -l h_vmem=100m -cwd -V -j y -b y -sync y'
...
--------------------------------------------------------------------------------


********************************************************************************
OPTIONAL PROGRAMS
********************************************************************************
BSmooth: http://rafalab.jhsph.edu/bsmooth/
QDMR: http://bioinfo.hrbmu.edu.cn/qdmr/
CpG_MPs: http://202.97.205.78/CpG_MPs/



********************************************************************************
TEST
********************************************************************************
1. Download a demo data (bspipe-demo.tgz) on http://sourceforge.net/projects/bspipe/files/

2. Uncompress and run the batch script
shell> tar zxf bspipe-demo.tgz
shell> cd demo

3. Read README



********************************************************************************
How to use
********************************************************************************
The following shows how to use bspipe. The following are examples used in commands.
bowtie2 - mapping program
hg19    - referece genome
project - project name

Use -par file that have all parameter settings.


################################################################################
1. Generate indexes for the reference sequences
################################################################################
shell> bspipe ind -mc bowtie2.conf -i hg19.fa -o index_dir

For RRBS, add -bn mspi -bs C-CGG

################################################################################
2. Generate a sample configuration file
################################################################################
shell> bspipe sam -i read_dir -o sample.conf

To automatically assign group names, use -regex file. Otherwise, modify group names.
--------------------------------------------------------------------------------
demo/sam2grp
--------------------------------------------------------------------------------
MU	Cancer
MCG	Normal
--------------------------------------------------------------------------------

To put the absolute path of FASTQ files into the configuration file, use -full.


################################################################################
3. MAPPING
################################################################################
shell> bspipe map -mc bowtie2.conf -rc index_dir/hg19.conf -sc sample.conf -o mapping_dir

If -full was used, use -sb read_dir to specify the directory of FASTQC files.

For RRBS, specify -bn mspi -bc rrbs.conf

For parallel execution with N multi-cores, use -t N
For parallel execution with N jobs in cluster, use -t N -gc sge.conf

If other patterns that CpG and CpH, use -nt.
For example) -nt CG:CHG:CHH   for CpG, CpHpG, and CpHpH

For cleaning, use -rlen, -qual5, -qual3, -trim5, and/or -trim3.

For high quality of cytosine methylation, adjust -bq, -mq, -e, and/or -clip.


################################################################################
4.a Compile DNA methyaltion of samples to one file
################################################################################
shell> bspipe met -sc sample.conf -rc index_dir/hg19.conf -i mapping_dir -p project -o cpg_dir

################################################################################
4.b Apply a sliding window apprach
################################################################################
shell> bspipe win -sc sample.conf -rc index_dir/hg19.conf -i mapping_dir -w WWW -st SSS -p project -o win_dir
WWW - windows size
SSS - step size

################################################################################
4.c Summarize DNA methylation based on features
################################################################################
shell> bspipe feat -sc sample.conf -rc index_dir/hg19.conf -fc hg19ucsc.conf -i mapping_dir -p project -o feat_dir

To apply a slinding window, use -w and -st.


################################################################################
4.d Use common parameters for 4.a, 4.b, and 4.c.
################################################################################
To use cytosines with N or more reads, use -cov N.

To filter out cytosines in MspI sites for RRBS, specify -bn mspi.

For other pattern than CpG, use -mode.
For example) -mode ch   for CpH


################################################################################
5. Measure sample correlation
################################################################################
shell> bspipe cor -sc sample.conf -i cpg_dir/project.cg.methy.gz -o cor_dir


################################################################################
5. Principle component analysis
################################################################################
shell> bspipe pca -sc sample.conf -i cpg_dir/project.cg.methy.gz -o cor_dir


################################################################################
6. Unsupervised clustering
################################################################################
shell> bspipe clu -sc sample.conf -i cpg_dir/project.cg.methy.gz -o cor_dir


################################################################################
7. Differential methylation
################################################################################
shell> bspipe dmr -sc sample.conf -i cpg_dir/project.cg.methy.gz     -o dmr_dir
shell> bspipe dmr -sc sample.conf -i wind_dir/project.cgwin.methy.gz -o dmr_dir
shell> bspipe dmr -sc sample.conf -i feat_dir/project.cgcgi.methy.gz -o dmr_dir

Without the following parameters, student T-test will be performed.
To perform Wilcoxon rank sum test,       add -wilcoxon.
To perform Kolmogorov–Smirnov test,      add -kstest.
To perform Rao-Scott chi-square test,    add -raoscott.
To perform ANOVA test,                   add -anova.
To perform pairwise Fisher's exaxt test, add -fisher.

For Rao-Scott and Fisher's exact test, use *.nread.gz instead of *.methy.gz.
e.g.) bspipe dmr -sc sample.conf -i cpg_dir/project.cg.nread.gz -o dmr_dir

To specify pairs for comparison, use -pc file.
--------------------------------------------------------------------------------
demo/pair.conf
--------------------------------------------------------------------------------
Cancer	Normal
--------------------------------------------------------------------------------

To filter sites or windows, use -cw.
To adjust signficancy, use -md, -pv and -fdr.
For heatmap, use -dist, -link, -width, and -height.

################################################################################
8. Annotating differentially methylated regions
################################################################################
shell> bspipe ann -ac hg19ann.conf -o ann_dir dmr_dir/project.cgwin.methy.ttest.cancer-normal.q0.05d0.25.csv

For hyper-methylated regions, specify -hyper N.N.
For hypo-methylated regions, specify -hypo  N.N.


################################################################################
9. Functional analysis using DAVID
################################################################################
shell> bspipe david -acc account -idtype GENBANK_ACCESSION -o david_dir ann_dir/project.cgwin.methy.ttest.cancer-normal.q0.05d0.25.hyper.refseq.list.5end

To see the supported ID types, run bspipe david -showid.
To see the supported databases, run bspipe david -showdb.

To adjust figures, use -nitem and -dval.


################################################################################
10. Generate input files for visualization
################################################################################
shell> bspipe vis -rc index_dir/hg19.conf -o vis_dir -sc sample.conf -i mapping_dir

To process BAM, specify -bam.
To generate wiggle files for read coverage, specify -cwig.
To process BAM and BigWig files, specify URL using -url http:// or file that includes the URL.

To use cytosines with N or more reads, use -cov N.
To generate BED files for methylation, use -mbed with cytosine patterns.
To generate wiggle files for methylation, use -mwig with cytosine patterns..
e.g.) -mbed cg   or   -mwig cg,ch

To process specified files, list files instead of -sc sample.conf -i mapping_dir.


################################################################################
11. Box plot
################################################################################
shell> bspipe box -i cpg_dir/project.cg.methy.gz -tar target1.txt -o box_dir
--------------------------------------------------------------------------------
target1.txt
--------------------------------------------------------------------------------
chr19	241776	241906	target1
chr19	246145	246540	target2
chr19	359277	361282	target3

If IDs are in the first column, use -id.
shell> bspipe box -i cpg_dir/project.cg.methy.gz -tar target2.txt -o box_dir -id
--------------------------------------------------------------------------------
target2.txt
--------------------------------------------------------------------------------
target1 chr19:241776-241906
target2 chr19:246145-246540
target3 chr19:359277-361282
--------------------------------------------------------------------------------

To draw stranded methylation, use -strand.
To show sites with no reads, use -show na.
To show positions, use -show pos.
To show the number of sites, use -show num.
To show legend, use -show met.


################################################################################
12. Quality control of BS-seq
################################################################################
shell> bspipe qc -rc index_dir/hg19.conf -o qc_dir -sc sample.conf -i mapping_dir

To process specified files, list files instead of -sc sample.conf -i mapping_dir.

