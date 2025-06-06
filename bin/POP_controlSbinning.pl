#!/usr/bin/perl -w

###############################################################################
# Script Name: gatk_pipeline_generator.pl
# Author: [Your Name]
# Description:
#   This script generates GATK-based variant calling pipeline shell scripts
#   for a set of DNA sequencing samples using provided cleaned data paths,
#   reference genomes, and configuration parameters.
#
# Purpose:
#   Automates genome profiling by creating:
#     - Shell scripts to run alignment (bowtie2), sorting (Picard), and coverage analysis
#     - Output directory structure for each sample-reference pair
#     - Parallel run support by splitting jobs based on configured thread counts
#
# Required Inputs:
#   -cl      <cleaned_data.txt>     : 2-column tab-separated file: sample_id \t fq1_path \t fq2_path
#   -R       <reference_list.txt>   : Tab-separated list of reference genome name and path (e.g., hg38\t/path/hg38.fa)
#   -cfg     <config.txt>           : Tab-delimited configuration file for tool paths and filtering thresholds
#
# Optional Inputs:
#   -ot_dir  <output_directory>     : Output base directory (default: current working directory)
#   -sh      <shell_prefix>         : Prefix for shell script names (default: gatk)
#
# Output Structure:
#   <output_dir>/
#     ├── Shell/               # All per-sample-per-reference .sh scripts
#     ├── gatk/                # Output directories for each sample/reference
#     ├── Main.gatk.sh         # Main runner that includes all shell script calls
#     ├── Main.gatk.FLAG*sh    # Splits of the main script for parallel runs
#
# Dependencies:
#   - Perl 5+
#   - Java (for Picard/GATK)
#   - samtools, bedtools, bowtie2, GATK, Picard installed and in $PATH
#
# Usage Example:
#   perl gatk_pipeline_generator.pl -cl clean_data.txt -R reference.txt -cfg config.txt -ot_dir ./output -sh gatk
#
###############################################################################
#!/usr/bin/perl -w
use strict;
use Cwd 'abs_path';
use File::Basename;
use Getopt::Long;


my %opt = qw();
GetOptions(\%opt,"cl:s","R:s","cfg:s","ot_dir:s","sh:s");
my $help = " This script is to generate shells for genome profiling
-cl     cleandatapath, 2 column file
-R	reference fasta file: fasName	fasFile
-cfg    configureation file
-ot_dir output directory, default is current directory
	./Shell
	./GATK
-sh     shell name, default is gatk
";
if (scalar keys %opt == 0){
print "$help";
exit;
}
#
#
my @error;
my $cleanpath;
my %reference;
my $cfg;
my $ot_dir;
my $shellname;
&CFG;
my $bowtie2thread;
my $gatkJar;
my $picardJar;
my $shell_splitNum;
my $process_number;
my $snpFilterExpression;
my $indelFilterExpression;
my $adRatio;
my $RGPL;
&PARAMETER;
#
my %sample_data;
my %data_fq2;
my %data_fq1;
&READ_CLEAN_PATH;
#generate shells
&MKSHELL;
#spit run
&SPLIT_RUN;
sub CFG{
	#clean data path
	if (exists $opt{"cl"}){
	$cleanpath = $opt{"cl"};
	push (@error,"cl:$cleanpath not exists!") unless (-e $cleanpath);
	}
	else{
	push (@error,"-cl must be intialized");
	}
	#reference
	if (exists $opt{"R"}){
		open R,$opt{"R"} or die "$opt{R} not E"; <R>;
		while (<R>){
		chomp;
		my @ar = split /\s+/,$_;
		$reference{$ar[0]} = $ar[1];
			push (@error,"R:$ar[1] not exists!") unless (-e $ar[1]);
			system("samtools faidx $ar[1]") unless -e "$ar[1].fai";
		       	my @R = split /\./,$ar[1];
		        system("gatk CreateSequenceDictionary -R $ar[1] --VERBOSITY ERROR") unless -e "$R[0].dict";
			system("bowtie2-build $ar[1] $ar[1]") unless -e "$ar[1].1.bt2";
		}
		close (R);
	}
	else{
	push (@error,"-R must be intialized");
	}
	#output directory
	if (exists $opt{"ot_dir"}){
#	$ot_dir = abs_path($opt{"ot_dir"});
	$ot_dir = $opt{"ot_dir"};
	`mkdir $ot_dir` unless (-e $ot_dir);
	}
	else{
	$ot_dir = abs_path("./");
	}
	#cfg
	if (exists $opt{"cfg"}){
	$cfg = $opt{"cfg"};
	push (@error,"cfg:$cfg not exists!") unless (-e $cfg);
	}
	else{
	push (@error,"-cfg must be intialized");
	}
	#shell name
	if (exists $opt{"sh"}){
	$shellname = $opt{"sh"};
	}
	else{
	$shellname = "gatk";
	}
if (scalar @error > 0){
my $error = join ("\n",@error);
die "$error\n######\n$help";
}

}

sub PARAMETER{
	my @filter_snp = ();
	my @filter_indel = ();
	open CF,"$cfg" or die "cfg:$cfg not E!";
	while (<CF>){
	chomp;
	my @cfg = split /\t/,$_;
	$bowtie2thread = $cfg[1] if $cfg[0] eq "bowtie2thread";
	$picardJar = $cfg[1] if $cfg[0] eq "picardJar";
	$gatkJar = $cfg[1] if $cfg[0] eq "gatkJar";
	#QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0
	push (@filter_snp,"QD < $cfg[1]") if $cfg[0] eq "SNP::QD";
	push (@filter_snp,"FS > $cfg[1]") if $cfg[0] eq "SNP::FS";
	push (@filter_snp,"MQ < $cfg[1]") if $cfg[0] eq "SNP::MQ";
	push (@filter_snp,"MQRankSum < $cfg[1]") if $cfg[0] eq "SNP::MQRankSum";
	push (@filter_snp,"ReadPosRankSum < $cfg[1]") if $cfg[0] eq "SNP::ReadPosRankSum";
	push (@filter_snp,"SOR > $cfg[1]") if $cfg[0] eq "SNP::SOR";
	push (@filter_snp,"DP < $cfg[1]") if $cfg[0] eq "SNP::DP";
	$adRatio = $cfg[1] if $cfg[0] eq "SNP::AD/DP";
	#QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0
	push (@filter_indel,"QD < $cfg[1]") if $cfg[0] eq "INDEL::QD";
	push (@filter_indel,"FS > $cfg[1]") if $cfg[0] eq "INDEL::FS";
	push (@filter_indel,"ReadPosRankSum < $cfg[1]") if $cfg[0] eq "INDEL::ReadPosRankSum";
	push (@filter_indel,"SOR > $cfg[1]") if $cfg[0] eq "INDEL::SOR";
	$process_number = $cfg[1] if $cfg[0] eq "process_number";
	$RGPL = $cfg[1] if $cfg[0] eq "RGPL";
	}
	close (CF);
$snpFilterExpression = join ("||",@filter_snp);
$indelFilterExpression = join ("||",@filter_indel);
}

sub READ_CLEAN_PATH{
	open CL,"$cleanpath"; <CL>;
	while (<CL>){
	chomp;
	my @array = split /\t/,$_;
	my $key = "$array[0]";
		if (exists $data_fq1{$key}){
		$data_fq1{$key} = "$data_fq1{$key},$array[1]";
		}
		else{
		$data_fq1{$key} = $array[1];
		}
		if (exists $data_fq2{$key}){
		$data_fq2{$key} = "$data_fq2{$key},$array[2]";
		}
		else{
		$data_fq2{$key} = $array[2];
		} 

	#	
		if (exists $sample_data{$array[0]}){
		$sample_data{$array[0]} = "$key\n$sample_data{$array[0]}";
		}
		else{
		$sample_data{$array[0]} = $key;
		}
	}
}

sub MKSHELL{
`mkdir $ot_dir/Shell`;
`mkdir $ot_dir/gatk`;
open MAIN,">$ot_dir/Main.gatk.sh" or die "Main.gatk.sh";
	foreach my $m (keys %sample_data){#/home/lia/flow/mPD/09_GATK_FLOW/gatk.parameter
	#mk sam for different inout fq
	my @files = split /\n/,$sample_data{$m}; my @fq1; my @fq2;
	foreach my $n (@files){push (@fq1,$data_fq1{$n});} foreach my $n (@files){push (@fq2,$data_fq2{$n});}
	my $fq1 = join (",",@fq1); 
	my $fq2 = join (",",@fq2);

	my $files_input =  "NA";
	   if ($fq2 eq "NA"){
           $files_input = "-U $fq1";
	   }
	   else{
	   $files_input = "-1 $fq1 -2 $fq2";
	   }
		foreach my $n (keys %reference){
		`mkdir $ot_dir/gatk/$m\_$n/`;
		my $reference = $reference{$n};
		my $otPrefix = "$ot_dir/gatk/$m\_$n/$m\_$n"; 
		#@ST-E00600:105:H75L3ALXX:1:1101:8332:1016 1:N:0:AGTCACTA
		#@S250077624L1C001R00100004416

		open SHELL,">$ot_dir/Shell/$m\_$n.sh" or die;
		print MAIN "sh $ot_dir/Shell/$m\_$n.sh > $ot_dir/Shell/$m\_$n.log\n";
#		my $head = `gzip -dc $fq1[0] | head -1`; my @head = split /\:/,$head; 
		my $machine = "HeNanGeneHospital"; my $lane = 1; my $flow= "Flow";
		#step0 build gatk index
		#step1 bowtie2 to generate sam
		print SHELL "#step1 bowtie2 generate sequences alignments map (SAM)\n";
		print SHELL "nohup bowtie2 -x $reference $files_input -S $otPrefix.step01.bowtie2.sam --no-unal -p $bowtie2thread > $otPrefix.step01.bowtie2Run.log\n";
		#step2 generate &sort bam
		print SHELL "#step2 Generate bam then sort\n";
		print SHELL "nohup java -jar $picardJar SortSam -INPUT $otPrefix.step01.bowtie2.sam -OUTPUT $otPrefix.step02.sort.bowtie2.sam -SORT_ORDER coordinate  -VALIDATION_STRINGENCY LENIENT > $otPrefix.step02.samSort.log \n";
		#step3 align matrics
		print SHELL "#step3  Collect Alignment & Insert Size Metrics\n";
		print SHELL "nohup java -jar $picardJar CollectAlignmentSummaryMetrics -R $reference -I $otPrefix.step02.sort.bowtie2.sam -O $otPrefix.step03.alignmengt_metrics.txt --VALIDATION_STRINGENCY LENIENT > $otPrefix.step03.collectSamMatrix.log\n";
		#print SHELL "samtools depth -a $otPrefix.step02.sort.bowtie2.sam > $otPrefix.step03.depth_out.txt\n";
		#step19 Compute Coverage Statistics
		print SHELL "samtools view  $otPrefix.step02.sort.bowtie2.sam -o $otPrefix.step.02.bam -O BAM\n";
		print SHELL "bedtools genomecov -bga -ibam  $otPrefix.step.02.bam > $otPrefix.genomecov.bedgraph\n";
		print SHELL "echo \'$m gatk is finished\' \n";
		close (SHELL);
		}
	}
close (MAIN);
}

sub SPLIT_RUN{
#$ot_dir/Main.gatk.sh
#
my $sample_number = `cat $ot_dir/Main.gatk.sh | wc -l`;
   $sample_number =~ s/\n//;
my $number = int ( $sample_number / $process_number );
$number = 1 if $number == 0;
`split -d -l $number  $ot_dir/Main.gatk.sh  $ot_dir/Main.gatk.FLAG`;
my $L = `ls $ot_dir/Main.gatk.FLAG*`;
my @L = split /\n/,$L;
	foreach my $m (@L){
	`mv $m $m.sh`;
#	`nohup sh $m.sh > $m.error &`;
	}
#`rm $ot_dir/Main.gatk.sh`;
#while (1){
#sleep(10);
#my $finishedNumber = `grep 'gatk is finished' $ot_dir/Shell/*.log | wc -l`;
#last if $finishedNumber == $sample_number;
#}
}

sub DATE{
my $DATE = `date "+%Y-%m-%d %H:%M:%S"`;
   $DATE =~ s/\n//;
return($DATE);
}

