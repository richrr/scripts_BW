use strict;
use POSIX qw(strftime);
# perl -MGetopt::Long -e 'print( "got it\n" )'
use warnings;
use Getopt::Long qw(GetOptions);


# not using this since a swarm job will run this several times and may create conflicts while writing.
# keep track of versions
# print the versions of the software used
sub trackVersionInfo {
	my $now_string = strftime "%a %b %e %H:%M:%S %Y", localtime;
	`echo "$now_string \n\ncutadapt --version" >> log.txt` ;
	`cutadapt --version >> log.txt` ;

	`echo "\nbbtools bbduk --version" >> log.txt` ;
	`bbtools bbduk --version &>> log.txt` ;

	`echo "\nbowtie2 --version" >> log.txt` ;
	`bowtie2 --version >> log.txt` ;

	`echo "\ntophat -v" >> log.txt` ;
	`tophat -v >> log.txt` ;

	`echo "\nsamtools --version" >> log.txt` ;
	`samtools --version >> log.txt` ;

	`echo "\nhtseq-count --help | tail -n 3" >> log.txt` ;
	`htseq-count --help | tail -n 3 >> log.txt` ;
}


# Quality trimming happens before adapter removal http://cutadapt.readthedocs.io/en/latest/guide.html#modifying-reads
# --max-n Discard reads containing more than COUNT N bases.
# -q trim the 5' and 3' ends with two comma-separated cutoffs (e.g. 15,10)
# -m --minimum-length
# my $cutadaptTemplate = "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC -m 60 --max-n=1 -q 20,20 -o output.fastq input.fastq.gz";



my $adpt_fwd = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" ;
my $adpt_rev = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT";


# page 11 http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_truseq/scriptseq-v2-rna-seq/scriptseq-v2-rna-seq-library-prep-guide.pdf
# AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
# AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT

# page 17 http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_truseq/scriptseq-complete/scriptseq-complete-kit-human-mouse-rat-low-input-library-prep-guide.pdf
# AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
# AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT



my $samtoolsBin = "samtools";
my $inputFile;
my $outputFolder = "./RNA-SEQ";
my $tophatProgram = "tophat";
my $threads = 8;


my $version;
my $ucsc;
my $ensembl;
my $human;
my $mouse;

#my %gf_hash;
#my %bt2_hash;

my $gf;
my $bowind;
my $paired;
my $help;

my $lexogen;

#my $servers_cores = "transkingdom_60";  #separate different machines by ';' # e.g. "transkingdom_10;transkingdom_10;samwise_15"
#my $sge_prefix = "SGE_Batch -c ". ." -m 50G -F 100G -P 10 -r log_mln_p1 -q transkingdom -M rodrrich@oregonstate.edu


#my $min_read_length;

GetOptions (
			"help"  => \$help,   # flag.
			"version=s" => \$version,    # string
			"threads=i" => \$threads,    # integer # number of threads or cores
            "file=s"   => \$inputFile,      # string
            "outputFolder=s"   => \$outputFolder,      # string
            "ucsc"  => \$ucsc,   # flag. # UCSC returns gene symbols
            "ensembl"  => \$ensembl,   # flag. # Ensembl returns ensembl ids
            "human"  => \$human,   # flag.
            "mouse"  => \$mouse,   # flag.
			"lexogen"  => \$lexogen,   # flag. #uses cutadapt by default
			#"sge" => \$sge, # submit jobs on sge
			#"servers_cores=s" => \$servers_cores, # split the full job to run on different machines with that many machine cores
			"paired"  => \$paired),   # flag. # default is single
or usage("Invalid commmand line options.");

if($help){ usage("See commmand line options."); }

usage("The organism must be specified.")
   unless (defined $human or defined $mouse);

usage("The namespace must be specified.")
   unless (defined $ucsc or defined $ensembl);

usage("The version must be specified.")
   unless defined $version;

#### to do: make a config file? ####

# usage: perl /data/rodriguesrr/scripts/perl/RNA-seq-pipeline-bw.pl --file ../MM12-july-19/MM12-july-19/MMP12-138637510/FASTQ_Generation_2019-07-19_07_43_04Z-190032047/62EAT-07-19_L003-ds.2ff961d256ff423e9e4a766fee05321a/62EAT-07-19_S62_L003_R1_001.fastq.gz --mouse --ensembl --version GRCm38 --outputFolder mmp12 --lexogen --threads 8


sub usage {
   my $message = $_[0];
   if (defined $message && length $message) {
      $message .= "\n"
         unless $message =~ /\n$/;
   }

   my $command = $0;
   $command =~ s#^.*/##;

   print STDERR (
      $message,
      "usage: $command --file inputFile [--outputFolder outputFolder] [--paired] --ensembl[|--ucsc] --human[|--mouse] [--lexogen] --version version\n" . #[--sge] [--servers_cores machines_cores;machine2_cores]
      "       Input file is required\n" .
			"       Organism is required\n" .
      "       Namespace is required\n" .
      "       Version of genome and annotation is required\n" . #https://genome.ucsc.edu/FAQ/FAQreleases.html
			"          use --ensembl --version GRCm38 or GRCh37\n" .
			"          OR  --ucsc --version mm10 (=GRCm38) or hg19 (=GRCh37) or hg38 (=GRCh38)\n" .
      "       Output folder defaults to ./RNA-SEQ\n" .
	    "	      Default 8 threads\n" .
	    "       Specify --lexogen\n" . #. Defaults to cutadapt.
      "       Defaults to single end.\n" # Specify --paired if needed\n" .
     # "       Use --sge with --servers_cores args to submit and split jobs on sge"
   );

   die("\n")
}


# find these indexes at https://hpc.cit.nih.gov/apps/db.php?f=Igenomes
# to do: add a way to map to gene symbols, or do both.


my $gen_idx_path = "/fdb/igenomes/" ;

if($human){
	$gen_idx_path .= "Homo_sapiens/" ;
} elsif($mouse){
	$gen_idx_path .= "Mus_musculus/" ;
}

if($ucsc){
	$gen_idx_path .= "UCSC/" ;
} elsif($ensembl){
	$gen_idx_path .= "Ensembl/" ;
}

$bowind = $gen_idx_path . $version . "/Sequence/Bowtie2Index/genome" ;
$gf = $gen_idx_path . $version . "/Annotation/Genes/genes.gtf";

#print $gf . "\n" . $bowind . "\n";

if (-e $gf) {
    print "$gf exists\n";
} else {
    usage("the gtf file does not exist!\n");
}

if (glob("$bowind*")) { # At least one file matches "$bowind*"
	print "$bowind exists\n";
} else {
	usage("the index files do not exist!\n");
}


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	#	work flow
	#	1. remove adapters
	#	2. remove low quality sequences(better performed before remove of "N",for some file ,every sequence has an "N" in a particular possition, if we remove "N" first, all sequences will be removed. but these "N"s may be removed in this quality filter step, so after it is removed, we do not have to worry about all sequences to be discarded)
	#	3. remove sequence with any "N"
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 														remove sequence with "N"====record change of file after each step
#																	use prinseq
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

my $filtertool = 'CutAdapt';
if ($lexogen) {$filtertool = 'Lexogen';}
print "\nUsing $filtertool\n";

my $outFolder = "$outputFolder.After$filtertool";
`mkdir $outFolder` if ( ! -d $outFolder);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

if(1==0){
sub cutadapt{
	my ($inputFile, $OutPutFile) = @_;

    my $cutadaptTemplate = "cutadapt -a $adpt_fwd -a $adpt_rev -m 60 --max-n=1 -q 20,20 -o output.fastq input.fastq.gz";
	$cutadaptTemplate =~s/input.fastq.gz/$inputFile/g;
	$cutadaptTemplate =~s/output.fastq/$OutPutFile/g;
	print $cutadaptTemplate,"\n\n"; #die;
	`$cutadaptTemplate` ;
}


# paired end
sub cutadapt_pe{
	my ($inputFileR1, $OutPutFileR1, $inputFileR2, $OutPutFileR2) = @_;

	my $cutadaptTemplate = "cutadapt -a $adpt_fwd -A $adpt_rev -m 60 --max-n=1 -q 20,20 -o out.1.fastq -p out.2.fastq reads.1.fastq.gz reads.2.fastq.gz";

	$cutadaptTemplate =~s/reads.1.fastq.gz/$inputFileR1/g;
	$cutadaptTemplate =~s/out.1.fastq/$OutPutFileR1/g;

	$cutadaptTemplate =~s/reads.2.fastq.gz/$inputFileR2/g;
	$cutadaptTemplate =~s/out.2.fastq/$OutPutFileR2/g;

	print $cutadaptTemplate,"\n\n"; #die;
	`$cutadaptTemplate` ;
}
}
#use sequencePreprocessing::FileHandle;


sub lexogen_cleaning{
	my ($inputFile, $OutPutFile) = @_;
	my $cmd="gunzip < $inputFile > $OutPutFile";
	`$cmd` ;

	my $cleanOutPutFile= "$OutPutFile._ft12trimmed_q15_clean" ;
	my $bshcmd="bbtools bbduk in=$OutPutFile out=$cleanOutPutFile ref=/data/rodriguesrr/bbtools_bbmap/resources/adapters.fa k=13 ktrim=r forcetrimleft=12 useshortkmers=t mink=5 qtrim=r trimq=15 minlength=20";
	`$bshcmd` ;
	return $cleanOutPutFile;
}



sub processFastq{

			my $fileNameR2 = ' ';
			my $InputFastqFileR1;
			my $OutPutFileR1 ;
			my $InputFastqFileR2;
			my $OutPutFileR2;
			my $cutAdaptorFileR2;

			# split the input filename to use the last string as the file name and the remaining is the folder
			my @spl = split('/', $inputFile);

			my $inputFolder = join('/', @spl[0..$#spl-1]);
			chomp $inputFolder;
	    my $fileNameR1 = $spl[-1] ;
			chomp $fileNameR1;

			#if($paired)	{
			#	$fileNameR2 = $fileNameR1;
			#	$fileNameR2 =~ s/_R1_001/_R2_001/; #$gzFilesR2[$i] ;
			#	chomp $fileNameR2 ;
			#}
			print "Input folder $inputFolder\n";
			print "File $fileNameR1 $fileNameR2\n"; # next;



			my $fileName = "$fileNameR1";
			#if($paired)	{
			#	$fileName = "$fileNameR1-$fileNameR2";
			#}

			$fileName =~s/_001.fastq.gz//g;


			if($lexogen){
				# do the uncompressing and suggested cleaning
				$InputFastqFileR1 = "$inputFolder/$fileNameR1";
				#$OutPutFileR1 = $InputFastqFileR1 ;
				$OutPutFileR1 = "$outFolder/$fileNameR1.$filtertool";
				$OutPutFileR1 = lexogen_cleaning($InputFastqFileR1, $OutPutFileR1);
			} #else {
				#$InputFastqFileR1 = "$inputFolder/$fileNameR1";
				#$OutPutFileR1 = "$outFolder/$fileNameR1.$filtertool";

				#if($paired)	{
				#	$InputFastqFileR2 = "$inputFolder/$fileNameR2";
				#	$OutPutFileR2 = "$outFolder/$fileNameR2.$filtertool";
				#	cutadapt_pe($InputFastqFileR1, $OutPutFileR1, $InputFastqFileR2, $OutPutFileR2) if ((! -e "$OutPutFileR1") && (! -e "$OutPutFileR2"));
				#} else {
				#	# run single end cutadapt
				#	cutadapt($InputFastqFileR1, $OutPutFileR1) if (! -e "$OutPutFileR1");
				#}
			#}

			my $cutAdaptorFileR1 = $OutPutFileR1;
			#if($paired)	{
			#	$cutAdaptorFileR2 = $OutPutFileR2;
			#}
			my $tophatOutFolder = "$outFolder.tophat/$fileName";
			`mkdir $outFolder.tophat` if (! -d "$outFolder.tophat");
			#if($paired)	{
			#	runTopHat_pe($cutAdaptorFileR1,$cutAdaptorFileR2,$tophatOutFolder) if (! -e "$tophatOutFolder/accepted_hits.bam");
			#} else {
				# run single end tophat
				runTopHat($cutAdaptorFileR1,$tophatOutFolder) if (! -e "$tophatOutFolder/accepted_hits.bam");
			#}


			my $samFolder = "$outFolder.tophat.Aligned.sam";
			`mkdir $samFolder` if (! -d $samFolder);
			my $acceptedBam = "$tophatOutFolder/accepted_hits.bam";
			# sort by name, convert to SAM for htseq-count
			print "$samtoolsBin sort -n --threads $threads -o $acceptedBam.sn.bam --output-fmt BAM $acceptedBam  \n";
			`$samtoolsBin sort -n --threads $threads -o $acceptedBam.sn.bam --output-fmt bam $acceptedBam ` if (! -e "$samFolder/$fileName");
			`$samtoolsBin view -h --threads $threads --output-fmt sam -o $samFolder/$fileName $acceptedBam.sn.bam` if (! -e "$samFolder/$fileName");
			`$samtoolsBin view --threads $threads --output-fmt sam -o $samFolder/view_wo-h/$fileName $acceptedBam.sn.bam` if (! -e "$samFolder/view_wo-h/$fileName");
# samtools view -?
# 3. SAM->BAM conversion: `samtools view -b in.sam.gz'.    -o FILE  output file name
#  4. BAM->SAM conversion: `samtools view -h in.bam'.  # do you need -h?



			my $countFolder = "$outFolder.tophat.Aligned.sam.htseqCount";
			`mkdir $countFolder` if (! -d $countFolder);
			#my $htseqCommand = "htseq-count --order=name -a 10 -f sam --stranded=no -o $countFolder/$fileName $samFolder/$fileName $gf"; # outputs SAM file
			my $htseqCommand = "htseq-count --order=name -a 10 -f sam --stranded=no $samFolder/$fileName $gf > $countFolder/$fileName.txt";

			print "\n\n$htseqCommand\n\n\n"; #die;
			`$htseqCommand\n` if (! -e "$countFolder/$fileName");


		#}
}





#sub runTopHat_pe {
#	my ($cutAdaptorFileR1, $cutAdaptorFileR2, $tophatOutFolder) = @_;
#	my $cmd = "$tophatProgram --no-coverage-search -G $gf -p $threads -o $tophatOutFolder $bowind $cutAdaptorFileR1 $cutAdaptorFileR2 ";
#	print ($cmd,"\n\n");
#	`$cmd`;
#}


sub runTopHat {
	my ($cutAdaptorFile, $tophatOutFolder) = @_;
	my $cmd = "$tophatProgram --no-coverage-search -G $gf -p $threads -o $tophatOutFolder $bowind $cutAdaptorFile";
	print ($cmd,"\n\n");
	`$cmd`;
}

# coverage search on
#sub runTopHat {
#	my ($cutAdaptorFile, $tophatOutFolder) = @_;
#	my $cmd = "$tophatProgram --coverage-search -G $gf -p 8 -o $tophatOutFolder $bowind $cutAdaptorFile";
#	print ($cmd,"\n\n");
#	`$cmd`;
#}

#trackVersionInfo;
processFastq;
