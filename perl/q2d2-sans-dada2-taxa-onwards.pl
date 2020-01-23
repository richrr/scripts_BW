#!/usr/local/bin/perl

# got this on oct 29 2019 from /data/MicrobiomeCore/scripts/q2d2.pl

# the code is basically to get the taxa and the trees after merging the tables and rep-seqs of some dada2 runs.

use strict;
use vars qw($opt_d $opt_c $opt_p $opt_m);
use Cwd;
use Getopt::Std;
use FindBin ();

getopts("d:c:p:m:");
my $usage = "\nUSAGE:
perl $0
	-d <name of input and output dir>
	-m <name of mapping file>
	[-p <number of processors> default=30]
	[-c <config file with configurable parameters> default=$FindBin::Bin/q2d2_default.config]


	This script will do the following:
	- creates output directory using provided name
	- use the config file for tweakable parameters for the various steps in the process

	Please note the following:
	- '$FindBin::Bin/q2d2_default.config' is used as the default config file
	- due to module dependencies; this script currently only functions when user is logged into the NIH Biowulf HPC system

Example:
perl $0
		-d /scratch/thovaraivv/test_out
		-m /scratch/thovaraivv/mapping.txt
		-c /data/MicrobiomeCore/scripts/q2d2_default.config
		-p 30
";

if(!(
	defined($opt_d) &&
	defined($opt_m))){
		print STDERR $usage;
		exit();
}

if(!(defined($opt_c))){
		$opt_c = "$FindBin::Bin/q2d2_default.config";
}

if(!(defined($opt_p))){
		$opt_p = 30;
}

######################
# Assigning variables


my $out_dir=$opt_d."/";
my $config_file=$opt_c;
my $map_file=$opt_m;
my $commands="commands.sh";
my $metrics_dir="metrics";
my $config_hash_ref = readConfig($config_file);


###########
# Set up

# check if output directory already exists
#if(-e $out_dir){
#    print STDERR "\nThe following output directory already exists.\n\t$out_dir\n";
#    print STDERR "Please enter a new directory name to avoid overwriting existing data.\nQuitting . . .\n";
#    die;
#}
#create_dir($out_dir);
chdir $out_dir or die "Can't chdir to $out_dir:$!\n" if $out_dir;
#create_dir("metrics");

print "\nWRITING TO:\n\t$out_dir\n";
print "\nUSING METADATA FILE:\n\t$map_file\n";
print "\nUSING CONFIG FILE:\n\t$config_file\n";


##########################
# build shell command file

# load modules
my $commands_string = "#!/bin/bash";
$commands_string .= "\nexport TMPDIR=/lscratch/\$SLURM_JOB_ID\n\nmodule unload qiime
					module load qiime/2-2019.4";


# qiime2 commands
$commands_string .= "\n\nqiime feature-table summarize \\
  --i-table table.qza \\
  --o-visualization table.qzv \\
  --m-sample-metadata-file $map_file
qiime feature-table tabulate-seqs \\
  --verbose \\
  --i-data rep-seqs.qza \\
  --o-visualization rep-seqs.qzv
qiime feature-classifier classify-consensus-vsearch \\
  --i-query rep-seqs.qza \\
  --i-reference-reads $config_hash_ref->{'ref_seq'} \\
  --i-reference-taxonomy $config_hash_ref->{'ref_tax'} \\
  --o-classification taxonomy.qza
qiime metadata tabulate \\
  --m-input-file taxonomy.qza \\
  --o-visualization taxonomy.qzv
qiime taxa barplot \\
  --i-table table.qza \\
  --i-taxonomy taxonomy.qza \\
  --m-metadata-file $map_file \\
  --o-visualization taxa-bar-plots.qzv
qiime alignment mafft \\
  --i-sequences rep-seqs.qza \\
  --o-alignment aligned-rep-seqs.qza
qiime alignment mask \\
  --i-alignment aligned-rep-seqs.qza \\
  --o-masked-alignment masked-aligned-rep-seqs.qza
qiime phylogeny fasttree \\
  --i-alignment masked-aligned-rep-seqs.qza \\
  --o-tree unrooted-tree.qza
qiime phylogeny midpoint-root \\
  --i-tree unrooted-tree.qza \\
  --o-rooted-tree rooted-tree.qza
qiime diversity alpha-rarefaction \\
  --i-table table.qza \\
  --i-phylogeny rooted-tree.qza \\
  --p-max-depth $config_hash_ref->{'rarefy_to'} \\
  --m-metadata-file $map_file \\
  --o-visualization prelim_alpha_rarefaction.qzv";

$commands_string=~s/[\t]+//g;

##########################
# run command shell file

open(CMD, ">$commands") || die "Could not open $commands ($!)\n";
print CMD "$commands_string";
close (CMD);
chmod 0755, $commands;


##########################
# set up sbatch job

my $cmd_str = "sbatch --cpus-per-task=$opt_p --time=20:00:00 --mail-type=END $commands";
print "\nLAUNCHING COMMAND SCRIPT ON $opt_p PROCESSORS:\n\t$out_dir$commands\n";
exec_cmd($cmd_str);


##############################################################################
# executes command string in unix shell

sub exec_cmd{
	my $cmd=shift;
	$cmd=~s/[\t]+/ /g;
	print STDERR "EXECUTING:\n\t$cmd\n";
	system($cmd);
}


##############################################################################
# creates dir

sub create_dir{

    my $dir=shift;

    if(-e $dir){
	    print STDOUT "Directory $dir already exists. Using it.\n";
    }else{
	    mkdir $dir;
	    if(-e $dir){
		    print STDERR "$dir has been created.\n";
	    }else{
		    die "Could not create $dir.\n";
	    }
    }
}

##############################################################################
# reads config file

sub readConfig{
    my $config = shift;
    my %config_hash;
	my $comment = "#";
    open(CF, "<$config") || die "Could not open $config($!)\n";
    while (<CF>){
		if((!/^$comment/) && (/(\S+)=(\S+)/)){
			$config_hash{$1} = $2;
		}
    }
    return \%config_hash;
}
