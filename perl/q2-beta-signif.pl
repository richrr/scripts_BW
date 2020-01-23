#!/usr/local/bin/perl

# got this on oct 29 2019 from /data/MicrobiomeCore/scripts/q2d2.pl

# the code is basically to get the taxa and the trees after merging the tables and rep-seqs of some dada2 runs.

use strict;
use vars qw($opt_d $opt_c $opt_p $opt_m $opt_b);
use Cwd;
use Getopt::Std;
use FindBin ();

getopts("d:c:p:m:b:");
my $usage = "\nUSAGE:
perl $0
	-d <name of input and output dir>
	-m <name of mapping file; arg for --m-metadata-file>
	-c <name of metadata column; arg for --m-metadata-column>
	-p <comma separated list of methods to use ('permanova,anosim'); arg for --p-method>
	-b <comma separated list of beta diversity distance matrices to use ('bray_curtis,jaccard,unweighted_unifrac,weighted_unifrac')>


	This script will do the following:
	- Determine whether groups of samples are significantly different from one
  another

	Please note the following:
	- due to module dependencies; this script currently only functions when user is logged into the NIH Biowulf HPC system

Example:
perl $0
		-d /scratch/thovaraivv/test_out
		-m /scratch/thovaraivv/mapping.txt
		-c group
		-p permanova,anosim
		-b bray_curtis,jaccard,unweighted_unifrac,weighted_unifrac
";

if(!(
	defined($opt_d) &&
	defined($opt_m))){
		print STDERR $usage;
		exit();
}

if(!(defined($opt_c))){
		print STDERR $usage;
		exit();
}

if(!(defined($opt_p))){
		$opt_p = "permanova,anosim";
}

if(!(defined($opt_b))){
		$opt_b = "bray_curtis,jaccard,unweighted_unifrac,weighted_unifrac";
}

######################
# Assigning variables


my $out_dir=$opt_d."/";
my $map_file=$opt_m;
my $time_stamp = get_time_stamp();
my $commands="commands.bsig.".$opt_c."_".$time_stamp.".sh";
my $vis_dir=$out_dir."beta_sig_".$opt_c."_".$time_stamp."/";


###########
# Set up

chdir $out_dir or die "Can't chdir to $out_dir:$!\n" if $out_dir;
create_dir($vis_dir);

print "\nWRITING TO:\n\t$vis_dir\n";
print "\nUSING METADATA FILE:\n\t$map_file\n";



##########################
# build shell command file

# load modules
my $commands_string = "#!/bin/bash";
$commands_string .= "\nexport TMPDIR=/lscratch/\$SLURM_JOB_ID\n\nmodule unload qiime
					module load qiime/2-2019.4";


my @betadists = split /,/, $opt_b;
my @methods = split /,/, $opt_p;

print "\nRunning SEVERAL significance calculations:\n";
foreach my $b (@betadists) {
  print $b,"\n";

	my $distmat = $out_dir.$b."_distance_matrix.qza" ;

	foreach my $p (@methods) {
	  print "    -",$p,"\n";

		my $outpfile = $vis_dir.$b."_".$p."_bsig.qzv" ;

		# qiime2 commands
		$commands_string .= "\n\nqiime diversity beta-group-significance \\
		--i-distance-matrix $distmat \\
		--m-metadata-file $map_file \\
		--m-metadata-column $opt_c \\
		--p-method $p \\
		--p-pairwise \\
		--o-visualization $outpfile
		";


	}

}

#print "$commands_string";
$commands_string=~s/[\t]+//g;

##########################
# run command shell file

open(CMD, ">$commands") || die "Could not open $commands ($!)\n";
print CMD "$commands_string";
close (CMD);
chmod 0755, $commands;


##########################
# set up sbatch job

my $cmd_str = "sbatch --cpus-per-task=2 --time=10:00:00 --mail-type=END $commands";
print "\nLAUNCHING COMMAND SCRIPT ON 2 PROCESSORS:\n\t$out_dir$commands\n";
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
# gets timestamp
sub get_time_stamp{
	my($day, $month, $year, $hour, $min, $sec) = (localtime)[3,4,5,2,1,0];
	$month = sprintf '%02d', $month+1;
	$day   = sprintf '%02d', $day;
	my $ds= join ('', $year+1900, $month, $day, $hour, $min, $sec);
	print "$ds\n";
	return($ds);
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
