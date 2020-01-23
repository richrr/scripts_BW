#!/usr/local/bin/perl

use strict;
use vars qw($opt_t $opt_s $opt_o $opt_m $opt_p $opt_e $opt_h);
use Cwd;
use Getopt::Std;
use FindBin ();
use List::Util qw(first);

# usage:
# split examples:
# perl q2sm.pl -t p1/table.qza -s p1/rep-seqs.qza -o somedir -m p1/p1-map.txt -p SAMPLE_TYPE,ALL
# perl q2sm.pl -t p1/table.qza -s p1/rep-seqs.qza -o somedir -m p1/p1-map.txt -p SAMPLE_TYPE,stool
# perl q2sm.pl -t p1/table.qza -o somedir -m p1/p1-map.txt -p SAMPLE_TYPE,ALL

# merge examples:
# perl q2sm.pl -t somedir/Oral_swab_culture.table.qza,somedir/stool.table.qza -s somedir/Oral_swab.rep-seqs.qza,somedir/stool.rep-seqs.qza -o somedir.m -m p1/p1-map.txt -e stool.oral.swab
# perl q2sm.pl -t somedir/Oral_swab_culture.table.qza,somedir/stool.table.qza -o somedir.m -m p1/p1-map.txt -e t.stool.oral.swab

getopts("t:s:o:m:p:e:h");
my $usage = "\nUSAGE:
perl $0
	-t <comma separated list of tables to be merged. Give only one table.qza for split>
	-s <comma separated list of sequences to be merged. Give only one sequence.qza for split>
	-o <name of output dir>
	-m <name of mapping file>
	-p <split using COLUMN,KEYWORD from metadata.
			COLUMN,ALL will split per unique entry (group) in the coulmn.
			COLUMN,KEYWORD will keep only the stuff that exactly matches the keyword. use if you want something specific>
	-e <name for the merged tables/seqs>

	This script needs:
	-> at least one of -t (table.qza) and -s (rep-seqs.qza)
	-> at least one of -p (split) and -e (merge)
			=> -p (split) requires -t (table.qza) with or without -s (rep-seqs.qza)
			=> -e (merge) can take either -t, -s, or both
	-> both -o (output) and -m (map file)


	This script can do the following:
	- split/merge tables/sequences
	- creates output directory using provided name

Example to build command:
	# to split
		perl $0 -t table.qza -s rep-seqs.qza -m map.txt -p Column,Group -o outdir
	# to merge
		perl $0 -t table1.qza,table2.qza -s rep-seqs1.qza,rep-seqs2.qza -o outdir -m map.txt -e outputPrefix
";


if(defined($opt_h)){
	print STDERR $usage;
	exit();
}

# needs mapping file and output directory
if(!(defined($opt_o) &&
	defined($opt_m))){
		print STDERR $usage;
		exit();
}

# needs table or sequence argument
if(!(
	defined($opt_t) ||
	defined($opt_s))){
		print STDERR $usage;
		exit();
}

##### figure this out ####
# needs to say whether split or merge
if((!defined($opt_p)) && (!defined($opt_e))){
		print STDERR $usage;
		exit();
}

######################
# Assigning variables
my $out_dir=$opt_o."/";
my $map_file=$opt_m;

my $time_stamp = get_time_stamp();
my $commands="commands.".$time_stamp.".";

my @BIG_CMD_ARRAY;



###########
# Set up

# check if output directory already exists
if(-e $out_dir){
    print STDERR "\nThe following output directory already exists.\n\t$out_dir\n";
    print STDERR "Please enter a new directory name to avoid overwriting existing data.\nQuitting . . .\n";
    die;
}
create_dir($out_dir);

print "\nWRITING TO:\n\t$out_dir\n";
print "\nUSING METADATA FILE:\n\t$map_file\n";


## merge ##
if(defined($opt_e)){

	$commands .= "merge.sh";

	## the table and sequence don't have ','
	if (index($opt_t, ',') == -1 && index($opt_s, ',') == -1) {
   print "Need ',' in atleast one of $opt_t (-t) or $opt_s (-s) when using with merge (-e)\n";
	 print STDERR $usage;
	 exit();
	}

	## the table
	if (index($opt_t, ',') != -1) {
	  my @tabs = split /,/, $opt_t;
		my $final_cmd = "qiime feature-table merge";
		foreach my $tbl (@tabs){
			 $final_cmd .= " --i-tables $tbl"
		}

		$final_cmd .= " --o-merged-table $out_dir$opt_e.merged.table.qza";
		print $final_cmd, "\n";
		push @BIG_CMD_ARRAY, $final_cmd;

		my $makeqzv = "qiime feature-table summarize --verbose --i-table $out_dir$opt_e.merged.table.qza --o-visualization $out_dir$opt_e.merged.table.qzv --m-sample-metadata-file $map_file";
		print $makeqzv, "\n\n";
		push @BIG_CMD_ARRAY, $makeqzv;
	}

	## the sequence
	if (index($opt_s, ',') != -1) {
		my @seqs = split /,/, $opt_s;
		my $final_cmd = "qiime feature-table merge-seqs";
		foreach my $seq (@seqs){
			 $final_cmd .= " --i-data $seq"
		}

		$final_cmd .= " --o-merged-data $out_dir$opt_e.merged.rep-seqs.qza";
		print $final_cmd, "\n";
		push @BIG_CMD_ARRAY, $final_cmd;

		my $makeqzv = "qiime feature-table tabulate-seqs --verbose --i-data $out_dir$opt_e.merged.rep-seqs.qza --o-visualization $out_dir$opt_e.merged.rep-seqs.qzv";
		print $makeqzv, "\n\n";
		push @BIG_CMD_ARRAY, $makeqzv;
	}

}



## split ##
if(defined($opt_p)){

	$commands .= "split.sh";

	## the table
	if (index($opt_t, ',') != -1) {
   print "Remove ',' from $opt_t (-t) when using with split (-p)\n";
	 print STDERR $usage;
	 exit();
	}

	## the sequence
	if (index($opt_s, ',') != -1) {
	 print "Remove ',' from $opt_s (-s) when using with split (-p)\n";
	 print STDERR $usage;
	 exit();
	}

	my ($column, $category) = split /,/, $opt_p;
	#print "$column\n", "$category\n";

	## if = ALL, read through the column of the file and get the different groups and build cmmd for each.
	if ($category =~ m/ALL/) {
    print "Creating one table per unique category (group) in the column $column\n";
		my @categories = getUniqCategs($map_file, $column);
		foreach my $categ (@categories){
			if(defined($opt_t)){
				q2table_qza_qzv($opt_t, $map_file, $column, $categ, $out_dir, $opt_s);
			}
		}
	} else{
			if(defined($opt_t)){
				q2table_qza_qzv($opt_t, $map_file, $column, $category, $out_dir, $opt_s);
			}
	}
}



# return the command to make the sub table (and seqs) qza and qzv
sub q2table_qza_qzv {
	my $opt_t = shift;
	my $map_file = shift;
	my $column = shift;
	my $category = shift;
	my $out_dir = shift;
	my $opt_s = shift ; #// "NA"; # use default if not defined

	# replace non alpha-numeric char with underscore
	my $sanitizedCategory = $category;
	$sanitizedCategory =~ s/[^A-Za-z0-9]/_/g;

	my $final_cmd = "qiime feature-table filter-samples --i-table $opt_t --m-metadata-file $map_file --p-where \"$column\='$category'\" --o-filtered-table $out_dir$sanitizedCategory.table.qza";
	print $final_cmd, "\n";

	my $makeqzv = "qiime feature-table summarize --verbose --i-table $out_dir$sanitizedCategory.table.qza --o-visualization $out_dir$sanitizedCategory.table.qzv --m-sample-metadata-file $map_file";
	print $makeqzv, "\n\n";

	push @BIG_CMD_ARRAY, $final_cmd;
	push @BIG_CMD_ARRAY, $makeqzv;

	if(defined($opt_s)){
		q2seqs_w_table_qza_qzv($opt_s, "$out_dir$sanitizedCategory");
	}

}


# requires table
# return the command to make the sub seqs qza and qzv
sub q2seqs_w_table_qza_qzv {
	my $opt_s = shift;
	my $name = shift;

	my $final_cmd = "qiime feature-table filter-seqs --i-data $opt_s --i-table $name.table.qza --o-filtered-data $name.rep-seqs.qza";
	print $final_cmd, "\n";

	my $makeqzv = "qiime feature-table tabulate-seqs --verbose --i-data $name.rep-seqs.qza --o-visualization $name.rep-seqs.qzv";
	print $makeqzv, "\n\n";

	push @BIG_CMD_ARRAY, $final_cmd;
	push @BIG_CMD_ARRAY, $makeqzv;

}



##########################
# build shell command file

# load modules
my $commands_string = "#!/bin/bash";
$commands_string .= "\nexport TMPDIR=/lscratch/\$SLURM_JOB_ID\n\nmodule unload qiime
					module load qiime/2-2019.4\n\n";


# qiime2 commands
$commands_string .= join("\n", @BIG_CMD_ARRAY);

$commands_string=~s/[\t]+//g;

##########################
# run command shell file

open(CMD, ">$commands") || die "Could not open $commands ($!)\n";
print CMD "$commands_string";
close (CMD);
chmod 0755, $commands;


##########################
# set up sbatch job

my $cmd_str = "sbatch --cpus-per-task=2 --time=2:00:00 --mail-type=END $commands";
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
# reads the specific column from the metadata file and returns the unique categories
sub getUniqCategs{
	my $column_separator = "\t";
	my $column_number = 0;
	my $flag = 0;
	my $file = shift;
	my $column_name = shift;
	my %hash;
	open(FILE,"<","$file");
	while (<FILE>) {
		my @columns = split(/$column_separator/);
		if($flag == 0){
			$column_number = first { $columns[$_] eq $column_name } 0..$#columns;
			$flag = 1;
		} else{
			my $key = $columns[$column_number] ;
			$hash{$key} = 1;
			#print $key, "\n";
		}
	}
	close FILE;
	return (keys %hash);
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
