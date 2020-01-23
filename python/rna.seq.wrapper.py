import sys
import os
from glob import glob
from utils import *

# list all the fastq gz files in the directory.
# keep track whether you want to do single or paired end.

# usage: python /data/rodriguesrr/scripts/python/rna.seq.wrapper.py ../MM12-july-19/MM12-july-19/MMP12-138637510/FASTQ_Generation_2019-07-19_07_43_04Z-190032047 --mouse --ensembl --version GRCm38 --outputFolder mmp12 --lexogen --threads 8
# usage: python /data/rodriguesrr/scripts/python/rna.seq.wrapper.py ../MM12-july-19/MM12-july-19/MMP12-138637510/FASTQ_Generation_2019-07-19_07_43_04Z-190032047 --human --ensembl --version GRCh37 --outputFolder mmp12 --lexogen --threads 8

cmd_prefix = 'perl /data/rodriguesrr/scripts/perl/RNA-seq-pipeline-bw.pl'

if '--help' in ' '.join(sys.argv):
    print "This script creates a swarm file where the following command is run for each fastq file.\n"
    os.system(cmd_prefix + ' --help')
    sys.exit(0)

inputfolder = sys.argv[1]
rem_cmd = ' '.join(sys.argv[2:]) #sys.argv[2]

cmd_prefix = cmd_prefix + ' --file'


#print inputfolder
#print rem_cmd


# check for flags:
## lexogen
## paired


files = ''
if '--paired' in rem_cmd:
    sys.exit("Can only process lexogen (single end) reads. Rerun without the --paired arg.")
elif '--lexogen' in rem_cmd:
    print "Lexogen kit used. Searching single end read files."
    files = [y for x in os.walk(inputfolder) for y in glob(os.path.join(x[0], '*_R1_001.fastq.gz'))]
    if len(files)==0:
        sys.exit("Found nothing here. Retry with different arguments!")
    print "Found " + str(len(files)) + " samples"
    #print '\n'.join(files)
else:
    sys.exit("Non-lexogen kit used. Code needs to be fixed so it can handle this case and paired end samples.")


print "Building the command for swarm jobs."
cmds_list = list()
for f in files:
    cmd = cmd_prefix + ' ' + f + ' ' + rem_cmd
    cmds_list.append(cmd)


writeLIST_to_file(cmds_list, 'swarm.job.sh')
os.chmod('swarm.job.sh', 0755)


# build the command to launch the swarm job and tell the user to check the command and then launch it.
swrmcmd = 'swarm -f swarm.job.sh -g 50 -t 8 --module bbtools/38.42,tophat/2.1.1,cutadapt/2.3,samtools/1.9,htseq/0.9.1 --job-name rnaseq --time 10:00:00 --sbatch "--mail-type=END" --logdir swarm.logs'
print "\nCheck and (copy/paste) run the following command:\n"
print swrmcmd, "\n"
