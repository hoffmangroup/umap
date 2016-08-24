# Bismap/Umap job submission commands follow.

JOBID=$(qsub -q hoffmangroup -t 1-3114 -N Bismap.UniqueKmers -terse -tc 120 -hold_jid 1 -cwd -b y -o /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UniqueKmers.LOG -e /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UniqueKmers.ERR python get_kmers.py /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/chrsize.tsv /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k24 k24 /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/chrs -var_id SGE_TASK_ID)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -t 1-3114 -N Bismap.RunBowtie -terse -tc 120 -hold_jid 1,$WAITID -cwd -b y -o /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.RunBowtie.LOG -e /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.RunBowtie.ERR python run_bowtie.py /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k24 /mnt/work1/software/bowtie/1.1.0 /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/genome Umap_bowtie.ind -var_id SGE_TASK_ID)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -t 1-24:1 -N Bismap.UnifyBowtie -terse -tc 120 -hold_jid $WAITID -cwd -b y -o /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UnifyBowtie.LOG -e /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UnifyBowtie.ERR python unify_bowtie.py /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k24 /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/chrsize.tsv -var_id SGE_TASK_ID)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -N Moving.Intermediary.Files -hold_jid $WAITID -terse -cwd -b y -o Bismap.FileMov.LOG -e Bismap.FileMov.ERR mv /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k24/*kmer* /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k24/TEMPs)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -N Moving.Intermediary.Files -hold_jid $WAITID -terse -cwd -b y -o Bismap.FileMov.LOG -e Bismap.FileMov.ERR mv /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k24/*bowtie* /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k24/TEMPs)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -t 1-3114 -N Bismap.UniqueKmers -terse -tc 120 -hold_jid $WAITID -cwd -b y -o /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UniqueKmers.LOG -e /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UniqueKmers.ERR python get_kmers.py /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/chrsize.tsv /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k36 k36 /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/chrs -var_id SGE_TASK_ID)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -t 1-3114 -N Bismap.RunBowtie -terse -tc 120 -hold_jid 1,$WAITID -cwd -b y -o /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.RunBowtie.LOG -e /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.RunBowtie.ERR python run_bowtie.py /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k36 /mnt/work1/software/bowtie/1.1.0 /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/genome Umap_bowtie.ind -var_id SGE_TASK_ID)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -t 1-24:1 -N Bismap.UnifyBowtie -terse -tc 120 -hold_jid $WAITID -cwd -b y -o /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UnifyBowtie.LOG -e /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UnifyBowtie.ERR python unify_bowtie.py /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k36 /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/chrsize.tsv -var_id SGE_TASK_ID)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -N Moving.Intermediary.Files -hold_jid $WAITID -terse -cwd -b y -o Bismap.FileMov.LOG -e Bismap.FileMov.ERR mv /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k36/*kmer* /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k36/TEMPs)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -N Moving.Intermediary.Files -hold_jid $WAITID -terse -cwd -b y -o Bismap.FileMov.LOG -e Bismap.FileMov.ERR mv /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k36/*bowtie* /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k36/TEMPs)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -t 1-3114 -N Bismap.UniqueKmers -terse -tc 120 -hold_jid $WAITID -cwd -b y -o /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UniqueKmers.LOG -e /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UniqueKmers.ERR python get_kmers.py /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/chrsize.tsv /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k50 k50 /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/chrs -var_id SGE_TASK_ID)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -t 1-3114 -N Bismap.RunBowtie -terse -tc 120 -hold_jid 1,$WAITID -cwd -b y -o /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.RunBowtie.LOG -e /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.RunBowtie.ERR python run_bowtie.py /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k50 /mnt/work1/software/bowtie/1.1.0 /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/genome Umap_bowtie.ind -var_id SGE_TASK_ID)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -t 1-24:1 -N Bismap.UnifyBowtie -terse -tc 120 -hold_jid $WAITID -cwd -b y -o /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UnifyBowtie.LOG -e /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UnifyBowtie.ERR python unify_bowtie.py /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k50 /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/chrsize.tsv -var_id SGE_TASK_ID)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -N Moving.Intermediary.Files -hold_jid $WAITID -terse -cwd -b y -o Bismap.FileMov.LOG -e Bismap.FileMov.ERR mv /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k50/*kmer* /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k50/TEMPs)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -N Moving.Intermediary.Files -hold_jid $WAITID -terse -cwd -b y -o Bismap.FileMov.LOG -e Bismap.FileMov.ERR mv /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k50/*bowtie* /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k50/TEMPs)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -t 1-3114 -N Bismap.UniqueKmers -terse -tc 120 -hold_jid $WAITID -cwd -b y -o /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UniqueKmers.LOG -e /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UniqueKmers.ERR python get_kmers.py /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/chrsize.tsv /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k76 k76 /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/chrs -var_id SGE_TASK_ID)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -t 1-3114 -N Bismap.RunBowtie -terse -tc 120 -hold_jid 1,$WAITID -cwd -b y -o /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.RunBowtie.LOG -e /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.RunBowtie.ERR python run_bowtie.py /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k76 /mnt/work1/software/bowtie/1.1.0 /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/genome Umap_bowtie.ind -var_id SGE_TASK_ID)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -t 1-24:1 -N Bismap.UnifyBowtie -terse -tc 120 -hold_jid $WAITID -cwd -b y -o /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UnifyBowtie.LOG -e /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UnifyBowtie.ERR python unify_bowtie.py /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k76 /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/chrsize.tsv -var_id SGE_TASK_ID)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -N Moving.Intermediary.Files -hold_jid $WAITID -terse -cwd -b y -o Bismap.FileMov.LOG -e Bismap.FileMov.ERR mv /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k76/*kmer* /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k76/TEMPs)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -N Moving.Intermediary.Files -hold_jid $WAITID -terse -cwd -b y -o Bismap.FileMov.LOG -e Bismap.FileMov.ERR mv /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k76/*bowtie* /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k76/TEMPs)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -t 1-3114 -N Bismap.UniqueKmers -terse -tc 120 -hold_jid $WAITID -cwd -b y -o /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UniqueKmers.LOG -e /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UniqueKmers.ERR python get_kmers.py /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/chrsize.tsv /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k100 k100 /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/chrs -var_id SGE_TASK_ID)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -t 1-3114 -N Bismap.RunBowtie -terse -tc 120 -hold_jid 1,$WAITID -cwd -b y -o /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.RunBowtie.LOG -e /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.RunBowtie.ERR python run_bowtie.py /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k100 /mnt/work1/software/bowtie/1.1.0 /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/genome Umap_bowtie.ind -var_id SGE_TASK_ID)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -t 1-24:1 -N Bismap.UnifyBowtie -terse -tc 120 -hold_jid $WAITID -cwd -b y -o /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UnifyBowtie.LOG -e /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UnifyBowtie.ERR python unify_bowtie.py /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k100 /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/chrsize.tsv -var_id SGE_TASK_ID)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -N Moving.Intermediary.Files -hold_jid $WAITID -terse -cwd -b y -o Bismap.FileMov.LOG -e Bismap.FileMov.ERR mv /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k100/*kmer* /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k100/TEMPs)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -N Moving.Intermediary.Files -hold_jid $WAITID -terse -cwd -b y -o Bismap.FileMov.LOG -e Bismap.FileMov.ERR mv /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k100/*bowtie* /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k100/TEMPs)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -t 1-3114 -N Bismap.UniqueKmers -terse -tc 120 -hold_jid $WAITID -cwd -b y -o /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UniqueKmers.LOG -e /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UniqueKmers.ERR python get_kmers.py /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/chrsize.tsv /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k152 k152 /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/chrs -var_id SGE_TASK_ID)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -t 1-3114 -N Bismap.RunBowtie -terse -tc 120 -hold_jid 1,$WAITID -cwd -b y -o /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.RunBowtie.LOG -e /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.RunBowtie.ERR python run_bowtie.py /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k152 /mnt/work1/software/bowtie/1.1.0 /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/genome Umap_bowtie.ind -var_id SGE_TASK_ID)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -t 1-24:1 -N Bismap.UnifyBowtie -terse -tc 120 -hold_jid $WAITID -cwd -b y -o /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UnifyBowtie.LOG -e /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.UnifyBowtie.ERR python unify_bowtie.py /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k152 /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/chrsize.tsv -var_id SGE_TASK_ID)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -N Moving.Intermediary.Files -hold_jid $WAITID -terse -cwd -b y -o Bismap.FileMov.LOG -e Bismap.FileMov.ERR mv /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k152/*kmer* /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k152/TEMPs)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -N Moving.Intermediary.Files -hold_jid $WAITID -terse -cwd -b y -o Bismap.FileMov.LOG -e Bismap.FileMov.ERR mv /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k152/*bowtie* /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/k152/TEMPs)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"

JOBID=$(qsub -q hoffmangroup -N CombineUmappedFiles -terse -t 1-24 -hold_jid $WAITID -cwd -b y -o /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.combine.LOG -e /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers/Bismap.combine.ERR python combine_umaps.py /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/kmers /mnt/work1/users/hoffmangroup/mkarimzadeh/2016/mappability/data/Umap/hg19/chrsize.tsv -var_id SGE_TASK_ID)
IDPARTS=($(echo $JOBID | tr "." "\n"))
WAITID=${IDPARTS[0]}
echo "Submitted JOB ID $WAITID"
