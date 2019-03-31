#!/bin/bash
touch hhb-cov60.running
echo "running hhblits job hhb-cov60.."
/storage/htc/bdm/tools/HHsuite/hhsuite-2.0.16/bin/hhblits -i T0953s2.fasta -d /storage/htc/bdm/tools/databases/uniprot20_2016_02/uniprot20_2016_02 -oa3m hhb-cov60.a3m -cpu 2 -n 3 -maxfilt 500000 -diff inf -e 0.001 -id 99 -cov 60 > hhb-cov60-hhblits.log
if [ ! -f "hhb-cov60.a3m" ]; then
   mv hhb-cov60.running hhb-cov60.failed
   echo "hhblits job hhb-cov60 failed!"
   exit
fi
egrep -v "^>" hhb-cov60.a3m | sed 's/[a-z]//g' > hhb-cov60.aln
if [ -f "hhb-cov60.aln" ]; then
   mv hhb-cov60.running hhb-cov60.done
   echo "hhblits hhb-cov60 job done."
   exit
fi
echo "Something went wrong! hhb-cov60.aln file not present!"
mv hhb-cov60.running hhb-cov60.failed
