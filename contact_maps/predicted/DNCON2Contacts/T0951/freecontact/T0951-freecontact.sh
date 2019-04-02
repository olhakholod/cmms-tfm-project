#!/bin/bash
touch freecontact.running
echo "running freecontact .."
/storage/htc/bdm/tools/freecontact-1.0.21/bin/freecontact < T0951.aln > T0951.freecontact.rr
if [ -s "T0951.freecontact.rr" ]; then
   mv freecontact.running freecontact.done
   echo "freecontact job done."
   exit
fi
echo "freecontact failed!"
mv freecontact.running freecontact.failed
