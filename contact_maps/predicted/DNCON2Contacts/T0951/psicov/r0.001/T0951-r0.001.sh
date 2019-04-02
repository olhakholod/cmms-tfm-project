#!/bin/bash
touch T0951-r0.001.running
echo "running T0951-r0.001 .."
date
/storage/htc/bdm/Collaboration/jh7x3/Contact_prediction_with_Tianqi/CAMEO_test/psicov/psicov_lotus/PSICOV/psicov -o -r 0.001 T0951.aln > T0951-r0.001.psicov
if [ -s "T0951-r0.001.psicov" ]; then
   mv T0951-r0.001.running T0951-r0.001.done
   echo "T0951-r0.001 job done."
   date
   exit
fi
mv T0951-r0.001.running T0951-r0.001.failed
echo "psicov job T0951-r0.001 failed!"
date
