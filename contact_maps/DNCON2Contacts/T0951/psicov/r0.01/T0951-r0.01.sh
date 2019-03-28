#!/bin/bash
touch T0951-r0.01.running
echo "running T0951-r0.01 .."
date
/storage/htc/bdm/Collaboration/jh7x3/Contact_prediction_with_Tianqi/CAMEO_test/psicov/psicov_lotus/PSICOV/psicov -o -r 0.01 T0951.aln > T0951-r0.01.psicov
if [ -s "T0951-r0.01.psicov" ]; then
   mv T0951-r0.01.running T0951-r0.01.done
   echo "T0951-r0.01 job done."
   date
   exit
fi
mv T0951-r0.01.running T0951-r0.01.failed
echo "psicov job T0951-r0.01 failed!"
date
