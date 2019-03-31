#!/bin/bash
touch T0951-d0.03.running
echo "running T0951-d0.03 .."
date
/storage/htc/bdm/Collaboration/jh7x3/Contact_prediction_with_Tianqi/CAMEO_test/psicov/psicov_lotus/PSICOV/psicov -o -d 0.03 T0951.aln > T0951-d0.03.psicov
if [ -s "T0951-d0.03.psicov" ]; then
   mv T0951-d0.03.running T0951-d0.03.done
   echo "T0951-d0.03 job done."
   date
   exit
fi
mv T0951-d0.03.running T0951-d0.03.failed
echo "psicov job T0951-d0.03 failed!"
date
