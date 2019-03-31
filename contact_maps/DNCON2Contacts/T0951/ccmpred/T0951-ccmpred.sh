#!/bin/bash
touch ccmpred.running
echo "running ccmpred .."
/storage/htc/bdm/Collaboration/jh7x3/Contact_prediction_with_Tianqi/CAMEO_test/psicov/psicov_lotus/CCMpred-0.1.0/bin/ccmpred -t 8 T0951.aln T0951.ccmpred > ccmpred.log
if [ -s "T0951.ccmpred" ]; then
   mv ccmpred.running ccmpred.done
   echo "ccmpred job done."
   exit
fi
echo "ccmpred failed!"
mv ccmpred.running ccmpred.failed
