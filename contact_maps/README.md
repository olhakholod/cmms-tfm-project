# Contact Maps
------
This directory ('cmms-tfm-project/contact\_maps') contains residue-residue contact maps stored in the 
[CASP RR format](http://predictioncenter.org/casprol/index.cgi?page=format#RR) and 
split into two subdirectories:

1. The [predicted](./predicted) subdirectory contains '.rr' files that were predicted directly
from the protein sequences, and are organized by the tool used for the prediction of a target. Each tool's subdirectory includes results for both of our targets. The corresponding references for each tool are included here:
  * [CMApro](https://www.ics.uci.edu/~baldig/index.html).
  * [DNCON2](https://github.com/multicom-toolbox/DNCON2).
  * [SVMcon](http://sysbio.rnet.missouri.edu/multicom_toolbox/tools.html) and 
  also [here](http://scratch.proteomics.ics.uci.edu/explanation.html), the reference paper can be found 
  [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1852326/)
  * [RaptorX-Contact](http://raptorx.uchicago.edu/ContactMap/), reference 
  paper [here](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005324)
  * [NeBcon](https://zhanglab.ccmb.med.umich.edu/NeBcon/)

2. The [true](./true) subdirectory contains .rr files that were derived from
the .pdb files of the native structures for our targets ([5uw2](../pdb/native/T0866/5uw2.pdb) 
and [5z82](../pdb/native/T0951_easy/5z82.pdb) for CASPs T0866 and T0951). The .rr files are also
organized by the tools used to determine the .rr files, starting from their respective PDB file:
  * [ConEva](http://iris.rnet.missouri.edu/coneva/)
  * [Our Python3 implementation](../scripts/3d_to_contact.py)
