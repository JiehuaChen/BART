COMPILE, LINK AND RUN DIRECTIONS:
1. Compile, link and run mbart.
2. Go to "src_prediction/" directory.
3. On command line, type:
[bart/src_prediction]$ make
[bart/src_prediction]$ ./predict


Explaination for predict program:
1. Input files: 
   The test data file: "MCMCresults/xdat.txt".
   The file of one trained forest: "MCMCresults/MCMC1100.txt".
2. Output file:
   The file for storing the predicted results of "xdat.txt": "predictedY.txt". The predicted Y have some errors compared with those in file "MCMCresults/Y.txt"
