2012-12-22

list of files changed by JC:

[~/Documents/programming/c/BayesTree]$ git diff --stat
 R/bart.R           |    2 +-
 src/BirthDeath.cpp |    2 +-
 src/CPriParams.h   |    1 -
 src/Node.cpp       |   12 ++++++------
 src/Node.h         |    1 -
 src/mbart.cpp      |   33 ++++++++++++++++++++++++++-------
 6 files changed, 34 insertions(+), 17 deletions(-)



2012-12-23

Note 1: BirthDeath.so is apparently not being generated when using the compile and run commands described in linkRwithC++.txt. Why?
2012-12-22

list of files changed by JC:

[~/Documents/programming/c/BayesTree]$ git diff --stat
 R/bart.R           |    2 +-
 src/BirthDeath.cpp |    2 +-
 src/CPriParams.h   |    1 -
 src/Node.cpp       |   12 ++++++------
 src/Node.h         |    1 -
 src/mbart.cpp      |   33 ++++++++++++++++++++++++++-------
 6 files changed, 34 insertions(+), 17 deletions(-)



2012-12-23

Note 1: BirthDeath.so is apparently not being generated when using the compile and run commands described in linkRwithC++.txt. Why?

Note 2: To run in R, you must use R64 and also have BirthDeath.so in the current working directory.


2013-01-08

Based on today's experiments it seems that Note 2 above may be wrong, which renders Note 1 above irrelevant.