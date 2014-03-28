=============================================
Bayesian Additive Regression Tree Project
README
Opened: 2013-01-07
=============================================

>> COMPILE AND LINK DIRECTIONS:

Requires 64-bit R to be installed on your system.

1. scp git-afsis@themathpath.com:~/bart.data/AfSIS-Core\ MIR\ first\ derivative.csv bart/data 
2. cd bart/src
3. on command line, type:

[BayesTree/src]$ R CMD COMPILE mbart.cpp BirthDeath.cpp ChangeRule.cpp Funs.cpp Lib.cpp Likelihood.cpp List.cpp MuS.cpp Node.cpp Prior.cpp Rlob.cpp Swap.cpp spLM.cpp util.cpp

4. on command line, type:

[BayesTree/src]$ R CMD SHLIB mbart.cpp BirthDeath.cpp ChangeRule.cpp Funs.cpp Lib.cpp Likelihood.cpp List.cpp MuS.cpp Node.cpp Prior.cpp Rlob.cpp Swap.cpp spLM.cpp util.cpp

This will build mbart.so in the bart/src directory. 
You need to repeat steps 2 through 4 if you modify any .cpp files and want to rebuild the programs.


>> EXECUTION DIRECTIONS:

1. cd bart/src
2. type R on command line, to enter R:
3. at R command prompt, run:
	
source("../R/bartcode_v2_101512.R")

you should then see numbers flying across the screen, and the following files written to MCMCresults:

MCMC100.txt
MCMC200.txt
MCMC300.txt
MCMC400.txt
MCMC500.txt
MCMC600.txt
MCMC700.txt
MCMC800.txt
MCMC900.txt
MCMC1000.txt
MCMC1100.txt

Each of these files currently contains the trees generated at different stages of the MCMC algorithm.


>> TROUBLESHOOTING

If upon executing the bartcode program in R, you see "Killed", it is possible that your computer does not have enough memory to run this program. 
Unfortunately I do not know the exact memory requirements for this program at the moment.
From experimentation, we have two data points: 256 MB is definitely not enough, but 8 GB is definitely enough.
A few gigabytes of memory is most likely preferred. 


2013-01-08

>> BART = Bayesian Additive Regression Tree

> Basic idea of BART

Suppose

	y = f(x) + e

where y is a scalar, x is a vector, e is noise, and f is an unknown function. Our goal is to estimate f, given training data (x_1,y_1),...,(x_N,y_N).

The BART approach to this problem assumes that 

[EQUATION (1)]	f(x) = g_1(x) + g_2(x) + g_3(x) + ... + g_m(x) 

where each g_i is a tree. That is, 

	g_i(x) = g(x,T_i,M_i), 
	
where T_i = ith tree, and M_i = leaf node parameters of ith tree ("mu values").

We treat all these trees as parameters, and use an iterative algorithm known as Markov Chain Monte Carlo (MCMC) to estimate the parameters. 

The code here (e.g. mbart.cpp) achieves all this. We supply it with the training data and desired number of trees "m", and it provides us with the T_i's and M_i's.

As the MCMC algorithm works, the trees are gradually refined. 
For now, our goals are to:

1. Save these forests of trees from different stages of the MCMC algorithm, and
2. Substitute these saved trees into [EQUATION (1)] to generate predictions for arbitrary values of x.



> References for further reading if desired:

http://wumath.wustl.edu/files/math/BART_Loeb_Lecture_St_Louis_11-3-11.pdf
http://math.acadiau.ca/chipmanh/talks/dal-bart.pdf
http://www-stat.wharton.upenn.edu/~edgeorge/Research_papers/BART%20June%2008.pdf
http://rss.acs.unt.edu/Rdoc/library/BayesTree/html/bart.html
http://www.stat.cmu.edu/~cshalizi/350-2006/lecture-10.pdf



> Summary of how program works:

1. The code generates a collection of trees -- a forest.
2. Each tree corresponds to a decision tree that operates on a vector.

Tree1
Depth:0node:1 n:361 TBN: 100 Var:832 ORDRule:(36)=-0.597998
Depth:1 node:1 n:107 TBN: 001 Var:625 ORDRule:(27)=-0.823074
Depth:2  node:1 n:22 TBN: 010
Depth:2  node:1 n:85 TBN: 010
Depth:1 node:1 n:254 TBN: 000 Var:733 ORDRule:(99)=0.944291
Depth:2  node:1 n:232 TBN: 010
Depth:2  node:1 n:22 TBN: 001 Var:1745 ORDRule:(44)=-0.385226
Depth:3   node:1 n:2 TBN: 010
Depth:3   node:1 n:20 TBN: 010

This Tree1 says:
	Compare 832nd entry of vector with -0.597998. 
	If < -0.597998, then take left branch; otherwise, take right branch.
	If you took left branch, now you are comparing 625nd entry of vector with -0.823704, and so on.
	You keep going down this decision tree until you hit a leaf (which has no ORDRule). All the leaves have a "mu" value associated with them, stored somewhere. The mu value of the leaf you end up in is the final output of this tree.

Then there is another tree:

Tree2
Depth:0node:1 n:361 TBN: 100 Var:1701 ORDRule:(92)=0.777852
Depth:1 node:1 n:275 TBN: 010
Depth:1 node:1 n:86 TBN: 001 Var:865 ORDRule:(21)=-1.030811
Depth:2  node:1 n:44 TBN: 010
Depth:2  node:1 n:42 TBN: 010

You do the same traversal on this tree, and end up with a mu value as well.

Summing up all the mu values gives you the final answer.



 William Wu, 2013-01-08
