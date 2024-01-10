Author: Brendan Gongol

###########################################################################################
Goal 

Explore the notion of chance and how it applies to significance measures. Experiment with simulation studies using theoretical 
normal distribution. Understand the role of assumptions made by specific tests (i.e. t-test).

###########################################################################################

Questions

1. Using simulation based on rnorm, study the performance of “unbiased” estimator of standard deviation (using n-1 in denominator) 
and two nearest “alternatives” that use “n” or “n-2” in the denominator instead:

- Go over samples ranging in size from 2 to 20;
- For each sample size n, draw multiple (say, N=1000) samples of that size from normal distribution with zero mean and standard deviation of one;
- for each such sample calculate (a) standard deviation within the sample; (b) biased sd alternative that uses n instead of n−1 in the denominator;
(c) sd alternative that uses n−2 in denominator;
- at this point, for each sample you have three estimators calculated, let’s call them s, s.n, s.n.1. We want to know 
how well each of them approximates the true underlying standard deviation. In order to accomplish this, calculate the mean and spread of each of these
three standard deviation estimates across all N samples for the current sample size n.
- Plot the calculated means and spreads of all three estimators of the standard deviation as functions of sample size n, 
observe how they approach standard deviation of the underlying distribution.

2. Study the behavior of t-test under the null hypothesis using numerical simulation:

- For a normal distribution with arbitrarily chosen mean and variance (e.g., 2 and 3) draw two samples of different sizes (e.g. 5 and 7),
calculate the difference between samples’ means and also the p-value of t-test run on the same two samples.
- Then repeat this large number of times (e.g. N=104). From the differences between means obtained in those resamplings, 
calculate the brute-force p-value of each of those N differences, using the rank of that particular difference among all N
differences obtained in the simulation.

- Plot:
- the distribution of differences observed in all the resamplings 
- the distribution of t-test p-values obtained for the pairs of samples in each resampling.
- the distribution of brute force p-values obtained for the pairs of samples in each resampling. 
- generate a scatter plot of t-test p-values versus the differences between the means in each pair of samples drawn. 
- a scatter plot of brute force p-values versus the differences between the means (in a scatterplot, each point should represent 
a single resampling experiment, x-axis=difference, y-axis=p-value). 
- a scatter plot of t-test p-values vs brute force p-values. 

Repeat this study for 4 times larger sample sizes

3. Using rnorm create two samples that result in insignificant t-test p-value, but represent two very different processes:

- Plot the two distributions you have come up with. 
- Repeat the following multiple times: draw one sample from each of the two distributions you chose, calculate t-test p-value
between such two samples. Plot the resulting distribution of p-values obtained across multiple trials, consider what it tells us and 
how consistent it is with what exactly t-test is testing for.

