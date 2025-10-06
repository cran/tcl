# tcl 
## Testing in Conditional Likelihood Context

An implementation of hypothesis testing in an extended Rasch modelling framework including sample
size planning procedures and power computations. Provides 4 statistical tests, i.e. gradient test
(GR),likelihood ratio test (LR), Rao score or Lagrange multiplier test (RS), and Wald test, for
testing a number of hypotheses referring to the Rasch model (RM), linear logistic test model (LLTM),
rating scale model (RSM), and partial credit model (PCM).Three types of functions for power and sample
size computations are provided. Firstly, functions to compute the sample size given a user-speciﬁed
(predetermined) deviation from the hypothesis to be tested, the level alpha, and the power of the
test. Secondly, functions to evaluate the power of the tests given a user-speciﬁed (predetermined)
deviation from the hypothesis to be tested and the level alpha of the test. Thirdly, functions to
evaluate the so-called post hoc power of the tests. This is the power of the tests given observed
deviation of the data from the hypothesis to be tested and a user-specified level alpha of the test.
Power and sample size computations are based on a Monte Carlo simulation approach. It is
computationally very efficient. The variance of the random error in computing power and sample size
arising from the simulation approach is analytically derived by the delta method.
