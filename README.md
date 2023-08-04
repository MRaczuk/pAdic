# padic
Library implementing p-adic numbers for python.

This is a package implementing p-adic numbers. Includes Padic class for p-adic numbers representation as well as method of finding roots of nupy.polynomial.Polynomial roots in Z_p via Hensel Lemma and some common functions: log, exp, sin, cos, binomial.

Each p-adic number is represented as an p-adic interval, namely a*p^v + O(p^N). This allows for quick computations with certainty of corectness of computed digits (which is not the case for standard floating point implementation of real numbers as computation errors may produce incorrect results of arithmetic operations).

Please note that optional argument N of log, exp, sin, cos, binomial refers to the number of terms of series defining respective functions NOT the number of calculated correct digits of the final result. This may change in the (near?) future.

Also note that series function gives easy way of defining one's own p-adic valued functions via power series.

Will add more info here one day.
