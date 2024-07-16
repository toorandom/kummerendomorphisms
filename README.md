# kummerendomorphisms
MAGMA Code to explicitly compute an endomorphism of the Jacobian of a curve y^2=x^5 + h in its associated Kummer Surface and an example on how to use one of this morphism to do primality testing in python.  This code can be modified (see the comments) to calculate the Kummer equivalent of any endomorphism of J.


Config:
MODIFY primetestkummerexample.py variables m and h to decide which m to use when testing primality of 4*m^2*5^n - 1 and which curve to use y^2 = x^5 + h

Run:

$ ./primetestkummerexample.py
