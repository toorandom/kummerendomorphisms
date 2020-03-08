# Calculates primes of the form 4*m^2*5^n -1.
# This is implemented for h \in { 2,3,31,10 } and m \in {1,3,7,11}
# using a Kummer surface endomorphism coming from 
# the Jacobian of y^2 = x^5 + h
# to change m modify that variable and the start_vector variable
# as described in the comment above it.
# To change the curve via h, change the f1,f2,f3,f4 with the desired h
# and the start_vector to match the the curve h value in its name
#
# Marc Paul Noordman & Eduardo Ruiz Duarte

#Formulas for sqrt5 for various values of h
import kumsqrt5endomorphisms as ksq5
#Starting vectors for each curve h and for each m
import kumstartvectors as ksvectors
from math import gcd, log, pow, floor

# Value for m
m = 3

# Point P = 4*m***Q on the Kummer surface,
# in standard coordinates. Use format start_vector_hX_mY
# where X and Y are the values of h and m respectively.
# it is implemented for h \in { 2,3,31,10 } and m \in {1,3,7,11}
start_vector = ksvectors.start_vector_h10_m3

# The polynomials giving multiplication by sqrt 5
# on the Kummer surface, in standard coordinates.
# Use the same h as in the starting vector.

f1 = ksq5.f1_h10
f2 = ksq5.f2_h10
f3 = ksq5.f3_h10
f4 = ksq5.f4_h10

mult_by_sqrt_5 = [f1,f2,f3,f4]

def lambda_mn(m,n):
  return 4*m**2 * 5**n - 1

def evaluate_mod(f, pt, N):
    '''Evaluate a list of functions f at a point pt and
    return the result modulo N. Assumes the results of 
    applying elements of f at pt are integers'''
    return [fi(*pt) % N for fi in f]
    
def test_primality(n, m=m,start_vector = start_vector, mult_by_sqrt_5 = mult_by_sqrt_5):     
    la = lambda_mn(m,n)
    curr = [x % la for x in start_vector]
    prev = [0,0,0,1]
    found_zero = False
    for r in range(0, 2*n+1):
        prev = curr
        curr = evaluate_mod(mult_by_sqrt_5, curr, la)
        if curr[0] == 0 and curr[1] == 0 and curr[2] == 0:
            found_zero = True
            break
    if not found_zero:
        return "Not prime"
    bound = 4*log(pow(la, 1/4) + 1)/log(5) 
    if r > bound:
        possible_divisors = [gcd(x, la) for x in prev]
        for d in possible_divisors:
            if d > 1 and d < la:
                return "Not prime, found divisor {}".format(d)
        return "Prime"
    return "Indeterminate, finished after {} steps (needed at least {} steps)".format(r,floor(bound)+1)


print ("n | Result")
print ("--|---------------")
for n in range(1,500):
    if n%2 == 1:
        print (n,'|', test_primality(n))
