# This is exactly primetestkummerexample.py but made multiprocessor
# note that this code was made by an LLM by prompting to add multiprocessing to the primetestkummerexample
# but the results should be the same.
# Eduardo Ruiz Duarte

import concurrent.futures
import os
import sys
from math import gcd, log, pow, floor

# Import the required modules (assuming kumsqrt5endomorphisms and kumstartvectors are available)
import kumsqrt5endomorphisms as ksq5
import kumstartvectors as ksvectors

print (sys.argv)
# Value for m (4*m^2 * 5^n - 1) and h for the hyperelliptic curve y^2= x^5 + h to be used.
m = 7
h = 10

# Coordinates for the starting point P = 4*m***Q on the Kummer surface
start_vector = getattr(ksvectors, 'start_vector_h' + str(h) + "_m" + str(m))

# The polynomials giving multiplication by sqrt 5 on the Kummer surface
f1 = getattr(ksq5, 'f1_h' + str(h))
f2 = getattr(ksq5, 'f2_h' + str(h))
f3 = getattr(ksq5, 'f3_h' + str(h))
f4 = getattr(ksq5, 'f4_h' + str(h))

mult_by_sqrt_5 = [f1, f2, f3, f4]

def lambda_mn(m, n):
    return 4 * m**2 * 5**n - 1

def evaluate_mod(f, pt, N):
    '''Evaluate a list of functions f at a point pt and return the result modulo N.
    Assumes the results of applying elements of f at pt are integers'''
    return [fi(*pt) % N for fi in f]

def test_primality(n, m=m, start_vector=start_vector, mult_by_sqrt_5=mult_by_sqrt_5): 
    process_id = os.getpid()  # Get the process ID   
    la = lambda_mn(m, n)
    curr = [x % la for x in start_vector]
    prev = [0, 0, 0, 1]
    found_zero = False
    for r in range(0, 2 * n + 1):
        prev = curr
        curr = evaluate_mod(mult_by_sqrt_5, curr, la)
        if curr[0] == 0 and curr[1] == 0 and curr[2] == 0:
            found_zero = True
            break
    if not found_zero:
        return f"4*{m}^2*5^{n}-1 | Not prime | Process ID: {process_id}"
    bound = 4 * log(pow(la, 1/4) + 1) / log(5)
    if r > bound:
        possible_divisors = [gcd(x, la) for x in prev]
        for d in possible_divisors:
            if d > 1 and d < la:
                return f"4*{m}^2*5^{n}-1 | Not prime, found divisor {d} | Process ID: {process_id}"
        return f"4*{m}^2*5^{n}-1 | Prime | Process ID: {process_id}"
    return f"4*{m}^2*5^{n}-1 | Indeterminate, finished after {r} steps (needed at least {floor(bound) + 1} steps) | Process ID: {process_id}"

def process_primality_check(n):
    if n % 2 == 1:
        return test_primality(n)
    return None

def main():
    print("n | Result | Process ID")
    print("--|---------------------")
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_primality_check, n) for n in range(1, 1000)]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result is not None:
                print(result)

if __name__ == '__main__':
    main()
