'''
In this problem, we will consider floating point numbers represented by two (unsigned) integers, sig for the significand and exp for the exponent. sig is the full significand (the leading 1 is included). Each significand has 16 bits, and is a bitwise represention of a number between 1 and 2. There is no sign bit.

You are given two numbers in this format: (sig1, exp1) and (sig2, exp2). Compute the representation (sigsum, expsum) of the sum of the two. To ensure you're covering all corner cases, we require that you test your code on multiple numbers. The loop code to do so is given, your code should go only inside the loop.

Hint: To help check your answer, it might be helpful to convert the values (and your sum) into Python's built-in floating point values.

INPUT (within the loop):

sig1, a 16-bit integer representing the significand of the first number.
exp1, an integer representing the exponent of the first number.
sig2, a 16-bit integer representing the significand of the second number.
exp2, an integer representing the exponent of the second number.
OUTPUT (within the loop):

sigsum, a 16-bit integer representing the significand of the sum of the two numbers.
expsum, an integer representing the exponent of the sum of the two numbers.
'''
sig1 = 61921 
sig2 = 59808 
exp1 = 3
exp2 = 0 
value1 = sig1/2**15 * 2**exp1
value2 = sig2/2**15 * 2**exp2
print(value1,value2, value1 + value2)
diffexp = abs(exp1 - exp2)
print(bin(sig1)[2:],bin(sig2)[2:])
if exp2 >= exp1:
    finalexp = exp2
    sig1 = sig1 >> diffexp
else:
    finalexp = exp1
    sig2 = sig2 >> diffexp
sigsum = sig1 + sig2
shift = sigsum.bit_length() - 16
"""
This follwing part is to deal with the overflow case
"""
if sigsum.bit_length() > 16:
    sigsum = sigsum >> shift
    finalexp += shift 

print(sig1,sig2,sigsum)
print(bin(sig1)[2:],bin(sig2)[2:],'\n')
print(bin(sigsum)[2:])
expsum = finalexp
print(expsum)
print(sigsum/2**15 * 2**expsum)
