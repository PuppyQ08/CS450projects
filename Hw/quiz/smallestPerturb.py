'''
Given a positive double precision float x, find the smallest (to within a factor of 2) positive y such that fl(x+y)>x.

INPUT:

x: positive number of type float .
OUTPUT:

y: minimum non-neglibile perturbation to x
'''
x = 123.434535545904894
idx = True
test = 0
i = 0.0
y = 1.0
while i < 1024:
    test = y*2**-1
    if (x + test > x) == False:
        break
    y *= 2**-1
    i+=1
    print(i)
print(i)
