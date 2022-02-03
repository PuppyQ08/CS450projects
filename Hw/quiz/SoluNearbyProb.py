import numpy as np
import math 
a = np.linspace(0,0.5,num=10000,endpoint=True)
xhat_approx = 0
leastdiff = 10
for ele in a:
    fxh = 1 + ele + 0.5*ele**2
    diff = abs(fxh - math.e**0.25)
    if diff < leastdiff:
        xhat_approx = ele
        leastdiff = diff

print(leastdiff, xhat_approx)
