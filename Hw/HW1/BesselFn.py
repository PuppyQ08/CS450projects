'''
The Bessel functions Jn(⋅) obey the recurrence relation
Jn+1(z)=(2n/z)Jn(z)−Jn−1(z).
This means that, given J0(z) and J1(z) for a fixed z∈ℝ∖{0}, Jn(z) for integers n≥2 can be computed--at least in theory.

Side note: If you're curious what Bessel functions are good for "in real life", they describe the behavior of waves in round geometries, along the 'radial' direction. See this video clip for an example of what that might look like.

You can evaluate Bessel functions with Jn(z)=𝚜𝚌𝚒𝚙𝚢.𝚜𝚙𝚎𝚌𝚒𝚊𝚕.𝚓𝚟(𝚗,𝚣).

Use the values for J0(z) and J1(z) (obtained from scipy) and the recurrence relation above to compute the values of J2(z),…,J49(z). The Bessel function values J2(z),⋯,J49(z) should be computed using the recurrence relation only.

Compute and store these values in a numpy array called jn. This array should contain the values J0(z),⋯,J49(z).

Compute the relative error compared to the value returned by the 𝚜𝚌𝚒𝚙𝚢 function you used in part 1. Store these values in a numpy array called e1_rel. You should observe that the recurrence leads to very inaccurate results start around n=30.

Now, compute the absolute error incurred at each single step of the recurrence. Use the scipy values of Jn(z) and Jn−1(z) to compute Jn+1(z) for 1≤n≤48. Compute the absolute error with respect to the Jn(z) given by scipy. Store these values for the error of J0(z),⋯,J49(z) into a numpy array called e2_abs (the first two elements should be zero).

Let en be your absolute error in computing Jn (what you store in e2_abs). Create a model for the relative error in part (1), by considering how the absolute error of each step of the recurrence becomes amplified throughout the recurrence. In particular, we can model the overall absolute error of Jn+1 calculated from the recurrence started from J0 and J1 as

En+1=|en+1+(2n/z)En+En−1|.
Calculate En for each n∈{0,…,49}, then use it to compute the relative error for Jn and store the result in e3_rel (the first two elements should be zero). Plot e1_rel and e3_rel, use a logarithmic scale for the y-axis. Format and label your plots appropriately.

INPUT:

z: The point at which the Bessel Functions are evaluated.
OUTPUT:

jn: Array of size 50 of Bessel function values up to n=49 computed via recurrence (part 1).
e1_rel: Array of size 50 of relative error Bessel function values computed via recurrence with respect to scipy (part 1).
e2_abs: Array of size 50 of absolute error incurred with at each step of the Bessel function recurrence (part 2).
e3_rel: Array of size 50 of model relative error for complete Bessel function recurrence (part 3).
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv
import math
def myjv(n0,n1,n,z,tgt,jn):
    if n == tgt -1 :
        return
    jn.append((2*n/z)*n1 - n0)
    myjv(jn[-2],jn[-1],n+1,z,tgt,jn)
    return
def myjv2(n0,n1,n,z,tgt,jn2):
    if n == tgt -1 :
        return
    jn2.append((2*n/z)*n1 - n0)
    myjv2(jv(n,z),jv(n+1,z),n+1,z,tgt,jn2)
    return
def err(n0,n1,e2_abs,n,z,tgt,en):
    if n == tgt -1 :
        return
    en.append(e2_abs[n+1] + (2*n/z)*n1 + n0)
    err(en[-2],en[-1],e2_abs,n+1,z,tgt,en)
    return



n = 50 
z = 1
jn = [jv(0,z),jv(1,z)]
jn_corr = jv(range(0,50),z)
myjv(jv(0,z),jv(1,z),1,z,50,jn)
jn = np.array(jn)
e1_rel = abs(jn-jn_corr)/jn_corr


jn2 = [jv(0,z),jv(1,z)]
myjv2(jv(0,z),jv(1,z),1,z,50,jn2)
jn2 = np.array(jn2)
e2_abs = abs(jn2-jn_corr)

en = [0,0]
err(0,0,e2_abs,1,z,50,en)
e3_rel = np.array(en)/jn_corr
print(e3_rel)
e3 = []
e1 = []
for i in range(3,50):
    e3.append(math.log10(e3_rel[i]))
    e1.append(math.log10(e1_rel[i]))
for j in range(47):
    print(e3[j]-e1[j])
plt.plot(list(range(47)),e1,label = "e1")
plt.plot(list(range(47)),e3,label = "e3")
plt.show()
