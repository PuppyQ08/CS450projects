'''
Safeguarded Methods
Newton's method while quadratically convergent may fail to converge if the starting guess is not close enough to the solution. On the other hand, the bisection method is guarranteed to converge but it does so rather slowly. In this problem you'll code what is know as a safeguarded method that combines the fast convergence of Newton's method with the safety of the bisection method.



Run one step of Newton's method from the previous approximate solution or from the middle of the current bracket

Check that the new Newton's method iterate remains inside the current bracket

If it doesn't, run one step of bisection, reducing the size of the bracket by 2, and taking the approximate solution to be the middle of the new interval
If it does, accept the Newton's iterate as the new approximate solution
You are given a function and its derivative as well as an initial bracket. Code the safeguarded method as explained (starting guess should be middle of bracket or result of previous Newton's method iteration) and save the approximate solution (including first guess) at every iteration. Stop when the function value is within the given tolerance from zero.
'''
import numpy as np

x_star = 99.0

def f(x):
    global x_star

    return 2.0 / (1.0 + np.exp(-(x - x_star))) - 1.0

def df(x):
    global x_star

    return (2.0 * np.exp(x + x_star)) / (np.exp(x) + np.exp(x_star)) ** 2.0

left, right = np.array([-100.0, 100.0])
tol = 1.0e-8


x = [(right+left)/2]

while np.abs(f(x[-1]))>tol:
    # Newton step
    x.append(x[-1] - f(x[-1]) / df(x[-1]))

    # Implement safeguarding here
    # ...
    
    if x[-1] < left or x[-1] > right:
        print("Warning: exited bracket")
        mid = (right + left)/2
        if np.sign(f(mid)) == np.sign(f(right)):
            right = mid
            x.append(mid)
        elif np.sign(f(mid)) == np.sign(f(left)):
            left = mid
            x.append(mid)
    else:
        pass


x = np.asarray(x)

