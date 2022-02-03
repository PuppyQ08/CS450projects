import math
import numpy as np
def std(x):
    op = 0
    tp = 0
    xavg = sum(x,0.0)/float(len(x))
    for xi in x:
        tp += (xi - xavg)**2
        op += xi**2
    tp = math.sqrt(tp/float(len(x)-1))
    op = math.sqrt((op - len(x)*(xavg**2))/float(len(x)-1))
    return op, tp
#sequence = np.array([1,3,55,74,6575467585676575.6765754654765,6,65,4,5,34,46,6,56,5,76,54,4,5,4,4])
sequence = np.linspace(0.0,1.0,num = 10000)
sequence += 10E4
var_seq_tp = 0
var_seq_op = 0
var_seq_op, var_seq_tp = std(sequence)
error = (var_seq_op - var_seq_tp)/var_seq_op
print(var_seq_tp,var_seq_op)
print(error)

