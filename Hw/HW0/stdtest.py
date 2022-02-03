import math
import numpy as np
def std(x):
    tp = 0
    op = 0
    xavg = sum(x,0.0)/float(len(x))
    for xi in x:
        tp += (xi - xavg)**2
        op += xi**2
    tp = math.sqrt(tp/float(len(x)-1))
    op = math.sqrt((op - len(x)*(xavg**2)))/float(len(x)-1)
    return op, tp
sequence = np.array([1,-353345,565,76586547657,464,-65754675856765756765765765765435433654654765,65765,655,4,5])
var_seq_op = 0 
var_seq_tp = 0
var_seq_op, var_seq_tp = std(sequence)
