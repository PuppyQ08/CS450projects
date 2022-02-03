import numpy as np
import numpy.linalg as la
from random import randrange
import matplotlib.pyplot as pt

nnodes = 10
circuit = [{i+1: randrange(1, 50) * 100}
            if i+1 < nnodes
            else {}
           for i in range(nnodes)]

nedges = int(nnodes**2 * 0.2)

for i in range(nedges):
    ijnode = randrange((nnodes-1)*nnodes/2)
    jnode = 0
    ijnode0 = 0
    while True:
        jnode += 1
        ijnode0 += jnode
        if ijnode0 > ijnode:
            break
    inode = ijnode - ijnode0 + jnode
    circuit[inode][jnode] = randrange(1, 50) * 100

fixed_voltages = {0: 5, nnodes-1: 0}

def net_to_dot(net, fixed_voltages, voltages=None):
    lines = []
    for from_node_nr in range(len(net)):
        connections = net[from_node_nr]
        for to_node_nr, item in connections.items():
            label = "%r Ohms" % item
            arrow = "none"
            color = "black"

            lines.append(
            'node%d -> node%d [label="%s",arrowhead="%s", fontcolor="%s"];'
            % (from_node_nr, to_node_nr, label, arrow, color))

    if voltages is not None:
        for node_nr in range(len(net)):
            lines.append(
            'node%d [label="node%d: %.2f V"];'
            % (node_nr, node_nr, voltages[node_nr]))

    for node_nr, voltage in fixed_voltages.items():
        lines.append('node%d -> voltage%d [arrowhead="none"];\n'
                 'voltage%d [label="%s",shape=box];\n'
                 % (node_nr, node_nr, node_nr, "%r V" % voltage))

    return "digraph network {\n%s\n}" % "\n".join(lines)

def plot_circuit(circuit, fixed_voltages, voltages=None):
    from tempfile import TemporaryDirectory
    from os.path import join
    from subprocess import check_call

    with TemporaryDirectory() as td:
        dot_name = "x.dot"
        with open(join(td, dot_name), "w") as dotf:
            dotf.write(net_to_dot(circuit, fixed_voltages, voltages=voltages))

        png_name = "x.png"
        with open(join(td, png_name), "wb") as pngf:
            check_call(["dot", "-Gdpi=200", dot_name, "-Tpng"], stdout=pngf, cwd=td)

        from scipy.misc import imread
        png_data = imread(join(td, png_name))

        pt.figure(figsize=(13,10))
        pt.imshow(png_data)
        pt.setp(pt.gca().get_xticklabels(), visible=False)
        pt.setp(pt.gca().get_yticklabels(), visible=False)
        pt.setp(pt.gca().get_xticklines(), visible=False)
        pt.setp(pt.gca().get_yticklines(), visible=False)


network_mat = np.zeros((nnodes,nnodes))
for i in range(nnodes):
    for node, R in circuit[i].items():
        network_mat[i][i]+= 1/R
        network_mat[i][node] = -1/R
        network_mat[node][i] = -1/R
        network_mat[node][node]+=1/R
rhs = np.zeros(nnodes)
mat = network_mat.copy()
for node, V in fixed_voltages.items():
    rhs[node] = V
    for i in range(nnodes):
        mat[node][i] = 0
    mat[node][node] = 1
voltages = la.solve(mat,rhs)
    


