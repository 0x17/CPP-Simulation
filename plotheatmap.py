from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt


def extract_points_from_file(filename):
	with open(filename, 'r') as fp:
		lines = fp.readlines()
		xs = list(map(lambda line: int(line.split(';')[0]), lines[1:]))
		ys = list(map(lambda line: int(line.split(';')[1]), lines[1:]))
		zs = list(map(lambda line: float(line.split(';')[2].strip()), lines[1:]))
		return xs, ys, zs


def figure_for_swarm_iteration(iteration_filename, ix):
	fig = plt.figure(ix)
	fig.clear()
	ax = fig.gca(projection='3d')
	ax.set_xlabel('b2')
	ax.set_ylabel('b3')
	ax.set_zlabel('obj')

	xs, ys, zs = extract_points_from_file('fullEnumObjectives.txt')
	ax.plot_trisurf(xs, ys, zs, cmap=cm.jet, linewidth=0.2, alpha=0.2)

	xs, ys, zs = extract_points_from_file(iteration_filename)
	ax.scatter(xs, ys, zs)

	ax.set_xbound(lower=0)
	ax.set_ybound(lower=0)
	ax.set_zbound(lower=0)

	ax.view_init(elev=50, azim=80)
	plt.savefig(iteration_filename[:-4] + '.pdf')


ix = 0
for iterationFilename in (['initialSwarm.txt'] + ['swarmIteration' + str(i) + '.txt' for i in range(1, 6)]):
	figure_for_swarm_iteration(iterationFilename, ix)
	ix += 1

#plt.show()
