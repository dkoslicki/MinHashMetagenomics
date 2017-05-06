# This will generate a figure depicting the general idea of the containment min hash
import shapely.geometry as sg
import matplotlib.pyplot as plt
import descartes
import numpy as np
import os
import subprocess

# global parameters
num_points_both = 100  # Number of points to sample from union (classical approach)
num_points_single = 50  # Number of points to sample from A (containment approach)

alpha = .4  # opacity of circles
x_big = -.5  # x-coordinate of bigger circle
radius_big = 1.4  # radius of bigger circle
x_small = 0.9  # x-coordinate of smaller circle
radius_small = .3  # radius of smaller circle
point_size_big = 0.01
point_size_small = 0.006

# Both sets together
plt.figure()
# create the circles with shapely
a = sg.Point(x_big, 0).buffer(radius_big)
b = sg.Point(x_small, 0).buffer(radius_small)

# Create random points
points = list()
# location of edges of circles
right_point_big = x_big + radius_big
left_point_big = x_big - radius_big
right_point_small = x_small + radius_small
left_point_small = x_small - radius_small
num_int_both = 0  # keep track of number of points in the intersection
for i in range(num_points_both):
	x = (right_point_small - left_point_big) * np.random.rand() + left_point_big
	y = 2 * radius_big * np.random.rand() - radius_big
	while (x - x_big)**2 + y**2 > radius_big**2 and (x - x_small)**2 + y**2 > radius_small**2:
		x = (right_point_small - left_point_big) * np.random.rand() + left_point_big
		y = 2 * radius_big * np.random.rand() - radius_big
	if (x - x_big) ** 2 + y ** 2 <= radius_big ** 2 and (x - x_small) ** 2 + y ** 2 <= radius_small ** 2:
		num_int_both += 1
	points.append(sg.Point(x, y).buffer(point_size_big))

# compute the 3 parts
left = a.difference(b)
right = b.difference(a)
middle = a.intersection(b)

# use descartes to create the matplotlib patches
font = {'family': 'serif',
		'color':  'black',
		'weight': 'normal',
		'size': 18,
		}
ax = plt.gca()
# Show the underlying sets
ax.add_patch(descartes.PolygonPatch(left, fc='r', ec='k', alpha=alpha))
ax.add_patch(descartes.PolygonPatch(right, fc='b', ec='k', alpha=alpha))
ax.add_patch(descartes.PolygonPatch(middle, fc='g', ec='k', alpha=alpha))

# Show the random points
for elem in points:
	ax.add_patch(descartes.PolygonPatch(elem, fc='k', ec='k', alpha=1))

# Label the points
ax.text(x_big, radius_big+0.02, r'$B$', horizontalalignment='left', verticalalignment='bottom', fontdict=font, usetex=True)
ax.text(x_small+radius_small, 0, r'$A$', horizontalalignment='left', verticalalignment='bottom', fontdict=font, usetex=True)

# control display
ax.set_xlim(x_big-radius_big-0.2, radius_big+0.2)
ax.set_ylim(-radius_big-0.3, radius_big+0.3)
ax.set_aspect('equal')
ax.text(-.1, 1, 'a)', horizontalalignment='left', verticalalignment='bottom', fontdict=font, transform=ax.transAxes)
plt.savefig('../Paper/Figs/ClassicalConceptual.png')


# Just the smaller set
plt.figure()
# create the circles with shapely
a = sg.Point(-.5, 0).buffer(1.4)
b = sg.Point(0.9, 0).buffer(.3)

# Create random points
points = list()
num_int_single = 0
for i in range(num_points_single):
	x = (right_point_small - left_point_small) * np.random.rand() + left_point_small
	y = 2 * radius_small * np.random.rand() - radius_small
	while (x-.9)**2 + y**2 > radius_small**2:
		x = (right_point_small - left_point_small) * np.random.rand() + left_point_small
		y = 2 * radius_small * np.random.rand() - radius_small
	if (x + .5) ** 2 + y ** 2 <= radius_big ** 2 and (x - .9) ** 2 + y ** 2 <= radius_small ** 2:
		num_int_single += 1
	points.append(sg.Point(x, y).buffer(point_size_small))

# compute the 3 parts
left = a.difference(b)
right = b.difference(a)
middle = a.intersection(b)

# use descartes to create the matplotlib patches
font = {'family': 'serif',
		'color':  'black',
		'weight': 'normal',
		'size': 18,
		}
ax = plt.gca()
# Show the underlying sets
ax.add_patch(descartes.PolygonPatch(right, fc='b', ec='k', alpha=alpha))
ax.add_patch(descartes.PolygonPatch(middle, fc='g', ec='k', alpha=alpha))

# Show the random points
for elem in points:
	ax.add_patch(descartes.PolygonPatch(elem, fc='k', ec='k', alpha=1))

# Label the points
ax.text(radius_small*np.cos(3.1415/2.+3.1415/4.)+x_small, radius_small*np.sin(3.1415/2.+3.1415/4.), r'$A\cap B$', horizontalalignment='right', verticalalignment='bottom', fontdict=font, usetex=True)
ax.text(x_small+radius_small, 0, r'$A$', horizontalalignment='left', verticalalignment='bottom', fontdict=font, usetex=True)

# control display
ax.set_xlim((x_small-radius_small-0.2)/((radius_small+0.2) - (-radius_small-0.2)), (x_small+radius_small+0.2)/((radius_small+0.2) - (-radius_small-0.2)))
ax.set_ylim(-radius_small-0.2, radius_small+0.2)
ax.set_aspect('equal')
ax.text(-.1, 1, 'b)', horizontalalignment='left', verticalalignment='bottom', fontdict=font, transform=ax.transAxes)
plt.savefig('../Paper/Figs/ContainmentConceptual.png')

print('Number in intersection when using both sets: %d' % num_int_both)
print('Number in intersection when using one set: %d' % num_int_single)
fid = open(os.path.abspath('../Paper/ClassicalConceptual.txt'), 'w')
fid.write("%d" % num_int_both)
fid.close()
fid = open(os.path.abspath('../Paper/ContainmentConceptual.txt'), 'w')
fid.write("%d" % num_int_single)
fid.close()

cmd = 'ls ../Paper/Figs/*.png | xargs -I{} convert {} -trim {}'
subprocess.check_output(cmd, shell = True)