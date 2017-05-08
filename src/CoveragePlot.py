import numpy as np
import matplotlib.pyplot as plt
import os
import getopt
import sys

bottom = 8
truncate_to = 2  # this needs to be a power of two since we are plotting on a sqrt scale
units = 'M'
num_labels = 8
figure_letter = ''
#input_coverage_file = '../data/SNAP/coverage_1000.txt'

try:
	opts, args = getopt.getopt(sys.argv[1:], "hi:o:b:t:u:n:f:", ["Help=", "InputFile=", "OutputFile=", "Bottom=", "TruncateTo=", "Units=", "NumLabel=", "FigureLetter="])
except getopt.GetoptError:
	print 'Unknown option, call using: python CoveragePlot.py -i <InputFile> -o <OutputFile> -b <BottomOfPlotFloat> -t <TruncateCoverageToPowerOfTwoInt> -u <Units "M" or "K"> -n <NumLabelInt> -f <FigureLetter>'
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print 'python CoveragePlot.py -i <InputFile> -o <OutputFile> -b <BottomOfPlotFloat> -t <TruncateCoverageToPowerOfTwoInt> -u <Units "M" or "K"> -n <NumLabelInt> -f <FigureLetter>'
		sys.exit(2)
	elif opt in ("-i", "--InputFile"):
		input_coverage_file = arg
	elif opt in ("-o", "--OutputFile"):
		output_file = arg
	elif opt in ("-b", "--Bottom"):
		bottom = int(arg)
	elif opt in ("-t", "--TruncateTo"):
		truncate_to = int(arg)
	elif opt in ("-u", "--Units"):
		units = str(arg)
	elif opt in ("-n", "--NumLabel"):
		num_labels = int(arg)
	elif opt in ("-f", "--FigureLetter"):
		figure_letter = str(arg)

# FYI, this is about 10,000X times easier to do by hand (as I did here) than to use a "pre-packaged"
# solution like, eg. Circleator (http://jonathancrabtree.github.io/Circleator/). Just give that a go and
# you'll be running back to a matplotlib solution faster than you can say an unsavory four-letter word ;)



# Helper rotate function
def rotate(l, n):
	return l[-n:] + l[:-n]


radii = list()
left_windows = list()
right_windows = list()
fid = open(os.path.abspath(input_coverage_file), 'r')
for line in fid.readlines():
	line = line.strip()
	line_split = line.split()
	left_windows.append(int(line_split[1]))
	right_windows.append(int(line_split[2]))
	radii.append(float(line_split[3]))

fid.close()
N = len(radii)
radii = rotate(radii, -int(np.floor(len(radii) / float(4))))  # rotate 90 degress counter clockwise
radii = np.flip(np.array(radii), 0)  # make it go clockwise
radii = np.sqrt(radii)  # do it on a sqrt scale
# Truncate the radii to the given truncate_to
radii_truncated = list()
for val in radii:
	if val > truncate_to:
		val = truncate_to
	radii_truncated.append(val)
radii_truncated = np.array(radii_truncated)
radii = radii_truncated  # replace with truncated value
max_height = int(np.ceil(max(radii)))  # outer ring location

# Create outer ring labels in Kb or Mb
genome_length = max(right_windows)
genome_locations = np.floor(np.linspace(0, genome_length, num_labels))
labels = list()
if units == 'K':
	for i in range(len(genome_locations)):
		labels.append('%.1f Kbp' % (genome_locations[i] / float(1000)))
elif units == 'M':
	for i in range(len(genome_locations)):
		labels.append('%.1f Mbp' % (genome_locations[i] / float(1000000)))
elif units == 'bp':
	for i in range(len(genome_locations)):
		labels.append('%d bp' % genome_locations[i])
else:
	raise Exception('Unknown unit type. Pick one of "bp", "K" or "M"')

labels.reverse()  # make labels go clockwise
labels = rotate(labels, int(np.floor(len(labels) / float(4))))  # rotate 90 degrees counter clockwise
labels.insert(0, labels[-1])  # start at 0
labels = labels[0:num_labels]  # make correct length


theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)  # location of bars (evenly across 0->2*Pi)
width = (2*np.pi) / N  # width of bars

# Do all the plotting and decorating
plt.figure(figsize=(8, 8))
ax = plt.subplot(111, polar=True)  # make plot
bars = ax.bar(theta, radii, width=width, bottom=bottom)  # plot data
ax.set_xticklabels(labels)  # put the labels on
ax.set_thetagrids(np.linspace(0, 360, num_labels, endpoint=False), frac=1.03)  # make labels not so close

# put labels on inner/outer rings
ax.set_rgrids([bottom, bottom + max_height],  # plotted on sqrt scale
			labels=['0 X', '%.1f X' % (max_height**2)],  # but the actual value is ^2 of what's given
			angle=90, va='top', ha='center', color='b')
ax.grid(linewidth=0)  # get rid of all the grid lines except the outer one (clear where it starts)
ax.spines['polar'].set_alpha(0.1)  # make the rings lighter

# add small lines to designate the location of the labels
inner_radius = bottom
outer_radius = bottom + max_height
for angle in np.linspace(0, 360, num_labels, endpoint=False):
	# Recall that we are in polar, so (theta, r)
	#plt.plot([angle * np.pi / 180., angle * np.pi / 180.], [bottom, bottom + max_height], color='k', alpha=0.1)
	plt.plot([angle * np.pi / 180., angle * np.pi / 180.], [bottom + 0.75*max_height, bottom + max_height], color='k', alpha=0.1)

# Set the ylim so it doesn't scale to the newly plotted lines
plt.ylim([0, bottom + max_height])

# Figure letter
font = {'family': 'serif',
		'color':  'black',
		'weight': 'normal',
		'size': 18,
		}
if figure_letter != '':
	ax.text(-.1, 1, '%s)' % figure_letter, horizontalalignment='left', verticalalignment='bottom', fontdict=font, transform=ax.transAxes)

#plt.savefig('../Paper/Figs/CoveragePlot.png')
plt.savefig(os.path.abspath(output_file))
#plt.show()









# Old code in case it's useful

# Use custom colors and opacity
#for r, bar in zip(radii, bars):
#	bar.set_facecolor(plt.cm.jet(r / 10.))
#	bar.set_alpha(0.8)

#plt.tight_layout()

# Attempt at curved text
#i = 0
#label = "Genome Coverage"
#start_angle = 75
#end_angle = 105
#thetas = np.linspace(end_angle*np.pi/180., start_angle*np.pi/180., len(label))
#for letter in label:
#	plt.text(thetas[i], .9*bottom, letter, color='b', va='top', rotation=(thetas[i]*180./np.pi - 90), ha='center')
#	i += 1
