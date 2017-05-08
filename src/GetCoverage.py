import sys

window_size = int(sys.argv[1])  # first argument is the window size (such as 10000)
input_file = sys.argv[2]  # second argument is input file (typically use pipes and /dev/fd/0)
output_file = sys.argv[3]  # third argument is output file (for use in MakeCoveragePlot.py input)

# Example call: samtools depth -a --reference PRJNA274798.fa out.sorted.sam | python get_coverage.py 10000 /dev/fd/0 /dev/fd/1 > coverage_10000.txt

val = 0
location_start = 1
location_end = 1
it = 1
out_fid = open(output_file, 'w')
for line in open(input_file, 'r').readlines():
	line_striped = line.strip()
	line_split = line_striped.split()
	header = line_split[0]
	location_end += 1
	val += int(line_split[-1])
	if it % window_size == 0:
		out_fid.write("%s\t%d\t%d\t%f\n" % (header, location_start, location_end, val / float(window_size)))
		location_start = location_end
		val = 0
	it += 1

if it < window_size:
	raise Exception("Window size is too large: genome is only %d bp long" % it)

out_fid.close()
