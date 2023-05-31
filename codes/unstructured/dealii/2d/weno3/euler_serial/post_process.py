# Open the files for reading and writing  

raw_file = open("plot_0177.dat", "r")
post_process_file = open("RHO_contour_0177.dat", "r")
final_file = open("plot_processed.dat", "w")

# Read the headers of the file and write them in the final file 

for i in range(0, 8):
    line = raw_file.readline()
    final_file.write(line)

# Now get the number of vertices and cells in the line 

wordList = line.split()

n_vertices = int((wordList[4]))
n_cells = int((wordList[8]))

# Now write all the x-coordinates of the vertices to the final file

for i in range(0, n_vertices):
    line = raw_file.readline()
    final_file.write(line)

# Read a blank line
line = raw_file.readline()
final_file.write(line)

# Now write all the x-coordinates of the vertices to the final file
for i in range(0, n_vertices):
    line = raw_file.readline()
    final_file.write(line)

# Read a blank line
line = raw_file.readline()
final_file.write(line)

# Now replace the node data with the new node data 
for i in range(0, n_vertices):
    line = raw_file.readline()
    line2 = post_process_file.readline()
    final_file.write(line2)

# Read a blank line
line = raw_file.readline()
final_file.write(line)

# Now add the cell vertex connectivity info 
for i in range(0, n_vertices):
    line = raw_file.readline()
    final_file.write(line)
 
# Close all the files  
raw_file.close()
post_process_file.close()
final_file.close()
