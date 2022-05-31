"""
Script to handle backward compatibility for gmsh output file.
"""
import sys
import time
from sys import argv
import random
import copy
import math

if __name__ == "__main__":
	inFilename  = argv[1] 
	outFilename = argv[2]
	
	inFile  = open(inFilename, "rb")
	outFile = open(outFilename, "wb")

	print "? reading file: %s..." % inFilename
	line = inFile.readline()
	while "$Nodes" not in line:
	    line = inFile.readline()

	outFile.write("$NOD\n")
	num_vertices = int(inFile.readline())
	outFile.write(str(num_vertices) + "\n")

	for i in xrange(0,num_vertices):
		line = inFile.readline()
		outFile.write(line)

	inFile.readline()
	outFile.write("$ENDNOD\n")
	inFile.readline()
	outFile.write("$ELM\n")

	num_elem = int(inFile.readline())
	outFile.write(str(num_elem) + "\n")

	group = 0 
	zone = 0
	no_vertices = 0

	for i in xrange(0,num_elem):
		s = inFile.readline().split()

		shift = 3 + int(s[2])

		group = int(s[3])
		zone = group + 1

		#line
		if int(s[1]) == 1 :
			v1 = s[shift + 0]
			v2 = s[shift + 1]
			no_vertices = 2

			outFile.write(s[0] + " " + s[1] + " " + str(group) + " " + str(zone) + " " + str(no_vertices) + " " + v1 + " " + v2 + "\n")

		#triangle
		elif int(s[1]) == 2 :
			v1 = s[shift + 0]
			v2 = s[shift + 1]
			v3 = s[shift + 2]
			no_vertices = 3

			outFile.write(s[0] + " " + s[1] + " " + str(group) + " " + str(zone) + " " + str(no_vertices) + " " + v1 + " " + v2 +  " " + v3 + "\n")

		#quad
		elif int(s[1]) == 3 :
			v1 = s[shift + 0]
			v2 = s[shift + 1]
			v3 = s[shift + 2]
			v4 = s[shift + 3]
			no_vertices = 4

			outFile.write(s[0] + " " + s[1] + " " + str(group) + " " + str(zone) + " " + str(no_vertices) + " " + v1 + " " + v2 +  " " + v3 + " " + v4 + "\n")

		#tet
		elif int(s[1]) == 4 :
			v1 = s[shift + 0]
			v2 = s[shift + 1]
			v3 = s[shift + 2]
			v4 = s[shift + 3]
			no_vertices = 4

			outFile.write(s[0] + " " + s[1] + " " + str(group) + " " + str(zone) + " " + str(no_vertices) + " " + v1 + " " + v2 +  " " + v3 + " " + v4 + "\n")

		#hex
		elif int(s[1]) == 5 :
			v1 = s[shift + 0]
			v2 = s[shift + 1]
			v3 = s[shift + 2]
			v4 = s[shift + 3]
			v5 = s[shift + 4]
			v6 = s[shift + 5]
			v7 = s[shift + 6]
			v8 = s[shift + 7]
			no_vertices = 8

			outFile.write(s[0] + " " + s[1] + " " + str(group) + " " + str(zone) + " " + str(no_vertices) + " " + v1 + " " + v2 +  " " + v3 + " " + v4 + " " + v5 + " " + v6 + " " + v7 + "\n")


	inFile.close()
	outFile.write("$ENDELM\n")

	outFile.close()

