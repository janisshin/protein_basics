#!/usr/bin/python
# to make spacefills
from pymol import cmd
from sys import argv

my_argv = argv[1:]

cmd.load("C:/Users/janis/Desktop/resurrection/VDB/%s.vdb" % my_argv[0])
cmd.create("spaceFill", "all within 70 of index %s" % my_argv[1])
# add to creation the individual atom

cmd.show("surface", "spaceFill")
cmd.save("spacefills/%s.pdb" % my_argv[2], "spaceFill")

