#!/usr/bin/python3
import sys
import random

def convert(fd, ofd):
    for l1 in fd:
        l2 = fd.readline()
        if (random.randint(0,1)==0):
            ofd.write(l1)
        else:
            ofd.write(l2)

if __name__ == '__main__':
    convert(sys.stdin, sys.stdout)
