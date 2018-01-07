#! /usr/bin/env python

#
# 1. Compute the central lines of fiber bundles without downsampling 
#    and store them in files.
# 2. Compute the forces of central line particles and store them in files.
#    The forces are summation of forces at fiber particles
#

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Compute the central line of a yarn and the forces long the yarn',
            epilog='Changxi Zheng (cxz@cs.columbia.edu)')
