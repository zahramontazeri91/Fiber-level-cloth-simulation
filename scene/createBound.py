
import argparse
import numpy as np
import struct, time, sys


def main(args):
    margin = args.margin
    size = 2*margin + 1
    eps = 1e-3

    with open(args.infile, mode='rb') as fin:
        assert fin.read(3) == 'VOL'
        assert fin.read(1) == '\x03'
        assert struct.unpack('I', fin.read(4)) == (1,)

        time.clock()
        sz = struct.unpack('3I', fin.read(12))
        assert struct.unpack('I', fin.read(4))[0] == 1

        aabbMin = np.array(struct.unpack('3f', fin.read(12)))
        aabbMax = np.array(struct.unpack('3f', fin.read(12)))

        sys.stdout.write('Loading: ')
        sys.stdout.flush()
        A0 = np.empty((sz[2], sz[1], sz[0]))
        for i in range(0, sz[2]):
            for j in range(0, sz[1]):
                A0[i, j, :] = np.asarray(struct.unpack('%df' % sz[0], fin.read(4*sz[0])))
            if i % 10 == 0:
                sys.stdout.write('\rLoading: %d/%d' %(i, sz[2]))
                sys.stdout.flush()
        sys.stdout.write('\rLoaded in %.2f secs.\n' % time.clock())

    A = np.zeros((sz[2], sz[1], sz[0]), dtype=np.int)
    for z in range(0, sz[2]):
        for y in range(0, sz[1]):
            for x in range(0, sz[0]):
                p = np.maximum([x - margin, y - margin, z - margin], [0, 0, 0])
                if np.max(A0[p[2] : p[2] + size, p[1] : p[1] + size, p[0] : p[0] + size]) > eps:
                    A[z, y, x] = 1

    extents = aabbMax - aabbMin
    unit = np.divide(extents, np.array(sz).astype(np.float))

    dx = [-1, 1, 0, 0, 0, 0]
    dy = [0, 0, -1, 1, 0, 0]
    dz = [0, 0, 0, 0, -1, 1]

    with open(args.o + '.obj', 'w') as fout:
        tot = 0
        for z in range(0, sz[2]):
            for y in range(0, sz[1]):
                for x in range(0, sz[0]):
                    if A[z, y, x] > 0:
                        for t in range(0, 6):
                            nx, ny, nz = x + dx[t], y + dy[t], z + dz[t]
                            if nx < 0 or nx >= sz[0] or ny < 0 or ny >= sz[1] or nz < 0 or nz >= sz[2] or A[nz, ny, nx] == 0:
                                corner = aabbMin + np.multiply([x, y, z], unit)

                                if t == 0:
                                    # -X
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0], corner[1], corner[2])
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0], corner[1], corner[2] + unit[2])
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0], corner[1] + unit[1], corner[2] + unit[2])
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0], corner[1] + unit[1], corner[2])
                                elif t == 1:
                                    # +X
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0] + unit[0], corner[1], corner[2])
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0] + unit[0], corner[1] + unit[1], corner[2])
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0] + unit[0], corner[1] + unit[1], corner[2] + unit[2])
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0] + unit[0], corner[1], corner[2] + unit[2])
                                elif t == 2:
                                    # -Y
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0], corner[1], corner[2])
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0] + unit[0], corner[1], corner[2])
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0] + unit[0], corner[1], corner[2] + unit[2])
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0], corner[1], corner[2] + unit[2])
                                elif t == 3:
                                    # +Y
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0], corner[1] + unit[1], corner[2])
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0], corner[1] + unit[1], corner[2] + unit[2])
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0] + unit[0], corner[1] + unit[1], corner[2] + unit[2])
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0] + unit[0], corner[1] + unit[1], corner[2])
                                elif t == 4:
                                    # -Z
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0], corner[1], corner[2])
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0], corner[1] + unit[1], corner[2])
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0] + unit[0], corner[1] + unit[1], corner[2])
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0] + unit[0], corner[1], corner[2])
                                else:
                                    assert t == 5
                                    # +Z
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0], corner[1], corner[2] + unit[2])
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0] + unit[0], corner[1], corner[2] + unit[2])
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0] + unit[0], corner[1] + unit[1], corner[2] + unit[2])
                                    print >> fout, 'v %.6f %.6f %.6f' % (corner[0], corner[1] + unit[1], corner[2] + unit[2])
                                tot += 4

        for i in range(0, tot, 4):
            print >> fout, 'f %d %d %d %d' % (i + 1, i + 2, i + 3, i + 4)

    sys.stdout.write('Processed in %.2f secs.\n' % time.clock())
    print('successfully converted volumes to obj files! \n')
    
if __name__ == "__main__":
    print ('generate obj mesh from volume AABB')
    parser = argparse.ArgumentParser(
            description='Create bounding meshes for single-channel Mitsuba volumes',
            epilog='Shuang Zhao (shz@ics.uci.edu)')

    parser.add_argument('infile', metavar='filename', type=str,
            help='Input Mitsuba volume in VOL format.')

    parser.add_argument('-o', metavar='filename', type=str, default='output',
            help="Output filename without extension. If not given, the default value is 'output'.")
    parser.add_argument('--margin', type=int, default=1,
            help='Margin. If not given, the default value is 1.')

    args = parser.parse_args()
    main(args)
