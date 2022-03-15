import numpy as np


def create_off_file(length, w):
    h = w * np.sin(60 * np.pi / 180)
    with open('base.off', 'r', encoding='utf-8') as inp:
        with open('base2.off', 'w', encoding='utf-8') as out:
            for line in inp:
                if line.find('W') != -1:
                    line = line.replace('W', str(w))
                if line.find('S') != -1:
                    line = line.replace('S', str(h))
                if line.find('L') != -1:
                    line = line.replace('L', str(length))
                if line.find('P') != -1:
                    line = line.replace('P', str(w/2))
                out.write(line)
