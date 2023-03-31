import time

def main():

    time_0 = time.time()

    file = 'traj.lammpstrj'
    file_new = 'traj.lammpstrj.new'
    interval = 20
    natoms = 384

    nline = natoms+7
    with open(file, 'r') as fp:
        with open(file_new, 'w') as fp_new:
            for line in fp:
                line = fp.readline()
                if (int(line)%interval == 0):
                    fp_new.write('ITEM: TIMESTEP\n' + line)
                    for i in range(nline):
                        fp_new.write(fp.readline())
                else:
                    for i in range(nline):
                        fp.readline()

    time_inter = time.time()-time_0
    print(time_inter)

main()
