
def main():

    file = 'traj.lammpstrj'
    interval = 1

    with open(file, 'r') as fp:
        fp.readline()
        frame0 = -interval
        for line in fp:
            frame1 = int(line.strip())
            if (frame1-frame0 != interval):
                print(frame0, frame1)
            for i in range(392):
                fp.readline()
            frame0 = frame1
main()
