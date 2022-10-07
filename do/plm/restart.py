
int_restartstep = 1000000

def read_timestep(
    str_lmpin: str
) -> float:
    with open(str_lmpin, 'r') as file_open:
        for str_line in file_open:
            list_line = str_line.split()
            if list_line and list_line[0] == 'timestep':
                float_timestep = float(list_line[1])
                break
    return float_timestep

def read_stride(
    str_plmin: str
) -> int:
    with open(str_plmin, 'r') as file_open:
        for str_line in file_open:
            list_line = str_line.split()
            if list_line and list_line[0] == 'PRINT':
                for str_arg in list_line:
                    if str_arg[:7]=='STRIDE=':
                        int_stride = int(str_arg[7:])
                        break
    return int_stride

def write_hills(
    str_old: str,
    str_new: str,
    float_time: float,
):
    with open(str_old, 'r') as file_old:
        with open(str_new, 'w') as file_new:
            for str_line in file_old:
                list_line = str_line.split()
                if list_line[0][0] != '#' and float(list_line[0]) >= float_time:
                    break
                file_new.write(str_line)

def write_colvar(
    str_old: str,
    str_new: str,
    float_time: float,
):
    with open(str_old, 'r') as file_old:
        with open(str_new, 'w') as file_new:
            for str_line in file_old:
                list_line = str_line.split()
                if list_line[0][0] != '#' and float(list_line[0]) >= float_time:
                    break
                file_new.write(str_line)

def write_log(
    str_old: str,
    str_new: str,
    int_restartstep: int
):
    with open(str_old, 'r') as file_old:
        with open(str_new, 'w') as file_new:
            for str_line in file_old:
                list_line = str_line.split()
                if list_line and list_line[0] == '-' and list_line[1][:1]=='[' and float(list_line[1][1:-1]) >= int_restartstep:
                    break
                file_new.write(str_line)

def write_traj(
    str_old: str,
    str_new: str,
    int_restartstep: int
):
    with open(str_old, 'r') as file_open:
        for str_line in file_open:
            if str_line == f'ITEM: TIMESTEP\n':
                int_line = 0
                break
        for str_line in file_open:
            int_line += 1
            if str_line == f'ITEM: TIMESTEP\n':
                break

    with open(str_old, 'r') as file_old:
        with open(str_new, 'w') as file_new:
            for str_line in file_old:
                int_step = int(file_old.readline())
                if int_step >= int_restartstep:
                    break
                file_new.write(str_line)
                file_new.write(f'{int_step}\n')
                for i in range(int_line-2):
                    file_new.write(file_old.readline())

float_timestep = read_timestep()
float_time = int_restartstep * float_timestep
write_traj( int_restartstep )
