import os
import shutil

def read_timestep(
    str_lmpin: str
) -> float:
    print('read:', str_lmpin)
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
    print('read:', str_plmin)
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
    print('read:', str_old)
    print('write:', str_new)
    with open(str_old, 'r') as file_old:
        with open(str_new, 'w') as file_new:
            for str_line in file_old:
                list_line = str_line.split()
                if list_line[0][0] != '#' and float(list_line[0]) > float_time:
                    break
                file_new.write(str_line)

def write_colvar(
    str_old: str,
    str_new: str,
    float_time: float,
):
    print('read:', str_old)
    print('write:', str_new)
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
    print('read:', str_old)
    print('write:', str_new)
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
    print('read:', str_old)
    with open(str_old, 'r') as file_open:
        for str_line in file_open:
            if str_line == f'ITEM: TIMESTEP\n':
                int_line = 0
                break
        for str_line in file_open:
            int_line += 1
            if str_line == f'ITEM: TIMESTEP\n':
                break
    print('write:', str_new)
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

int_restartstep = 23000000
str_dir = 'bk'
list_filecp = [
    'restart.1000000'
]
#---------------------------------------------
list_filecp.extend([
    f'restart.{int_restartstep}',
    'lmp.in',
    'plm.in'
])

for str_file in list_filecp:
    shutil.copy(os.path.join(str_dir, str_file), '.')
float_timestep = read_timestep('lmp.in')
float_time = int_restartstep * float_timestep
write_hills(
    str_old = os.path.join(str_dir, 'HILLS'),
    str_new = 'HILLS',
    float_time = float_time,
)
write_colvar(
    str_old = os.path.join(str_dir, 'COLVAR'),
    str_new = 'COLVAR',
    float_time = float_time,
)
write_log(
    str_old = os.path.join(str_dir, 'log'),
    str_new = 'log',
    int_restartstep = int_restartstep
)
write_traj(
    str_old = os.path.join(str_dir, 'traj.lammpstrj'),
    str_new = 'traj.lammpstrj',
    int_restartstep = int_restartstep
)

