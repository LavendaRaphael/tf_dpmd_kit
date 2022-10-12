import dpdata
import glob

def post_fp_pwscf():

    '''
    Modified from: https://github.com/deepmodeling/dpgen/blob/master/dpgen/generator/run.py#L3180-L3223
    '''

    fp_tasks = glob.glob( 'task.*')
    fp_tasks.sort()
    if len(fp_tasks) == 0 :
        return

    system_index = []
    for ii in fp_tasks :
        system_index.append(ii.split('.')[1])
    system_index.sort()
    set_tmp = set(system_index)
    system_index = list(set_tmp)
    system_index.sort()

    for ss in system_index :
        sys_output = glob.glob("task.%s.*/output"%ss)
        sys_input = glob.glob("task.%s.*/input"%ss)
        sys_output.sort()
        sys_input.sort()

        flag=True
        for ii,oo in zip(sys_input,sys_output) :
            if flag:
                _sys = dpdata.LabeledSystem(oo, fmt = 'qe/pw/scf')
                if len(_sys)>0:
                   all_sys=_sys
                   flag=False
                else:
                   pass
            else:
                _sys = dpdata.LabeledSystem(oo, fmt = 'qe/pw/scf')
                if len(_sys)>0:
                   all_sys.append(_sys)

        sys_data_path = 'data.%s'%ss
        print(sys_data_path, len(all_sys))
        all_sys.to_deepmd_raw(sys_data_path)
        all_sys.to_deepmd_npy(sys_data_path, set_size = len(sys_output))

post_fp_pwscf()
