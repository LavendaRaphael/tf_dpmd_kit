import os
import json
import tf_dpmd_kit.convert

list_mark = [
    #"CO2_H2O_127"               ,
    #"CO3_H3O_2_far_H2O_124"     ,
    #"CO3_H3O_2_near_H2O_124"    ,
    #"CO3_H3O_2_near_far_H2O_124",
    "H2CO3_CC_H2O_126"          ,
    "H2CO3_CT_H2O_126"          ,
    "H2CO3_TT_H2O_126"          ,
    "H3O_OH_far_H2O_126"        ,
    "H3O_OH_near_H2O_126"       ,
    "HCO3_H3O_far_H2O_125"      ,
    "HCO3_H3O_near_H2O_125"     ,
    #"H2O_128"                  ,
]

def def_task( str_type ):
    
    str_in = str_type + '.in'
    str_log = '$ESPRESSO_TMPDIR/' + str_type + '.log'

    dict_task_temp = {
        "command": "if [ -f ../please.stop ]; then true; else mpirun qe.cpbo_libxc.cp.x < "+str_in+" && if [ -f ../please.continue ]; then false; fi; fi",
        "forward_files": [
            str_in
        ],
        "backward_files": [],
        "outlog": str_log,
        "errlog": str_log
    }
    
    for str_dir in list_mark:
        os.chdir(str_dir)
        dict_task = dict_task_temp.copy()
        dict_task['task_work_path'] = str_dir
        dict_task['command'] = 'export ESPRESSO_TMPDIR=/public/spst/home/tianff/dpdispatcher/cpbo/TMPDIR/' + str_dir + '; '+ dict_task['command']
        with open('task.json', 'w') as open_json:
            json.dump( dict_task, fp=open_json, indent=4 )
        os.chdir('..')

def def_poscar2vccp_gs():
    for str_dir in list_mark:
        os.chdir(str_dir)
        tf_dpmd_kit.convert.def_poscar2vccp_gs()
        os.chdir('..')
def def_poscar2vccp_cpbo():
    for str_dir in list_mark:
        os.chdir(str_dir)
        tf_dpmd_kit.convert.def_poscar2vccp_cpbo()
        os.chdir('..')
def def_poscar2cp_gs():
    for str_dir in list_mark:
        os.chdir(str_dir)
        tf_dpmd_kit.convert.def_poscar2cp_gs()
        os.chdir('..')
def def_poscar2cp_cpbo():
    for str_dir in list_mark:
        os.chdir(str_dir)
        tf_dpmd_kit.convert.def_poscar2cp_cpbo()
        os.chdir('..')

os.chdir('..')
#def_poscar2cp_gs()
#def_task('gs')
def_poscar2cp_cpbo()
def_task('cp')
