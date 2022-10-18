from matplotlib import pyplot as plt
import numpy as np
from matplotlib import rc

def gen_filesave(
    tup_elements: tuple,
) -> list:

    return f'rdf.{tup_elements[0]}.{tup_elements[1]}.converge.pdf'

def gen_filename(
     tup_elements: tuple,
     tup_id: tuple,
) -> list:

    return f'rdf.{tup_elements[0]}.{tup_elements[1]}.{tup_id[0]:07d}_{tup_id[-1]:07d}.csv'

def gen_label(
    float_timestep: float,
    tup_id: tuple,
):
    int_start = int(tup_id[0]*float_timestep)
    int_end = int(tup_id[-1]*float_timestep)
    return f'{int_start}-{int_end}ps'

def gen_listfile(
    tup_elements: tuple,
    list_id: list,
    float_timestep: float
) -> list:

    list2d_file = []
    for tup_id in list_id:
        list2d_file.append([gen_filename(tup_elements, tup_id), gen_label(float_timestep, tup_id)])
    return list2d_file

def def_plt(
    list2d_data,
    tup_xrange,
    tup_yrange,
    str_save,
    str_title
) -> None:

    rc('font',**{'size':15, 'family':'sans-serif','sans-serif':['Arial']})

    fig, ax = plt.subplots()
    
    for list_data in list2d_data:
        array_rdf = np.loadtxt( 
            fname = list_data[0],
        )
        ax.plot(
            array_rdf[:,0], 
            array_rdf[:,1], 
            label = list_data[1],
        )
    ax.legend()
    ax.set_title(str_title)
    ax.set_xlabel('r (Å)')
    ax.set_ylabel('RDF')
    ax.set_xlim(tup_xrange)
    ax.set_ylim(tup_yrange)
    fig.set_size_inches(8, 6)
    if str_save:
        fig.savefig(str_save, bbox_inches='tight')

# setup

list_id = []
#'''
list_id.append((      0, 50000))
list_id.append((  50000,100000))
list_id.append(( 100000,150000))
list_id.append(( 150000,200000))
#'''
float_timestep = 0.005

'''
list_id.append((1000,46000))
list_id.append((1000,55000))
float_timestep = 4.837769e-04
#'''

def_plt(
    list2d_data = gen_listfile(('c','o_w'), list_id, float_timestep),
    tup_xrange = (2,6),
    tup_yrange = (0,2.5),
    str_save = gen_filesave(('c','o_w')),
    str_title = r'C-O$_W$'
)
def_plt(
    list2d_data = gen_listfile(('o_1','h_w'), list_id, float_timestep),
    tup_xrange = (1,6),
    tup_yrange = (0,2),
    str_save = gen_filesave(('o_1','h_w')),
    str_title = r'$^=$O-H$_W$'
)
def_plt(
    list2d_data = gen_listfile(('o_0_2','h_w'), list_id, float_timestep),
    tup_xrange = (1,6),
    tup_yrange = (0,2),
    str_save = gen_filesave(('o_0_2','h_w')),
    str_title = r'O$_H$-H$_W$'
)
def_plt(
    list2d_data = gen_listfile(('h_0_1','o_w'), list_id, float_timestep),
    tup_xrange = (1,6),
    tup_yrange = (0,3),
    str_save = gen_filesave(('h_0_1','o_w')),
    str_title = r'H$_O$-O$_W$'
)
def_plt(
    list2d_data = gen_listfile(('o_w','o_w'), list_id, float_timestep),
    tup_xrange = (2,6),
    tup_yrange = (0,4),
    str_save = gen_filesave(('o_w','o_w')),
    str_title = r'O$_W$-O$_W$'
)
#'''
plt.show()

