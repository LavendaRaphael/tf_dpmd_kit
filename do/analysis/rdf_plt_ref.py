from matplotlib import pyplot as plt
import matplotlib
import numpy as np
import os
from matplotlib import rc

# def

def def_rdf_plt_compare(
    list2d_data,
    tup_xrange,
    tup_yrange,
    str_title,
    str_save: str = None
) -> None:

    rc('font',**{'size':15, 'family':'sans-serif','sans-serif':['Arial']})

    list_linestyle = [
        'solid',
        'dashdot',
        'dotted',
        'dashed',
        (5,(10,3))
    ] 

    fig, ax = plt.subplots()

    ax.plot([], linestyle='', label=str_title)
    for list_data, str_linestyle in zip(list2d_data, list_linestyle):
        array_rdf = np.loadtxt(list_data[0])
        ax.plot(
            array_rdf[:,0],
            array_rdf[:,1],
            label = list_data[1],
            linewidth = 2,
            linestyle = str_linestyle
        )
    ax.legend()
    #ax.set_title(str_title)
    ax.set_xlabel('r (Å)')
    ax.set_ylabel('RDF')
    ax.set_xlim(tup_xrange)
    ax.set_ylim(tup_yrange)
    if str_save:
        #fig.set_size_inches(6, 5)
        fig.savefig(str_save, bbox_inches='tight')

# data

list_ref_2009_JPCB = ['/home/faye/research/202203_MDCarbonicAcid/record/ref/2009_JPCB_KumarPPadma/', 'BLYP 310K (2009 JPCB)']
list_ref_2010_JPCA = ['/home/faye/research/202203_MDCarbonicAcid/record/ref/2010_JPCA_LiuXiaodong/', 'BLYP 305K (2010 JPCA)']
list_ref_2014_ChemComm = ['/home/faye/research/202203_MDCarbonicAcid/record/ref/2014_ChemComm_SandeepKReddy/', 'BLYP+vdW (2014 ChemComm)']

str_label_aimd = 'SCAN 330K (AIMD)'
str_path_aimd = '/home/faye/research_d/202203_MDCarbonicAcid/server/01.init/H2CO3_TT_H2O_126/rdf/'
str_range = '0001000_0057877'

str_pair = 'c.o_w'
def_rdf_plt_compare(
    list2d_data = [
        [str_path_aimd+f'rdf.{str_pair}.{str_range}.csv', str_label_aimd],
        [list_ref_2009_JPCB[0]+'Fig_2.a.TT.csv', list_ref_2009_JPCB[1]]
    ],
    tup_xrange = (2,6),
    tup_yrange = (0,3.0),
    str_save = f'rdf.{str_pair}.ref.pdf',
    str_title = r'C-O$_W$'
)

str_pair = 'o_1.h_w'
def_rdf_plt_compare(
    list2d_data = [
        [str_path_aimd+f'rdf.{str_pair}.{str_range}.csv', str_label_aimd],
        [list_ref_2009_JPCB[0]+'Fig_2.b.TT.csv', list_ref_2009_JPCB[1]],
        [list_ref_2014_ChemComm[0]+'Fig_6.a.TT.csv', list_ref_2014_ChemComm[1]]
    ],
    tup_xrange = (1,6),
    tup_yrange = (0,2.5),
    str_save = f'rdf.{str_pair}.ref.pdf',
    str_title = r'$^=$O-H$_W$'
)

str_pair = 'o_0_2.h_w'
def_rdf_plt_compare(
    list2d_data = [
        [str_path_aimd+f'rdf.{str_pair}.{str_range}.csv', str_label_aimd],
        [list_ref_2009_JPCB[0]+'Fig_2.c.TT.csv', list_ref_2009_JPCB[1]],
        [list_ref_2010_JPCA[0]+'Fig_2.b.TT.csv', list_ref_2010_JPCA[1]]
    ],
    tup_xrange = (1,6),
    tup_yrange = (0,2.5),
    str_save = f'rdf.{str_pair}.ref.pdf',
    str_title = r'O$_H$-H$_W$'
)

str_pair = 'h_0_1.o_w'
def_rdf_plt_compare(
    list2d_data = [
        [str_path_aimd+f'rdf.{str_pair}.{str_range}.csv', str_label_aimd],
        [list_ref_2009_JPCB[0]+'Fig_2.d.TT.csv', list_ref_2009_JPCB[1]],
        [list_ref_2010_JPCA[0]+'Fig_2.a.TT.csv', list_ref_2010_JPCA[1]],
        [list_ref_2014_ChemComm[0]+'Fig_6.b.TT.csv', list_ref_2014_ChemComm[1]]
    ],
    tup_xrange = (1,6),
    tup_yrange = (0,3),
    str_save = f'rdf.{str_pair}.ref.pdf',
    str_title = r'H$_O$-O$_W$'
)

str_pair = 'o_w.o_w'
def_rdf_plt_compare(
    list2d_data = [
        [str_path_aimd+f'rdf.{str_pair}.{str_range}.csv', str_label_aimd],
        [list_ref_2009_JPCB[0]+'Fig_3.TT.csv', list_ref_2009_JPCB[1]]
    ],
    tup_xrange = (2,6),
    tup_yrange = (0,3.5),
    str_save = f'rdf.{str_pair}.ref.pdf',
    str_title = r'O$_W$-O$_W$'
)

plt.show() 
