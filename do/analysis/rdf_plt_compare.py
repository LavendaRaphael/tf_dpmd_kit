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
    str_save,
    str_title
) -> None:

    rc('font',**{'size':15, 'family':'sans-serif','sans-serif':['Arial']})

    fig, ax = plt.subplots()
    
    for list_data in list2d_data:
        array_rdf = np.loadtxt( 
            fname = list_data[0],
            delimiter = ','
        )
        ax.plot(
            array_rdf[:,0], 
            array_rdf[:,1], 
            label = list_data[1],
            linewidth = 2
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

# data

list_ref_2009_JPCB = ['/home/faye/research/202203_MDCarbonicAcid/record/ref/2009_JPCB_KumarPPadma', 'BLYP 310K 46Mol (2009 Kumar P. Padma)']
list_ref_2014_ChemComm = ['/home/faye/research/202203_MDCarbonicAcid/record/ref/2014_ChemComm_SandeepKReddy', 'BLYP+vdW 64Mol (2014 Sandeep K. Reddy)']
str_label_aimd = 'SCAN 330K 127Mol (Mine AIMD)'
str_label_dpmd = 'SCAN 330K 127Mol (Mine DPMD)'
str_path_aimd = '/home/faye/research_d/202203_MDCarbonicAcid/server/01.init/H2CO3_TT_H2O_126/rdf'
str_path_dpmd = '/home/faye/research_d/202203_MDCarbonicAcid/server/03.train/01.init_model/02.md/H2CO3_TT_H2O_126_CSVR/rdf'

#'''
def_rdf_plt_compare(
    list2d_data = [
        [os.path.join(str_path_aimd, 'rdf.c.o_w.0001000_0046000.csv'), str_label_aimd],
        [os.path.join(str_path_dpmd, 'rdf.c.o_w.0050000_0200000.csv'), str_label_dpmd],
        [os.path.join(list_ref_2009_JPCB[0], 'Fig_2.a.TT.csv'), list_ref_2009_JPCB[1]]
    ],
    tup_xrange = (2,6),
    tup_yrange = (0,2.5),
    str_save = f'rdf.c.o_w.ref.pdf',
    str_title = r'C-O$_W$'
)
def_rdf_plt_compare(
    list2d_data = [
        [os.path.join(str_path_aimd, 'rdf.o_1.h_w.0001000_0046000.csv'), str_label_aimd],
        [os.path.join(str_path_dpmd, 'rdf.o_1.h_w.0050000_0200000.csv'), str_label_dpmd],
        [os.path.join(list_ref_2009_JPCB[0], 'Fig_2.b.TT.csv'), list_ref_2009_JPCB[1]],
        [os.path.join(list_ref_2014_ChemComm[0], 'Fig_6.a.TT.csv'), list_ref_2014_ChemComm[1]]
    ],
    tup_xrange = (1,6),
    tup_yrange = (0,2),
    str_save = 'rdf.o_1.h_w.ref.pdf',
    str_title = r'$^=$O-H$_W$'
)
def_rdf_plt_compare(
    list2d_data = [
        [os.path.join(str_path_aimd, 'rdf.o_0_2.h_w.0001000_0046000.csv'), str_label_aimd],
        [os.path.join(str_path_dpmd, 'rdf.o_0_2.h_w.0050000_0200000.csv'), str_label_dpmd],
        [os.path.join(list_ref_2009_JPCB[0], 'Fig_2.c.TT.csv'), list_ref_2009_JPCB[1]]
    ],
    tup_xrange = (1,6),
    tup_yrange = (0,2),
    str_save = 'rdf.o_0_2.h_w.ref.pdf',
    str_title = r'O$_H$-H$_W$'
)
def_rdf_plt_compare(
    list2d_data = [
        [os.path.join(str_path_aimd, 'rdf.h_0_1.o_w.0001000_0046000.csv'), str_label_aimd],
        [os.path.join(str_path_dpmd, 'rdf.h_0_1.o_w.0050000_0200000.csv'), str_label_dpmd],
        [os.path.join(list_ref_2009_JPCB[0], 'Fig_2.d.TT.csv'), list_ref_2009_JPCB[1]],
        [os.path.join(list_ref_2014_ChemComm[0], 'Fig_6.b.TT.csv'), list_ref_2014_ChemComm[1]]
    ],
    tup_xrange = (1,6),
    tup_yrange = (0,3),
    str_save = 'rdf.h_0_1.o_w.ref.pdf',
    str_title = r'H$_O$-O$_W$'
)
def_rdf_plt_compare(
    list2d_data = [
        [os.path.join(str_path_aimd, 'rdf.o_w.o_w.0001000_0046000.csv'), str_label_aimd],
        [os.path.join(str_path_dpmd, 'rdf.o_w.o_w.0050000_0200000.csv'), str_label_dpmd],
        [os.path.join(list_ref_2009_JPCB[0], 'Fig_3.TT.csv'), list_ref_2009_JPCB[1]]
    ],
    tup_xrange = (2,6),
    tup_yrange = (0,4),
    str_save = 'rdf.o_w.o_w.ref.pdf',
    str_title = r'O$_W$-O$_W$'
)
#'''

plt.show() 
