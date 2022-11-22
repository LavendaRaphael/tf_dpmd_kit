from MDAnalysis.analysis.rdf import InterRDF
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc

def plt_subplots(
    dict_subplot: dict,
    dict_data: dict,
    str_xlabel: str,
    str_ylabel: str,
    tup_xlim: tuple,
    str_save: str = None,
    tup_size: tuple = None,
    bool_legend: bool = True,
) -> None:
    
    rc('font',**{'size':15, 'family':'sans-serif','sans-serif':['Arial']})

    int_nplot = len(dict_subplot)
    fig, axs = plt.subplots(int_nplot, 1, sharex='all')
    if int_nplot==1:
        axs = [axs]

    for ax, str_title in zip(axs, dict_subplot):
        tup_ylim = dict_subplot[str_title]
        ax.set_ylim(tup_ylim)
        ax.set_ylabel(str_ylabel)
        ax.text(
            x=0.9,
            y=0.9,
            s = str_title,
            horizontalalignment = 'right',
            verticalalignment = 'top',
            transform=ax.transAxes
        )
        for str_label,dict_tmp in dict_data.items():
            if not (str_title in dict_tmp):
                continue
            str_file = dict_tmp[str_title]
            np_data = np.loadtxt(str_file)
            ax.plot( np_data[:,0], np_data[:,1], label=str_label)
    if bool_legend:
        axs[0].legend(loc='upper center') 
    axs[-1].set_xlabel(str_xlabel)
    axs[0].set_xlim(tup_xlim)
    if str_save:
        if not tup_size is None:
            fig.set_size_inches(tup_size)
        fig.savefig(str_save, bbox_inches='tight', dpi=300)

def rdf_plt_compare(
    dict_data: dict,
    tup_xrange: tuple = None,
    tup_yrange: tuple = None,
    str_title: str = None,
    str_save: str = None,
    list_linestyle: list = [
        'solid',
        'dashdot',
        'dotted',
        'dashed',
        (5,(10,3))
    ],
    tup_colormap: tuple =  None,
) -> None:

    rc('font',**{'size':15, 'family':'sans-serif','sans-serif':['Arial']})

    fig, ax = plt.subplots()

    if tup_colormap:
        sm = plt.cm.ScalarMappable(cmap='coolwarm', norm=plt.Normalize(vmin=tup_colormap[0], vmax=tup_colormap[1]))

    if not (str_title is None):
        ax.plot([], linestyle='', label=str_title)
    for int_id,str_label in enumerate(dict_data):
        array_rdf = np.loadtxt(dict_data[str_label])

        if list_linestyle is None:
            str_linestyle = None
        else:
            str_linestyle = list_linestyle[int_id]

        if tup_colormap:
            color = sm.to_rgba(int(str_label[:3]))
        else:
            color = None
        
        ax.plot(
            array_rdf[:,0],
            array_rdf[:,1],
            label = str_label,
            linewidth = 2,
            linestyle = str_linestyle,
            color = color
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

def rdf(
    mda_atomgroup_0,
    mda_atomgroup_1,
    tup_snaprange: tuple,
    str_save: str,
    tup_rrange: tuple,
    int_nbins: int = 200  
) -> None:

    print(str_save)
    
    mda_rdf = InterRDF( 
        mda_atomgroup_0,
        mda_atomgroup_1, 
        nbins = int_nbins,
        range = tup_rrange,
    )

    mda_rdf.run(
        start = tup_snaprange[0],
        stop = tup_snaprange[1],
        verbose = True,
    )
    array_final = np.zeros( shape=(int_nbins), dtype=[('r(Angstrom)', 'f4'), ('RDF', 'f4')])
    array_final['r(Angstrom)'] = mda_rdf.results.bins
    array_final['RDF'] = mda_rdf.results.rdf
    
    np.savetxt(
        fname = str_save,
        X = array_final,
        header = ' '.join(array_final.dtype.names)
    )
