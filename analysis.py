from MDAnalysis.analysis.rdf import InterRDF
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc

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
    ]
) -> None:

    rc('font',**{'size':15, 'family':'sans-serif','sans-serif':['Arial']})

    fig, ax = plt.subplots()

    if not (str_title is None):
        ax.plot([], linestyle='', label=str_title)
    for int_id,str_label in enumerate(dict_data):
        array_rdf = np.loadtxt(dict_data[str_label])

        if list_linestyle is None:
            str_linestyle = None
        else:
            str_linestyle = list_linestyle[int_id]
        
        ax.plot(
            array_rdf[:,0],
            array_rdf[:,1],
            label = str_label,
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
