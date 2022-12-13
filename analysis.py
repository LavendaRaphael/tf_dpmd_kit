import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
from tf_dpmd_kit import plm
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
import MDAnalysis.analysis.msd

def msd(
    mda_u,
    str_select: str,
    float_dt: float,
    str_save: str,
    tup_snaprange: tuple = None,
):

    mda_msd = MDAnalysis.analysis.msd.EinsteinMSD(
        u = mda_u,
        select = str_select
    )

    if tup_snaprange is None:
        tup_snaprange = (0, len(mda_u.trajectory))

    mda_msd.run(
        start = tup_snaprange[0],
        stop = tup_snaprange[1],
        verbose = True,
    )
    
    int_nframes = mda_msd.n_frames
    np_final = np.zeros( shape=int_nframes, dtype=[('time', 'f4'), ('msd', 'f4')])
    np_final['time'] = np.arange(int_nframes)*float_dt # make the lag-time axis
    np_final['msd'] = mda_msd.results.timeseries
    np.savetxt(
        fname = str_save,
        X = np_final,
        header = ' '.join(np_final.dtype.names)
    )

def hbonds(
    universe,
    hydrogens_sel: str,
    acceptors_sel: str,
    str_save: str,
    donors_sel: str = None,
    update_selections: bool = False,
    float_ave: float = 1.0,
    d_a_cutoff: float = 3.0,
    d_h_a_angle_cutoff: float = 150,
    int_sperate: int = 10000,
) -> None:

    print(str_save)

    mda_hba = HydrogenBondAnalysis(
        universe = universe,
        hydrogens_sel = hydrogens_sel,
        acceptors_sel = acceptors_sel,
        donors_sel = donors_sel,
        d_a_cutoff = d_a_cutoff,
        d_h_a_angle_cutoff = d_h_a_angle_cutoff,
        update_selections = update_selections,
    )

    tup_snaprange = (0, len(universe.trajectory))

    np_final = np.zeros( shape=tup_snaprange[1], dtype=[('snap', 'i4'), ('nhbond', 'f4')])

    for start in range(tup_snaprange[0], tup_snaprange[1], int_sperate):
        stop = start + int_sperate
        if stop > tup_snaprange[1]:
            stop = tup_snaprange[1]
        print(f'{start} -> {stop}/{tup_snaprange[1]}')
        mda_hba.run(
            start = start,
            stop = stop,
            verbose = True,
        )

        np_final['snap'][start:stop] = mda_hba.times
        np_final['nhbond'][start:stop] = mda_hba.count_by_time()/float_ave

    np.savetxt(
        fname = str_save,
        X = np_final,
        header = ' '.join(np_final.dtype.names)
    )

def plt_subplots(
    dict_title: dict,
    dict2d_data: dict[dict],
    str_xlabel: str,
    str_ylabel: str,
    tup_xlim: tuple = (None, None),
    dict_ylim: dict = None,
    str_save: str = None,
    tup_size: tuple = None,
    bool_legend: bool = True,
    bool_scatter: bool = False,
) -> None:
    
    rc('font',**{'size':15, 'family':'sans-serif','sans-serif':['Arial']})

    int_nplot = len(dict_title)
    fig, axs = plt.subplots(int_nplot, 1, sharex='all')
    if int_nplot==1:
        axs = [axs]

    for ax, str_key in zip(axs, dict_title):
        if not(dict_ylim is None) and (str_key in dict_ylim):
            ax.set_ylim(dict_ylim[str_key])
        ax.set_ylabel(str_ylabel)
        ax.text(
            x=0.9,
            y=0.9,
            s = dict_title[str_key],
            horizontalalignment = 'right',
            verticalalignment = 'top',
            transform=ax.transAxes
        )
        for str_label,dict_data in dict2d_data.items():
            if not (str_key in dict_data):
                continue
            str_file = dict_data[str_key]
            np_data = np.loadtxt(str_file)
            if bool_scatter:
                ax.scatter( np_data[:,0], np_data[:,1], label=str_label, s=0.5)
            else:
                ax.plot( np_data[:,0], np_data[:,1], label=str_label)
    if bool_legend:
        axs[0].legend(loc='upper center') 
    axs[-1].set_xlabel(str_xlabel)
    axs[0].set_xlim(tup_xlim)
    if str_save:
        if not tup_size is None:
            fig.set_size_inches(tup_size)
        fig.savefig(str_save, bbox_inches='tight', dpi=300)

def plt_compare_text(
    dict_data: dict,
    tup_xlim: tuple = None,
    tup_ylim: tuple = None,
    str_title: str = None,
    str_save: str = None,
    str_xlabel: str = None,
    str_ylabel: str = None,
    tup_colormap: tuple = None,
    dict_temperature: dict = None,
    list_linestyle: list = None,
) -> None:

    rc('font',**{'size':15, 'family':'sans-serif','sans-serif':['Arial']}) 

    fig, ax = plm.plt_compare(
        dict_data = dict_data,
        str_xlabel = str_xlabel,
        str_ylabel = str_ylabel,
        str_save = None,
        tup_xlim = tup_xlim,
        tup_ylim = tup_ylim,
        tup_colormap = tup_colormap,
        dict_temperature = dict_temperature,
        list_linestyle = list_linestyle,
    )

    ax.text(
        x=0.9,
        y=0.9,
        s = str_title,
        horizontalalignment = 'right',
        verticalalignment = 'top',
        transform=ax.transAxes
    )
    ax.legend(loc='upper center')

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
