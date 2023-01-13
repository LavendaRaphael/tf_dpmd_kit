import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from matplotlib.lines import Line2D
import os
import matplotlib as mpl

def set_rcparam():

    mpl.rcParams['font.size'] = 7
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['font.family'] = 'Arial'

def inset_img(
    fig,
    ax,
    dict_img: dict,
    dict_arrow: dict = None,
    str_save: str = None,
    tup_size: tuple = None,
) -> None:

    set_rcparam()

    for str_img, tup_pos in dict_img.items():
        image = plt.imread(str_img)
        axin = ax.inset_axes(tup_pos)
        im = axin.imshow(image)
        #axin.axis('off')
        axin.set_xticks([])
        axin.set_yticks([])

    if not(dict_arrow is None):
        for tup_xy, tup_xytext in dict_arrow.items():
            ax.annotate(
                text = '',
                xy = tup_xy,
                xytext = tup_xytext,
                textcoords = ax.transAxes,
                arrowprops = dict(
                    arrowstyle= '-|>',
                    linestyle='--',
                    color = 'orange',
                    #linewidth = 2.0
                )
            )
    
    fig.set_tight_layout(True)
    if str_save:
        if tup_size:
            fig.set_size_inches(tup_size)
        fig.savefig(str_save, dpi=600)

def plt_subplots(
    dict_title: dict,
    dict2d_data: dict[dict],
    str_xlabel: str,
    str_ylabel: str,
    tup_xlim: tuple = (None, None),
    dict_ylim: dict = None,
    str_save: str = None,
    tup_size: tuple = None,
    dict_legend: dict = None,
    bool_scatter: bool = False,
    bool_error: bool = False,
    str_yaxisformat: str = None,
) -> None:
    
    set_rcparam()
    mpl.rcParams['font.size'] = 6

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
        for str_label, dict_data in dict2d_data.items():
            if not (str_key in dict_data):
                continue
            str_file = dict_data[str_key]
            np_data = np.loadtxt(str_file)
            if bool_scatter:
                ax.scatter( np_data[:,0], np_data[:,1], label=str_label, s=0.5)
            else:
                line, = ax.plot( np_data[:,0], np_data[:,1], label=str_label, lw=1)
                if bool_error:
                    ax.fill_between( np_data[:,0], np_data[:,1]-np_data[:,2], np_data[:,1]+np_data[:,2], alpha=0.5)

        if not(str_yaxisformat is None):
            ax.yaxis.set_major_formatter(plt.FormatStrFormatter(str_yaxisformat))

        if (dict_legend is None) or (str_key in dict_legend):
            ax.legend(
                frameon = False,
                handlelength = 1
            )

    axs[-1].set_xlabel(str_xlabel)
    axs[0].set_xlim(tup_xlim)
    
    fig.set_tight_layout(True)
    if str_save:
        if not tup_size is None:
            fig.set_size_inches(tup_size)
        fig.savefig(str_save, dpi=600)

def plt_compare_text(
    dict_data: dict,
    tup_xlim: tuple = None,
    tup_ylim: tuple = None,
    str_title: str = None,
    str_save: str = None,
    str_xlabel: str = None,
    str_ylabel: str = None,
    tup_size: tuple = None,
    dict_color: dict = None,
    bool_error: bool = False,
    float_lw: float = None,
) -> None:

    set_rcparam()

    fig, ax = plt_compare(
        dict_data = dict_data,
        str_xlabel = str_xlabel,
        str_ylabel = str_ylabel,
        str_save = None,
        tup_xlim = tup_xlim,
        tup_ylim = tup_ylim,
        dict_color = dict_color,
        bool_error = bool_error,
        float_lw = float_lw,
    )

    ax.text(
        x=0.35,
        y=0.9,
        s = str_title,
        transform=ax.transAxes
    )
    ax.legend(
        frameon = False,
        handlelength = 1.0,
        #labelspacing = 0.1,
    )

    fig.set_tight_layout(True)
    if str_save:
        if not(tup_size is None):
            fig.set_size_inches(tup_size)
        fig.savefig(str_save, dpi=600)

def dict_color_temperature(
    tup_colormap: tuple,
    dict_temperature: dict,
) -> dict:

    sm = plt.cm.ScalarMappable(cmap='coolwarm', norm=plt.Normalize(vmin=tup_colormap[0], vmax=tup_colormap[1]))
    dict_color = {}
    for str_key in dict_temperature:
        dict_color[str_key] = sm.to_rgba(dict_temperature[str_key])
        
    return dict_color

def plt_compare(
    dict_data: dict,
    str_xlabel: str,
    str_ylabel: str,
    str_save: str = None,
    tup_xlim: tuple = None,
    tup_ylim: tuple = None,
    bool_minzero: bool = False,
    bool_maxzero: bool = False,
    bool_minus: bool = False, 
    bool_legend: bool = True,
    tup_size: tuple = None,
    dict_color: dict = None,
    bool_error: bool = False,
    float_lw: float = None,
) -> None:

    set_rcparam()

    fig, ax = plt.subplots()

    if dict_color is None:
        dict_color = {}

    for str_label in dict_data:
        str_file = dict_data[str_label]
        if not os.path.isfile(str_file):
            print(str_file, 'Not found')
            continue
        print(str_file)

        if str_label in dict_color:
            color = dict_color[str_label]
        else:
            color = None

        np_data = np.loadtxt(str_file)

        np_data_y = np_data[:,1]
        if bool_minus:
            np_data_y *= -1

        if bool_minzero:
            np_data_y -= min(np_data_y)

        if bool_maxzero:
            np_data_y -= max(np_data_y)

        ax.plot( np_data[:,0], np_data_y, label=str_label, color=color, linewidth=float_lw)

        if bool_error:
            ax.fill_between( np_data[:,0], np_data_y-np_data[:,2], np_data_y+np_data[:,2], color=color, alpha=0.5)
    
    if bool_legend:
        ax.legend(
            frameon = False,
            handlelength = 1.0,
            labelspacing = 0.1,
        )
    ax.set_xlabel(str_xlabel)
    ax.set_ylabel(str_ylabel)
    ax.set_xlim(tup_xlim)
    ax.set_ylim(tup_ylim)

    fig.set_tight_layout(True)
    if str_save:
        if not(tup_size is None):
            fig.set_size_inches(tup_size)
        fig.savefig(str_save, dpi=600)

    return fig, ax

def plt_error(
    dict_data: dict,
    str_xlabel: str,
    str_ylabel: str,
    str_save: str = None,
    tup_size: tuple = None,
    tup_ylim: tuple = None,
    tup_xlim: tuple = None,
    str_title: str = None,
) -> None:

    set_rcparam()
    mpl.rcParams['font.size'] = 6

    fig, ax = plt.subplots()

    for str_label, str_marker in zip(dict_data, Line2D.filled_markers):
        str_file = dict_data[str_label]
        np_data = np.loadtxt(str_file, ndmin=2)
        if np_data.shape[1] > 2:
            yerr = np_data[:,2]
        else:
            yerr = None
        ax.errorbar(np_data[:,0], np_data[:,1], yerr=yerr, linestyle=':', marker=str_marker, label=str_label, capsize=2, lw=0.75, markersize=2)

    ax.set_xlabel(str_xlabel)
    ax.set_ylabel(str_ylabel)
    ax.set_xlim(tup_xlim)
    ax.set_ylim(tup_ylim)

    ax.legend(
        frameon = False,
        handlelength = 1.0,
        labelspacing = 0.2,
        title = str_title
    )

    fig.set_tight_layout(True)
    if not(tup_size is None):
        fig.set_size_inches(tup_size)
    if str_save:
        fig.savefig(str_save, dpi=600)
