import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from matplotlib.lines import Line2D
import os
import matplotlib as mpl
import pandas as pd
from matplotlib import patches

def set_rcparam():

    mpl.rcParams['font.size'] = 6
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['figure.dpi'] = 600
    #mpl.rcParams["figure.autolayout"] = True
    mpl.rcParams['figure.constrained_layout.use'] = True

def set_lw(
    ax,
    float_lw,
):
    for spin in ax.spines.values():
        spin.set_linewidth(float_lw)
    ax.tick_params(width=float_lw)

def plt_hist(
    dict_data: dict,
    str_xlabel: str,
    str_ylabel: str,
    tup_xlim: tuple = None,
    tup_ylim: tuple = None,
    bool_legend: bool = True,
    dict_color: dict = None,
    bool_error: bool = False,
    f_lw: float = None,
    float_scale: float = 1.0,
) -> None:

    fig, ax = plt.subplots()

    if dict_color is None:
        dict_color = {}

    for str_label in dict_data:
        str_file = dict_data[str_label]
        print(str_file)
        if not os.path.isfile(str_file):
            print('Not found')
            continue

        if str_label in dict_color:
            color = dict_color[str_label]
        else:
            color = None

        df_data = pd.read_csv(str_file)
        print(df_data)
        ax.hist(
            df_data*float_scale,
            label = str_label,
            bins = 'auto',
            density = True,
            color = color,
            histtype = 'step',
            lw = f_lw
        )

    set_lw(ax, f_lw)

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

    return fig, ax

def save(
    fig,
    file_save: str = None,
    tup_size: tuple = None,
    list_type: list = None,
) -> None:

    if not tup_size is None:
        fig.set_size_inches(tup_size)

    if file_save:
        if list_type is None:
            fig.savefig(file_save)
        else:
            for str_type in list_type:
                fig.savefig(f'{file_save}.{str_type}')

def add_text(
    ax,
    dict_text: dict,
    **kwargs,
) -> None:

    for tup_xy, str_text in dict_text.items():
        ax.text(
            x = tup_xy[0],
            y = tup_xy[1],
            s = str_text,
            **kwargs,
        )

def inset_img(
    ax,
    dict_img: dict,
    dict_spinecolor: dict = None,
    dict_spinels: dict = None,
    float_lw: float = None,
    bool_rot90: bool = False,
    bool_axis: bool = True,
    **kw
) -> None:

    if dict_spinecolor is None:
        dict_spinecolor = {}
    if dict_spinels is None:
        dict_spinels = {}

    axs = []
    for str_img, tup_pos in dict_img.items():
        axin = ax.inset_axes(tup_pos, **kw)
        axs.append(axin)
        image = plt.imread(str_img)
        if bool_rot90:
            image = np.rot90(np.array(image), k=3)
        axin.imshow(image)
        if not bool_axis:
            axin.axis('off')
        axin.set_xticks([])
        axin.set_yticks([])
        for spine in axin.spines.values():
            spine.set_linewidth(float_lw)
            if str_img in dict_spinecolor:
                spine.set_edgecolor(dict_spinecolor[str_img])
            if str_img in dict_spinels:
                spine.set_linestyle(dict_spinels[str_img])
    return axs

def add_line(
    ax,
    dict_line: dict = None,
    **kwargs,
):

    for key, list_pos in dict_line.items():
        ax.plot(
            [list_pos[0][0], list_pos[1][0]],
            [list_pos[0][1], list_pos[1][1]],
            **kwargs
        )
    
def add_arrow(
    ax,
    list_arrow: list = None,
    **kwargs,
):

    for list_pos in list_arrow:
        arrow = patches.FancyArrowPatch(
            list_pos[0],
            list_pos[1],
            **kwargs
        )
        ax.add_patch(arrow)
    
def plt_subplots(
    axs,
    dict_title: dict,
    dict2d_data: dict[dict],
    str_xlabel: str,
    str_ylabel: str,
    tup_xlim: tuple = (None, None),
    dict_ylim: dict = None,
    dict_legend: dict = None,
    bool_scatter: bool = False,
    bool_error: bool = False,
    str_yaxisformat: str = None,
    float_lw: float = None,
) -> None:
    
    int_nplot = len(dict_title)

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
                line, = ax.plot( np_data[:,0], np_data[:,1], label=str_label, lw=float_lw)
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
    ax,
    dict_data: dict,
    str_xlabel: str,
    str_ylabel: str,
    tup_xlim: tuple = None,
    tup_ylim: tuple = None,
    bool_minzero: bool = False,
    bool_maxzero: bool = False,
    bool_minus: bool = False, 
    bool_legend: bool = True,
    dict_color: dict = None,
    bool_error: bool = False,
    float_lw: float = None,
) -> None:

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
    
    set_lw(ax, float_lw)

    if bool_legend:
        ax.legend(
            frameon = False,
        )
    ax.set_xlabel(str_xlabel)
    ax.set_ylabel(str_ylabel)
    ax.set_xlim(tup_xlim)
    ax.set_ylim(tup_ylim)

def plt_error(
    dict_data: dict,
    str_xlabel: str,
    str_ylabel: str,
    tup_ylim: tuple = None,
    tup_xlim: tuple = None,
    str_title: str = None,
    float_lw: float = None,
    bool_error: bool = True,
    float_scale: float = 1.0,
    dict_ls: dict = None,
) -> None:

    fig, ax = plt.subplots()

    if dict_ls is None:
        dict_ls = {}

    for str_label, str_marker in zip(dict_data, Line2D.filled_markers):
        str_file = dict_data[str_label]
        np_data = np.loadtxt(str_file, ndmin=2)
        if bool_error:
            yerr = np_data[:,2]*float_scale
        else:
            yerr = None
        if str_label not in dict_ls:
            dict_ls[str_label] = ''
        ax.errorbar(np_data[:,0], np_data[:,1]*float_scale, yerr=yerr, linestyle=dict_ls[str_label], marker=str_marker, label=str_label, capsize=2, lw=float_lw, markersize=2)

    set_lw(ax, float_lw)

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

    return fig, ax
