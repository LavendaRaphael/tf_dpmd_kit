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
    #mpl.rc('text', usetex=True)
    #mpl.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"

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

    texts = []
    for tup_xy, str_text in dict_text.items():
        text = ax.text(
            x = tup_xy[0],
            y = tup_xy[1],
            s = str_text,
            **kwargs,
        )
        texts.append(text)
    return texts

def inset_img(
    ax,
    dict_img: dict,
    spinecolor: dict = None,
    spinels: dict = None,
    spinealpha: float = None,
    spinelw: float = None,
    img_rot90: bool = False,
    axin_axis: bool = True,
    **kw
) -> None:

    axs = []
    for img, tup_pos in dict_img.items():
        axin = ax.inset_axes(tup_pos, **kw)
        axs.append(axin)
        image = plt.imread(img)
        if img_rot90:
            image = np.rot90(np.array(image), k=3)
        axin.imshow(image)
        if not axin_axis:
            axin.axis('off')
        axin.set_xticks([])
        axin.set_yticks([])
        for spine in axin.spines.values():
            spine.set_linewidth(spinelw)
            if not(spinecolor is None):
                if isinstance(spinecolor, dict):
                    if img in spinecolor:
                        spine.set_edgecolor(spinecolor[img])
                else:
                    spine.set_edgecolor(spinecolor)
            if not(spinels is None):
                if isinstance(spinels, dict):
                    if img in spinels:
                        spine.set_linestyle(spinels[img])
                else:
                    spine.set_linestyle(spinels)
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
    
def dict_color_temperature(
    colormap: tuple,
    dict_temperature: dict,
) -> dict:

    sm = plt.cm.ScalarMappable(cmap='coolwarm', norm=plt.Normalize(vmin=colormap[0], vmax=colormap[1]))
    dict_color = {}
    for key in dict_temperature:
        dict_color[key] = sm.to_rgba(dict_temperature[key])
        
    return dict_color

def plt_compare(
    ax,
    dict_data: dict,
    bool_minzero: bool = False,
    bool_maxzero: bool = False,
    bool_minus: bool = False, 
    dict_color: dict = None,
    bool_error: bool = False,
    float_lw: float = None,
) -> None:

    if dict_color is None:
        dict_color = {}

    for str_label in dict_data:
        if bool_error:
            x, y, e = dict_data[str_label]
        else:
            x, y = dict_data[str_label]

        if str_label in dict_color:
            color = dict_color[str_label]
        else:
            color = None

        if bool_minus:
            y *= -1

        if bool_minzero:
            y -= min(y)

        if bool_maxzero:
            y -= max(y)

        ax.plot( x, y, label=str_label, color=color, linewidth=float_lw)

        if bool_error:
            ax.fill_between( x, y-e, y+e, color=color, alpha=0.5)
    
    set_lw(ax, float_lw)

def plt_error(
    ax,
    dict_data: dict,
    float_lw: float = None,
    bool_error: bool = True,
    float_scale: float = 1.0,
    dict_ls: dict = None,
) -> None:

    if dict_ls is None:
        dict_ls = {}

    for label, marker in zip(dict_data, Line2D.filled_markers):
        if bool_error:
            x, y, e = dict_data[label]
            e *= float_scale
        else:
            x, y = dict_data[label]
            e = None
        if label not in dict_ls:
            dict_ls[label] = ''
        ax.errorbar(x, y*float_scale, e, linestyle=dict_ls[label], marker=marker, label=label, capsize=2, lw=float_lw, markersize=4)

    set_lw(ax, float_lw)

