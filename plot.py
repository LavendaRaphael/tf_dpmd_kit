import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from matplotlib.lines import Line2D
import os
import matplotlib as mpl
import pandas as pd
from matplotlib.patches import Arc
from matplotlib.transforms import IdentityTransform, TransformedBbox, Bbox

class AngleAnnotation(Arc):
    """
    Draws an arc between two vectors which appears circular in display space.
    """
    def __init__(self, xy, p1, p2, size=75, unit="points", ax=None,
                 text="", textposition="inside", text_kw=None, **kwargs):
        """
        Parameters
        ----------
        xy, p1, p2 : tuple or array of two floats
            Center position and two points. Angle annotation is drawn between
            the two vectors connecting *p1* and *p2* with *xy*, respectively.
            Units are data coordinates.

        size : float
            Diameter of the angle annotation in units specified by *unit*.

        unit : str
            One of the following strings to specify the unit of *size*:

            * "pixels": pixels
            * "points": points, use points instead of pixels to not have a
              dependence on the DPI
            * "axes width", "axes height": relative units of Axes width, height
            * "axes min", "axes max": minimum or maximum of relative Axes
              width, height

        ax : `matplotlib.axes.Axes`
            The Axes to add the angle annotation to.

        text : str
            The text to mark the angle with.

        textposition : {"inside", "outside", "edge"}
            Whether to show the text in- or outside the arc. "edge" can be used
            for custom positions anchored at the arc's edge.

        text_kw : dict
            Dictionary of arguments passed to the Annotation.

        **kwargs
            Further parameters are passed to `matplotlib.patches.Arc`. Use this
            to specify, color, linewidth etc. of the arc.

        """
        self.ax = ax or plt.gca()
        self._xydata = xy  # in data coordinates
        self.vec1 = p1
        self.vec2 = p2
        self.size = size
        self.unit = unit
        self.textposition = textposition

        super().__init__(self._xydata, size, size, angle=0.0,
                         theta1=self.theta1, theta2=self.theta2, **kwargs)

        self.set_transform(IdentityTransform())
        self.ax.add_patch(self)

        self.kw = dict(ha="center", va="center",
                       xycoords=IdentityTransform(),
                       xytext=(0, 0), textcoords="offset points",
                       annotation_clip=True)
        self.kw.update(text_kw or {})
        self.text = ax.annotate(text, xy=self._center, **self.kw)

    def get_size(self):
        factor = 1.
        if self.unit == "points":
            factor = self.ax.figure.dpi / 72.
        elif self.unit[:4] == "axes":
            b = TransformedBbox(Bbox.unit(), self.ax.transAxes)
            dic = {"max": max(b.width, b.height),
                   "min": min(b.width, b.height),
                   "width": b.width, "height": b.height}
            factor = dic[self.unit[5:]]
        return self.size * factor

    def set_size(self, size):
        self.size = size

    def get_center_in_pixels(self):
        """return center in pixels"""
        return self.ax.transData.transform(self._xydata)

    def set_center(self, xy):
        """set center in data coordinates"""
        self._xydata = xy

    def get_theta(self, vec):
        vec_in_pixels = self.ax.transData.transform(vec) - self._center
        return np.rad2deg(np.arctan2(vec_in_pixels[1], vec_in_pixels[0]))

    def get_theta1(self):
        return self.get_theta(self.vec1)

    def get_theta2(self):
        return self.get_theta(self.vec2)

    def set_theta(self, angle):
        pass

    # Redefine attributes of the Arc to always give values in pixel space
    _center = property(get_center_in_pixels, set_center)
    theta1 = property(get_theta1, set_theta)
    theta2 = property(get_theta2, set_theta)
    width = property(get_size, set_size)
    height = property(get_size, set_size)

    # The following two methods are needed to update the text position.
    def draw(self, renderer):
        self.update_text()
        super().draw(renderer)

    def update_text(self):
        c = self._center
        s = self.get_size()
        angle_span = (self.theta2 - self.theta1) % 360
        angle = np.deg2rad(self.theta1 + angle_span / 2)
        r = s / 2
        if self.textposition == "inside":
            r = s / np.interp(angle_span, [60, 90, 135, 180],
                                          [3.3, 3.5, 3.8, 4])
        self.text.xy = c + r * np.array([np.cos(angle), np.sin(angle)])
        if self.textposition == "outside":
            def R90(a, r, w, h):
                if a < np.arctan(h/2/(r+w/2)):
                    return np.sqrt((r+w/2)**2 + (np.tan(a)*(r+w/2))**2)
                else:
                    c = np.sqrt((w/2)**2+(h/2)**2)
                    T = np.arcsin(c * np.cos(np.pi/2 - a + np.arcsin(h/2/c))/r)
                    xy = r * np.array([np.cos(a + T), np.sin(a + T)])
                    xy += np.array([w/2, h/2])
                    return np.sqrt(np.sum(xy**2))

            def R(a, r, w, h):
                aa = (a % (np.pi/4))*((a % (np.pi/2)) <= np.pi/4) + \
                     (np.pi/4 - (a % (np.pi/4)))*((a % (np.pi/2)) >= np.pi/4)
                return R90(aa, r, *[w, h][::int(np.sign(np.cos(2*a)))])

            bbox = self.text.get_window_extent()
            X = R(angle, r, bbox.width, bbox.height)
            trans = self.ax.figure.dpi_scale_trans.inverted()
            offs = trans.transform(((X-s/2), 0))[0] * 72
            self.text.set_position([offs*np.cos(angle), offs*np.sin(angle)])

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
    str_save: str = None,
    tup_size: tuple = None,
    list_type: list = None,
) -> None:

    if not tup_size is None:
        fig.set_size_inches(tup_size)

    if str_save:
        if list_type is None:
            fig.savefig(str_save)
        else:
            for str_type in list_type:
                fig.savefig(f'{str_save}.{str_type}')

def add_text(
    ax,
    dict_text: dict,
) -> None:

    for str_text, tup_xy in dict_text.items():
        ax.text(
            x = tup_xy[0],
            y = tup_xy[1],
            s = str_text,
            horizontalalignment = 'center',
        )

def inset_img(
    ax,
    dict_img: dict,
    dict_spinecolor: dict = None,
    dict_spinels: dict = None,
    float_lw: float = None,
    bool_rot90: bool = False,
) -> None:

    if dict_spinecolor is None:
        dict_spinecolor = {}
    if dict_spinels is None:
        dict_spinels = {}

    for str_img, tup_pos in dict_img.items():
        axin = ax.inset_axes(tup_pos)
        image = plt.imread(str_img)
        if bool_rot90:
            image = np.rot90(np.array(image), k=3)
        axin.imshow(image)
        #axin.axis('off')
        axin.set_xticks([])
        axin.set_yticks([])
        for spine in axin.spines.values():
            spine.set_linewidth(float_lw)
            if str_img in dict_spinecolor:
                spine.set_edgecolor(dict_spinecolor[str_img])
            if str_img in dict_spinels:
                spine.set_linestyle(dict_spinels[str_img])

def add_arrow(
    ax,
    dict_arrow: dict = None,
    float_lw: float = None,
    str_arrowstyle: str = '-|>',
    str_color: str = None,
    str_xycoords: str = 'data',
    str_textcoords: str = None,
    str_ls: str = None,
):

    if str_textcoords is None:
        str_textcoords = str_xycoords

    for key, list_pos in dict_arrow.items():
        ax.annotate(
            text = '',
            xy = list_pos[1],
            xytext = list_pos[0],
            xycoords = str_xycoords,
            textcoords = str_xycoords,
            arrowprops = dict(
                arrowstyle = str_arrowstyle,
                linestyle = str_ls,
                color = dict_color[key],
                linewidth = float_lw,
                mutation_scale = 10,
            )
        )
    
def plt_subplots(
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
    
    return fig, ax

def add_text(
    ax,
    dict_text: dict,
):

    for str_text, tup_pos in dict_text.items():
        ax.text(
            x = tup_pos[0],
            y = tup_pos[1],
            s = str_text,
            transform=ax.transAxes
        )

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
    
    set_lw(ax, float_lw)

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
