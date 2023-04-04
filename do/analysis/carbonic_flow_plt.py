import matplotlib.pyplot as plt
from tf_dpmd_kit import analysis
from tf_dpmd_kit import plot
import matplotlib as mpl
import numpy as np

def run(ax):
    ax.axis('off')

    # horizon fig
    dx = 8/7*0.3
    dy = 5/6*0.3
    # vertical fig
    dx1 = dy*6/7
    dy1 = dx*7/6

    # right, top
    r = np.array([ dx/2,    0])
    t = np.array([    0,   dy/2])
    r1 = np.array([ dx1/2,    0])
    t1 = np.array([    0,   dy1/2])
    # shift
    sx = np.array([0.05,   0])
    sy = np.array([   0,0.05])

    # pos
    p_cc = np.array([0.4, dy/2])
    p_ct = np.array([0.4, 0.5])
    p_tt = np.array([0.4, 1-dy/2])
    p_xx = np.array([1-dx1/2, 0.5])
    # middle
    p_cc_ct = (p_cc+p_ct)/2
    p_ct_tt = (p_ct+p_tt)/2
    p_xx_tt = (p_xx-r1+t1+p_tt+r)/2
    p_xx_ct = (p_xx-r1+p_ct+r)/2
    p_xx_cc = (p_xx-r1-t1+p_cc+r)/2

    # img
    dir_cp  = '/home/faye/research_d/202203_MDCarbonicAcid/server/01.init/H2CO3_TT_H2O_126/plm/'
    dir_cc  = '/home/faye/research_d/202203_MDCarbonicAcid/server/04.md_nvt_velocity/330K/CC/plm/'
    [axin3] = plot.inset_img(
        ax,
        dict_img = {
            dir_cp +   '60147.png': (p_xx[0]-0.5*dx1, p_xx[1]-dy1/2, dx1, dy1),
        },
        dict_spinecolor = {
            dir_cp +   '60147.png': 'tab:grey',
        },
        bool_rot90 = True,
    )
    [axin0, axin1, axin2] = plot.inset_img(
        ax,
        dict_img = {
            dir_cc +'1.100001.png': (p_tt[0]-0.5*dx, p_tt[1]-0.5*dy, dx, dy),
            dir_cp +   '60281.png': (p_ct[0]-0.5*dx, p_ct[1]-0.5*dy, dx, dy),
            dir_cc +'0.003922.png': (p_cc[0]-0.5*dx, p_cc[1]-0.5*dy, dx, dy),
        },
        dict_spinecolor = {
            dir_cc +'1.100001.png': 'tab:blue',
            dir_cp +   '60281.png': 'tab:orange', 
            dir_cc +'0.003922.png': 'tab:green',
        }
    )

    # arrow
    plot.add_arrow(
        ax,
        list_arrow = [
            [p_tt-t-sx, p_ct+t-sx],
            [p_tt+r +sy,   p_xx+t1],
        ],
        arrowstyle = 'simple, head_length=2, head_width=2, tail_width=0.2',
        color = 'tab:blue',
    )
    plot.add_arrow(
        ax,
        list_arrow = [
            [p_tt-r, p_cc-r],
        ],
        arrowstyle = 'simple, head_length=2, head_width=2, tail_width=0.2',
        connectionstyle = 'arc3, rad=0.5',
        color = 'tab:blue',
    )
    plot.add_arrow(
        ax,
        list_arrow = [
            [p_ct+t+sx, p_tt-t+sx],
            [p_ct-t-sx, p_cc+t-sx],
            [p_ct+r+sy, p_xx-r1+sy],
        ],
        arrowstyle = 'simple, head_length=2, head_width=2, tail_width=0.2',
        color = 'tab:orange',
    )
    plot.add_arrow(
        ax,
        list_arrow = [
            [p_ct-r+[0,0.1], p_ct-r-[0,0.1]],
        ],
        arrowstyle = 'simple, head_length=2, head_width=2, tail_width=0.2',
        connectionstyle = 'arc3, rad=0.4',
        color = 'tab:orange',
    )
    plot.add_arrow(
        ax,
        list_arrow = [
            [p_cc+t+sx, p_ct-t+sx],
            [p_cc+r+sy, p_xx-r1-sy*3],
        ],
        arrowstyle = 'simple, head_length=2, head_width=2, tail_width=0.2',
        color = 'tab:green',
    )
    plot.add_arrow(
        ax,
        list_arrow = [
            [p_xx-r1+sy*3, p_tt+r-sy],
            [p_xx-r1-sy,   p_ct+r-sy],
            [p_xx-t1, p_cc+r-sy],
        ],
        arrowstyle = 'simple, head_length=2, head_width=2, tail_width=0.2',
        color = 'tab:grey',
    )

    # text in img
    plot.add_text(
        axin0,
        dict_text = {
            (0.02, 0.95): '0.46',
        },
        transform = axin0.transAxes,
        va = 'top',
        ha = 'left',
        color = 'white',
        fontweight = 'bold',
        bbox = dict(boxstyle='round', fc='tab:blue', lw=0)
    )
    plot.add_text(
        axin1,
        dict_text = {
            (0.02, 0.95): '0.82',
        },
        transform = axin1.transAxes,
        va = 'top',
        ha = 'left',
        color = 'white',
        fontweight = 'bold',
        bbox = dict(boxstyle='round', fc='tab:orange', lw=0)
    )
    plot.add_text(
        axin2,
        dict_text = {
            (0.02, 0.95): '0.12',
        },
        transform = axin2.transAxes,
        va = 'top',
        ha = 'left',
        color = 'white',
        fontweight = 'bold',
        bbox = dict(boxstyle='round', fc='tab:green', lw=0)
    )
    plot.add_text(
        axin3,
        dict_text = {
            (0.95, 0.5): '1.13',
        },
        transform = axin3.transAxes,
        va = 'center',
        ha = 'right',
        color = 'white',
        fontweight = 'bold',
        bbox = dict(boxstyle='round', fc='tab:grey', lw=0)
    )

    # text in arrow
    plot.add_text(
        ax,
        dict_text = {
            tuple(p_ct_tt-sx): '0.101',
            tuple(p_xx_tt+sy): '0.333',
            tuple(p_ct-r-3*sx): '0.008',
        },
        va = 'center',
        ha = 'center',
        bbox = dict(boxstyle='round', ec='tab:blue', fc='white')
    )
    plot.add_text(
        ax,
        dict_text = {
            tuple(p_cc_ct-sx): '0.023',
            tuple(p_ct_tt+sx): '0.093',
            tuple(p_xx_ct+sy): '0.690',
            tuple(p_ct-r-sx): '0.008',
        },
        va = 'center',
        ha = 'center',
        bbox = dict(boxstyle='round', ec='tab:orange', fc='white')
    )
    plot.add_text(
        ax,
        dict_text = {
            tuple(p_cc_ct+sx): '0.023',
            tuple(p_xx_cc+sy): '0.085',
        },
        va = 'center',
        ha = 'center',
        bbox = dict(boxstyle='round', ec='tab:green', fc='white')
    )
    plot.add_text(
        ax,
        dict_text = {
            tuple(p_xx_ct-sy): '0.698',
            tuple(p_xx_tt-sy): '0.333',
            tuple(p_xx_cc-sy): '0.085',
        },
        va = 'center',
        ha = 'center',
        bbox = dict(boxstyle='round', ec='tab:grey', fc='white')
    )

def main():

    plot.set_rcparam()
    cm = 1/2.54
    mpl.rcParams['figure.dpi'] = 300
    fig, ax = plt.subplots(figsize=(7*cm, 6*cm))

    run(ax)

    plot.save(
        fig,
        file_save = 'carbonic_flow',
        list_type = ['pdf', 'svg']
    )

    plt.show()

main()

