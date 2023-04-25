import matplotlib.pyplot as plt
from tf_dpmd_kit import analysis
from tf_dpmd_kit import plot
import matplotlib as mpl
import numpy as np

def run(ax):

    f_cc = 0.451
    f_ct = 0.747
    f_tt = 0.131
    f_xx = 1.063
    f_cc_ct = 0.114
    f_cc_xx = 0.329
    f_cc_tt = 0.004
    f_ct_cc = 0.076
    f_ct_ct = 0.008
    f_ct_xx = 0.629
    f_ct_tt = 0.025
    f_xx_cc = 0.363
    f_xx_ct = 0.599
    f_xx_tt = 0.101
    f_tt_ct = 0.025
    f_tt_xx = 0.105

    ax.axis('off')

    # horizon fig
    w, h = ax.bbox.width, ax.bbox.height
    dx = 0.27
    dy = 5/8*w/h*dx
    # vertical fig
    dx1 = dy*h/w
    dy1 = dx*w/h

    # right, top
    r = np.array([ dx/2,    0])
    t = np.array([    0,   dy/2])
    r1 = np.array([ dx1/2,    0])
    t1 = np.array([    0,   dy1/2])
    # shift
    sx = np.array([dx/5,   0])
    sy = np.array([   0,dy/5])

    # pos
    p_tt = np.array([1.2*dx, dy/2])
    p_ct = np.array([1.2*dx, 0.5])
    p_cc = np.array([1.2*dx, 1-dy/2])
    p_xx = np.array([1-dx1/2, 0.5])
    # middle
    p_tt_ct = (p_tt+p_ct)/2
    p_ct_cc = (p_ct+p_cc)/2
    p_xx_cc = (p_xx-r1+t1+p_cc+r)/2
    p_xx_ct = (p_xx-r1+p_ct+r)/2
    p_xx_tt = (p_xx-r1-t1+p_tt+r)/2

    # img
    dir_cp  = '/home/faye/research_d/202203_MDCarbonicAcid/server/01.init/H2CO3_CC_H2O_126/plm/'
    dir_tt  = '/home/faye/research_d/202203_MDCarbonicAcid/server/04.md_nvt_velocity/330K/TT/plm/'
    [axin3] = plot.inset_img(
        ax,
        dict_img = {
            dir_cp +   '60147.png': (p_xx[0]-0.5*dx1, p_xx[1]-dy1/2, dx1, dy1),
        },
        dict_spinecolor = {
            dir_cp +   '60147.png': 'tab:purple',
        },
        bool_rot90 = True,
    )
    [axin0, axin1, axin2] = plot.inset_img(
        ax,
        dict_img = {
            dir_tt +'1.100001.png': (p_cc[0]-0.5*dx, p_cc[1]-0.5*dy, dx, dy),
            dir_cp +   '60281.png': (p_ct[0]-0.5*dx, p_ct[1]-0.5*dy, dx, dy),
            dir_tt +'0.003922.png': (p_tt[0]-0.5*dx, p_tt[1]-0.5*dy, dx, dy),
        },
        dict_spinecolor = {
            dir_tt +'1.100001.png': 'tab:blue',
            dir_cp +   '60281.png': 'tab:orange',
            dir_tt +'0.003922.png': 'tab:green',
        }
    )

    headstyle = 'head_length=4, head_width=4'
    # CC
    plot.add_arrow(
        ax,
        list_arrow = [
            [p_cc-t-sx, p_ct+t-sx],
            [p_cc+r+sy, p_xx+t1],
        ],
        arrowstyle = f'simple, {headstyle}, tail_width=0.2',
        color = 'tab:blue',
    )
    plot.add_arrow(
        ax,
        list_arrow = [
            [p_cc-r, p_tt-r],
        ],
        arrowstyle = f'simple, {headstyle}, tail_width=0.2',
        connectionstyle = 'arc3, rad=0.5',
        color = 'tab:blue',
    )
    plot.add_text(
        axin0,
        dict_text = {
            (0.02, 0.95): f_cc,
        },
        transform = axin0.transAxes,
        va = 'top',
        ha = 'left',
        color = 'white',
        fontweight = 'bold',
        bbox = dict(boxstyle='round', fc='tab:blue', lw=0)
    )
    plot.add_text(
        ax,
        dict_text = {
            tuple(p_ct_cc-sx+0.3*sy): f_cc_ct,
            tuple(p_xx_cc+sy): f_cc_xx,
            tuple(p_ct_cc-1.7*r): f_cc_tt,
        },
        va = 'center',
        ha = 'center',
        bbox = dict(boxstyle='round', ec='tab:blue', fc='white')
    )
    # CT
    plot.add_arrow(
        ax,
        list_arrow = [
            [p_ct+t+sx, p_cc-t+sx],
            [p_ct-t-sx, p_tt+t-sx],
            [p_ct+r+sy, p_xx-r1+sy],
        ],
        arrowstyle = f'simple, {headstyle}, tail_width=0.2',
        color = 'tab:orange',
    )
    plot.add_arrow(
        ax,
        list_arrow = [
            [p_ct-r+[0,0.1], p_ct-r-[0,0.1]],
        ],
        arrowstyle = f'simple, {headstyle}, tail_width=0.2',
        connectionstyle = 'arc3, rad=0.4',
        color = 'tab:orange',
    )
    plot.add_text(
        axin1,
        dict_text = {
            (0.02, 0.95): f_ct,
        },
        transform = axin1.transAxes,
        va = 'top',
        ha = 'left',
        color = 'white',
        fontweight = 'bold',
        bbox = dict(boxstyle='round', fc='tab:orange', lw=0)
    )
    plot.add_text(
        ax,
        dict_text = {
            tuple(p_tt_ct-sx+0.3*sy): f_ct_tt,
            tuple(p_ct_cc+sx-0.3*sy): f_ct_cc,
            tuple(p_xx_ct+sy): f_ct_xx,
            tuple(p_ct-r-sx): f_ct_ct,
        },
        va = 'center',
        ha = 'center',
        bbox = dict(boxstyle='round', ec='tab:orange', fc='white')
    )

    # TT
    plot.add_arrow(
        ax,
        list_arrow = [
            [p_tt+t+sx, p_ct-t+sx],
            [p_tt+r+sy, p_xx-r1-sy*3],
        ],
        arrowstyle = f'simple, {headstyle}, tail_width=0.2',
        color = 'tab:green',
    )
    plot.add_text(
        axin2,
        dict_text = {
            (0.02, 0.95): f_tt,
        },
        transform = axin2.transAxes,
        va = 'top',
        ha = 'left',
        color = 'white',
        fontweight = 'bold',
        bbox = dict(boxstyle='round', fc='tab:green', lw=0)
    )
    plot.add_text(
        ax,
        dict_text = {
            tuple(p_tt_ct+sx-0.3*sy): f_tt_ct,
            tuple(p_xx_tt+sy): f_tt_xx,
        },
        va = 'center',
        ha = 'center',
        bbox = dict(boxstyle='round', ec='tab:green', fc='white')
    )
    # HCO3
    plot.add_arrow(
        ax,
        list_arrow = [
            [p_xx-r1+sy*3, p_cc+r-sy],
            [p_xx-r1-sy,   p_ct+r-sy],
            [p_xx-t1, p_tt+r-sy],
        ],
        arrowstyle = f'simple, {headstyle}, tail_width=0.2',
        color = 'tab:purple',
    )
    plot.add_text(
        axin3,
        dict_text = {
            (0.95, 0.5): f_xx,
        },
        transform = axin3.transAxes,
        va = 'center',
        ha = 'right',
        color = 'white',
        fontweight = 'bold',
        bbox = dict(boxstyle='round', fc='tab:purple', lw=0)
    )
    plot.add_text(
        ax,
        dict_text = {
            tuple(p_xx_ct-sy): f_xx_ct,
            tuple(p_xx_cc-sy): f_xx_cc,
            tuple(p_xx_tt-sy): f_xx_tt,
        },
        va = 'center',
        ha = 'center',
        bbox = dict(boxstyle='round', ec='tab:purple', fc='white')
    )

    plot.add_text(
        ax,
        dict_text = {
            (0.9, 0.95): r'Frequency (ns$^{-1}$)',
        },
        va = 'top',
        ha = 'right',
    )

def main():

    plot.set_rcparam()
    cm = 1/2.54
    mpl.rcParams['figure.dpi'] = 300
    fig, ax = plt.subplots(figsize=(8.6*cm, 6*cm))

    run(ax)

    plot.save(
        fig,
        file_save = 'carbonic_flow',
        list_type = ['pdf', 'svg']
    )

    plt.show()

main()

