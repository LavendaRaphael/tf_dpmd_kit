from matplotlib import pyplot as plt
from tf_dpmd_kit import train
from tf_dpmd_kit import plot

cm = 1/2.54

def run_e():
    fig, ax = plt.subplots(figsize=(4.3*cm, 3.7*cm))
    train.dptest_parity_plt(
        ax,
        str_file = 'dptest.e.out',
        int_natoms = 384,
        float_lw = 0.75,
        list_ticks = [-5, 0, 5],
    )
    plot.save(
        fig,
        file_save = 'dptest.e.parity',
        list_type = ['pdf', 'svg']
    )

def run_f():
    fig, ax = plt.subplots(figsize=(4.3*cm, 3.7*cm))
    train.dptest_parity_plt(
        ax,
        str_file = 'dptest.f.out',
        float_lw = 0.75,
        list_ticks = [-5, 0, 5],
    )
    plot.save(
        fig,
        file_save = 'dptest.f.parity',
        list_type = ['pdf', 'svg']
    )

def main():
    plot.set_rcparam()

    run_e()
    run_f()
    plt.show()

main()
