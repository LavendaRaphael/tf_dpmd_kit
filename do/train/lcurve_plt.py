from tf_dpmd_kit import train
from tf_dpmd_kit import plot
from matplotlib import pyplot as plt

def main():
    plot.set_rcparam()
    cm = 1/2.54
    
    fig, ax = plt.subplots(figsize = (8.6*cm, 5*cm))
    train.plt_lcurve(ax)
    
    cm = 1/2.54
    plot.save(
        fig,
        file_save = 'lcurve.png',
    )
    
    plt.show()

main()
