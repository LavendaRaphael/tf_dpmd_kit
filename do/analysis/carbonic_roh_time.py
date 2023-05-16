import pandas as pd
import matplotlib.pyplot as plt

def main():

    df = pd.read_csv('carbonic.csv')
    fig, ax = plt.subplots()
    ax.plot(df['frame'], df['roh0(ang)'])
    ax.plot(df['frame'], df['roh1(ang)'])
    plt.show()

main()
