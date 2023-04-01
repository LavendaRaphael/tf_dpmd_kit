from tf_dpmd_kit import analysis

def main():

    analysis.carbonic_statistic_mean(list_file = [
        '../TT/carbonic/carbonic_statistic.csv',
        '../CT/carbonic/carbonic_statistic.csv',
        '../CC/carbonic/carbonic_statistic.csv',
    ])

main()
