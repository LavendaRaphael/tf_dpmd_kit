import json
from tf_dpmd_kit import analysis

def main():
    with open('timelong.json', 'r') as fp:
        dict_timelong = json.load(fp)
        print(dict_timelong)

    analysis.carbonic_statistic(
        time_tot = dict_timelong['timelong(ps)'],
        file_save = 'carbonic_statistic.csv',
    )

main()
