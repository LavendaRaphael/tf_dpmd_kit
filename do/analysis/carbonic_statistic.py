import json
from tf_dpmd_kit import analysis

def main():
    with open('timelong.json', 'r') as fp:
        time_tot = json.load(fp)['timelong(ps)']
        print(time_tot)
    with open('volume.json', 'r') as fp:
        volume = json.load(fp)['volume(ang3)']
        print(volume)

    analysis.carbonic_statistic(
        time_tot = time_tot,
        volume = volume,
        file_save = 'carbonic_statistic.csv',
    )

main()
