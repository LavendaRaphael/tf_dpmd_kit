import pandas as pd
import matplotlib.pyplot as plt
import json

def main():
    file_list = [
        '../8/log',
        '../9/log',
        '../10/log',
        '../11/log',
        '../12/log',
        '../13/log',
        '../14/log',
        '../15/log',
        '../16/log',
        '../17/log',
        '../18/log',
        '../19/log',
    ]
    #file_list = ['test.log']
    data = []
    for file in file_list:
        print(file)
        with open(file) as f:
            for line in f:
                if line[:5] == '  - [':
                    data.append(line[5:-4].split(', '))
                elif line[:9] == 'keywords:':
                    columns = line[12:-5].split("', '")
                    print(columns)
    data = pd.DataFrame(data, columns=columns, dtype='float64')
    print(data)
    v_mean = data['Volume'].mean()
    with open('volume.json', 'w') as fp:
        json.dump({'volume(ang3)': v_mean}, fp, indent=4)
    print(v_mean)
    fig, ax = plt.subplots()
    for x in [10000]:
        df_x = data.rolling(x, min_periods=1, center=True, step=x).mean()
        ax.plot(df_x['Step'], df_x['Volume'], label=x)
    ax.hlines(v_mean, data['Step'][0], data['Step'].iloc[-1], color='tab:green', label=f'{v_mean:.2f}')
    ax.legend()
    plt.show()

main()
