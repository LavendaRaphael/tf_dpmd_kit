import pandas as pd

def main():

    list_file = [
        '../0/carbonic.csv',
        '../1/carbonic.csv',
        '../2/carbonic.csv',
        '../3/carbonic.csv',
        '../4/carbonic.csv',
        '../5/carbonic.csv',
        '../6/carbonic.csv',
        '../7/carbonic.csv',
        '../8/carbonic.csv',
        '../9/carbonic.csv',
        '../10/carbonic.csv',
    ]
    file_save = 'carbonic.csv'

    df = pd.DataFrame()
    for file in list_file:
        print(file)
        df_tmp = pd.read_csv(file, index_col=0)
        print(df_tmp)
        df = pd.concat([df, df_tmp], ignore_index=True)

    df.index.name = 'frame'
    print(file_save)
    print(df)
    df.to_csv(file_save)

main()
