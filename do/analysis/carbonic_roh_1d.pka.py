import pandas as pd
from tf_dpmd_kit import plm

def main():

    list_key = []
    list_df = []
    for temperature in [290, 310, 330, 350]:
        df_data = pd.read_csv(f'../{temperature}K/carbonic/carbonic_roh_1d.csv', index_col='roh0(ang)')
        r = df_data.index
        deltag = df_data[(r>3.5) & (r<4.5)].min()-df_data.min()
        pka = plm.deltag_to_pka(deltag=deltag, temperature=temperature)
        df_tmp = pd.DataFrame({'deltag(kcal/mol)': deltag, 'pka': pka})
        df_tmp.index.name = 'state'
        list_key.append(temperature)
        list_df.append(df_tmp)
    df_save = pd.concat(list_df, keys=list_key, names=['temperature(K)'])
    file_save = f'carbonic_roh_1d.pka.csv'
    print(file_save)
    print(df_save)
    df_save.to_csv(file_save)

main()
