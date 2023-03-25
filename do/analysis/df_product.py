import pandas as pd

def run(
    str_file: str,
    tup_select: tuple,
    str_save: str,
):

    print(str_file)
    df_data = pd.read_csv(str_file)
    print(df_data)

    df_data = df_data.iloc[tup_select[0]:tup_select[1]]
    print(str_save)
    print(df_data)
    df_data.to_csv(str_save, index=False)

run(
    str_file = 'carbonic.csv',
    tup_select = (200000, 5000000),
    str_save = 'carbonic.product.csv',
)
run(
    str_file = 'carbonic_state.csv',
    tup_select = (200000, 5000000),
    str_save = 'carbonic_state.product.csv',
)
