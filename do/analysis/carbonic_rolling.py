import pandas as pd

def run(
    int_window: int = 10,
    str_file: str = 'carbonic_state.csv',
    str_save: str = 'carbonic_rolling.csv',
):
    df_data = pd.read_csv(str_file)
    print(df_data)
    df_data = df_data.where(df_data.notnull(), 0)
    df_new = df_data.rolling(int_window, min_periods=1, center=True, step=int_window).mean()
    df_new = df_new.where( df_new!=0, None)
    print(df_new)
    df_new.to_csv(str_save, index=False)

#'''
run(
    int_window = 200,
)
#'''
'''
run(
    int_window = 100,
)
#'''
'''
run(
    int_window = 800,
)
#'''
