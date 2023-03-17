import pandas as pd

def run(
    int_window: int = 10,
):
    df_data = pd.read_csv('carbonic_state.csv')
    print(df_data)
    df_data = df_data.where(df_data.notnull(), 0)
    df_new = df_data.rolling(int_window, min_periods=1, center=True, step=int_window).mean()
    df_new = df_new.where( df_new!=0, None)
    print(df_new)
    df_new.to_csv('carbonic_rolling.csv', index=False)

run(
    int_window = 800,
    #int_window = 80,
)
