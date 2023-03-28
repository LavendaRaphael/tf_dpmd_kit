import pandas as pd

run(
    file_data = 'carbonic.csv',
    tup_select = (200000, 5000000),
    file_save = 'carbonic.product.csv',
)
run(
    file_data = 'carbonic_state.csv',
    tup_select = (200000, 5000000),
    file_save = 'carbonic_state.product.csv',
)
