import pandas as pd

# Input parameters
df = '../../gather_data/data/data.csv'  # Path to data file

df = pd.read_csv(df)

print(df)
