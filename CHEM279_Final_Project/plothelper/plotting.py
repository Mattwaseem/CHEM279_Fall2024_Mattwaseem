import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

mass_H = 1.6735575e-27

# Read the CSV file into a DataFrame
df = pd.read_csv('log.txt')
df = df[df['t'] > 1733972000]
df['T_v'] = df['T_v']*1e20
df['f_LJ'] = df['f_LJ']/mass_H*1e20
df['f_C'] = df['f_C']/mass_H*1e10
df['E_LJ'] = df['E_LJ']
grouped = df.groupby('t').mean()
print(grouped['E_LJ'])
plt.plot(grouped.index, grouped['E_LJ'])
# plt.plot(df['t'], df['E_LJ'])
# Plot each column individually
# for column in grouped.columns:
#     plt.plot(grouped.index, grouped[column], label=column)


# plt.yscale('log')
# Add title and labels
plt.title('Average Energy Over Time')
plt.xlabel('Time')
plt.ylabel('Energy (J)')

# Show legend
plt.legend()

# Display the plot
plt.show()


df['T_v'] = df['T_v']*mass_H/1e20/1e3*6.023e23
df['f_LJ'] = df['f_LJ']*mass_H/1e20/1e3*6.023e23
df['f_C'] = df['f_C']*mass_H/1e10/1e3*6.023e23
df['E_LJ'] = df['E_LJ']*mass_H/1e10/1e3*6.023e23
df.to_csv('converted_to_kJ_per_mol.txt', index=False)