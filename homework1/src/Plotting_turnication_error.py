#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np

# Loading data from my truncation_error text file
data = np.loadtxt('../truncation_error_Problem2.txt', delimiter=',', skiprows=1)

# parse and find h values and errors
h_values = data[:, 0]
forward_errors = data[:, 1]
central_errors = data[:, 2]

# ignore duplicates
unique_h_values, indices = np.unique(h_values, return_index=True)
forward_errors = forward_errors[indices]
central_errors = central_errors[indices]

#plot
plt.figure()
plt.loglog(unique_h_values, forward_errors, 'o-', label="Forward Difference Error", color='blue')  # Removed redundant marker
plt.loglog(unique_h_values, central_errors, 'x-', label="Central Difference Error", color='orange')  # Removed redundant marker
plt.xlabel("Step size (h)")
plt.ylabel("Truncation Error")
plt.legend()

# Save the plot as a PNG file
plt.savefig("truncation_error_plot-v2.png")
print("Plot saved to truncation_error_plot-v2.png")

# Uncomment to display the plot
# plt.show()
