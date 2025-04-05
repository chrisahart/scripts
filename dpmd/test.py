import numpy as np

# Example initialization (replace with your actual data)
coordinates = np.random.rand(2700, 3, 120)  # Example array

pos_index = [0, 1, 1, 3, 4, 5, 6, 6, 7]  # Example indices with duplicates

# Step 1: Identify Unique Indices while maintaining order

unique_indices = []
seen = set()

for idx in pos_index:
    if idx not in seen:
        unique_indices.append(idx)
        seen.add(idx)

# Step 2: Use these unique indices to filter the coordinates array
print(unique_indices)
coordinates_no_duplicates = coordinates[unique_indices]

# Output the shape of the filtered array to verify

print(coordinates_no_duplicates.shape)