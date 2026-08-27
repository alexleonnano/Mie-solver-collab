import pandas as pd

def load(path):
	# Read using a whitespace-aware separator to handle mixed tabs/spaces
	df = pd.read_csv(path, sep=r"\s+", engine="python", comment='#', header=0, names=["wl", "n", "k"]) 
	# Coerce columns to numeric 
	for col in ["wl", "n", "k"]:
		df[col] = pd.to_numeric(df[col], errors='coerce')

	print("Material file dtypes:\n", df.dtypes)		# Print the data types for verification
	# Error handling for whitespaces or NaN values in the data
	if df[["wl", "n", "k"]].isnull().any(axis=1).any():
		print("Warning: NaNs detected in material file. Rows with issues:")
		print(df[df[["wl", "n", "k"]].isnull().any(axis=1)].to_string(index=False))

	# Drop rows with missing data for clean data
	df = df.dropna(subset=["wl", "n", "k"]) 
	return df