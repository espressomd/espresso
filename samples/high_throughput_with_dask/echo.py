import dask_espresso as de

data = de.get_data_from_stdin()
data.update(processed=True)

print(de.encode_transport_data(data))
