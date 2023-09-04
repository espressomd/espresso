import dask_espresso as de
import numpy as np


def test_encode_decode(): 
    data = {"a": 1, "b": np.random.random((3, 3))}
    encoded = de.encode_transport_data(data)
    decoded = de.decode_transport_data(encoded)
    for k, v in decoded.items():
        assert np.all(data[k]) == np.all(v)
    assert list(data.keys()) == list(decoded.keys())


def test_espresso_runner():
    data = {"hello": "world"}
    result = de.dask_espresso_task("python3", "echo.py", **data).compute()
    assert result == {"hello": "world", "processed": True}
