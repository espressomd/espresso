if (HDF5_FOUND)
    configure_file(h5md.py h5md.py COPYONLY)
    configure_file(__init__.py __init__.py COPYONLY)
endif(HDF5_FOUND)
