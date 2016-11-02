#!/bin/bash

if [ -d ${CVMFS_PIPEDIR}/bayeswave ]
then
    export BAYESWAVE_DEPS=${CVMFS_PIPEDIR}/bayeswave
    echo "Setting BAYESWAVE_DEPS=${BAYESWAVE_DEPS}"
else
    echo "CVMFS_PIPEDIR is unset, exiting"
    return
fi

# Default lalsuite version for BayesWave
# You can override this by setting LALSUITE_PREFIX manually
if [ ! -z ${LALSUITE_PREFIX} ]
then
    echo "LALSUITE_PREFIX=${LALSUITE_PREFIX}"
else
    export LALSUITE_PREFIX=${BAYESWAVE_DEPS}/lscsoft/lalsuite-6.38
    echo "Setting LALSUITE_PREFIX=${LALSUITE_PREFIX}"
fi

# -------------------------------------------
#               PYTHON & FRIENDS            #
# -------------------------------------------

# Python
case "$(python --version 2>&1)" in
    *" 2.7"*)
        echo "Local:"
        echo "python --version"
        python --version
        ;;
    *)
        echo "Local:"
        echo "python --version"
        python --version
        echo "using ${BAYESWAVE_DEPS}/non-lsc/Python-2.7.5"

		export PATH=${BAYESWAVE_DEPS}/non-lsc/Python-2.7.5/bin:${PATH}
		export LD_LIBRARY_PATH=${BAYESWAVE_DEPS}/non-lsc/Python-2.7.5/lib:${LD_LIBRARY_PATH}
		export LIBRARY_PATH=${BAYESWAVE_DEPS}/non-lsc/Python-2.7.5/lib:${LD_LIBRARY_PATH}
		export PKG_CONFIG_PATH=${BAYESWAVE_DEPS}/non-lsc/Python-2.7.5/lib/pkgconfig:${PKG_CONFIG_PATH}
		export PYTHONPATH=${BAYESWAVE_DEPS}/non-lsc/Python-2.7.5/lib/python2.7/site-packages:${PYTHONPATH}


        # Numpy
        export PATH=${BAYESWAVE_DEPS}/non-lsc/numpy-1.9.1/bin:${PATH}
        export PYTHONPATH=${BAYESWAVE_DEPS}/non-lsc/numpy-1.9.1/lib/python2.7/site-packages:${PYTHONPATH}

        # Scipy
        export PATH=${BAYESWAVE_DEPS}/non-lsc/scipy-0.12.1/bin:${PATH}
        export PYTHONPATH=${BAYESWAVE_DEPS}/non-lsc/scipy-0.12.1/lib/python2.7/site-packages:${PYTHONPATH}

        # Matplotlib
        export PATH=${BAYESWAVE_DEPS}/non-lsc/matplotlib-1.2.0/bin:${PATH}
        export PYTHONPATH=${BAYESWAVE_DEPS}/non-lsc/matplotlib-1.2.0/lib/python2.7/site-packages:${PYTHONPATH}
        export PYTHONPATH=${BAYESWAVE_DEPS}/non-lsc/basemap-1.0.7/lib/python2.7/site-packages:${PYTHONPATH}

        ;;
esac

# GEOS
export GEOS_DIR=${BAYESWAVE_DEPS}/non-lsc/geos-3.5.0
export LD_LIBRARY_PATH=${GEOS_DIR}/lib:${LD_LIBRARY_PATH}

# Misc dependencies
export PYTHONPATH=${BAYESWAVE_DEPS}/non-lsc/setuptools-28.6.1/lib/python2.7/site-packages:${PYTHONPATH}
export PYTHONPATH=${BAYESWAVE_DEPS}/non-lsc/healpy-1.9.1/lib/python2.7/site-packages:${PYTHONPATH}
export PYTHONPATH=${BAYESWAVE_DEPS}/non-lsc/acor-1.1.1/lib/python2.7/site-packages:${PYTHONPATH}
export PYTHONPATH=${BAYESWAVE_DEPS}/non-lsc/h5py-2.6.0/lib/python2.7/site-packages:${PYTHONPATH}
export PYTHONPATH=${BAYESWAVE_DEPS}/non-lsc/skyarea-0.1/lib/python2.7/site-packages:${PYTHONPATH}
export PYTHONPATH=${BAYESWAVE_DEPS}/non-lsc/M2Crypto-0.25.1/lib/python2.7/site-packages/:${PYTHONPATH}

# -------------------------------------------
#           LALSUITE DEPENDENCIES           #
# -------------------------------------------
# Note: not all of these are required at run-time but all are required for
# building (e.g., swig)


# SWIGPATH (only necessary for building)
export SWIGPATH=${BAYESWAVE_DEPS}/non-lsc/swig-3.0.7
export PATH=${SWIGPATH}/bin:${PATH}

# FFTW
export FFTW=${BAYESWAVE_DEPS}/non-lsc/fftw-3.3.5
export PATH=${FFTW}/bin:${PATH}
export PKG_CONFIG_PATH=${FFTW}/lib/pkgconfig:${PKG_CONFIG_PATH}
export LD_LIBRARY_PATH=${FFTW}/lib:${LD_LIBRARY_PATH}

# GSL
export GSL=${BAYESWAVE_DEPS}/non-lsc/gsl-1.15
export PATH=${GSL}/bin:${PATH}
export PKG_CONFIG_PATH=${GSL}/lib/pkgconfig:${PKG_CONFIG_PATH}
export LD_LIBRARY_PATH=${GSL}/lib:${LD_LIBRARY_PATH}

# LIBFRAME
export LIBFRAME=${BAYESWAVE_DEPS}/lscsoft/libframe-8.30
export PATH=${LIBFRAME}/bin:${PATH}
export PKG_CONFIG_PATH=${LIBFRAME}/lib/pkgconfig:${PKG_CONFIG_PATH}
export LD_LIBRARY_PATH=${LIBFRAME}/lib:${LD_LIBRARY_PATH}

# LIBMETAIO
export LIBMETAIO=${BAYESWAVE_DEPS}/lscsoft/libmetaio-8.4.0
export PATH=${LIBMETAIO}/bin:${PATH}
export PKG_CONFIG_PATH=${LIBMETAIO}/lib/pkgconfig:${PKG_CONFIG_PATH}
export LD_LIBRARY_PATH=${LIBMETAIO}/lib:${LD_LIBRARY_PATH}

# HDF5
export HDF5=$BAYESWAVE_DEPS/non-lsc/hdf5-1.8.16
export PATH=${HDF5}/bin:${PATH}
export LD_LIBRARY_PATH=${HDF5}/lib:${LD_LIBRARY_PATH}
export C_INCLUDE_PATH=${HDF5}/include:${C_INCLUDE_PATH} # for building

# CFITSIO
export CFITSIO=${BAYESWAVE_DEPS}/non-lsc/cfitsio-3.39
export LD_LIBRARY_PATH=${CFITSIO}/lib:${LD_LIBRARY_PATH}
export PKG_CONFIG_PATH=${CFITSIO}/lib/pkgconfig:${PKG_CONFIG_PATH}

# HEALPIX
export HEALPIX=${BAYESWAVE_DEPS}/non-lsc/Healpix_3.30
export LD_LIBRARY_PATH=${HEALPIX}/lib:${LD_LIBRARY_PATH}
export PKG_CONFIG_PATH=${HEALPIX}/lib:${PKG_CONFIG_PATH}


