# We will install our software on Alpine Linux
FROM python:3.7-alpine

# Install stuff
RUN apk add --no-cache --virtual build-dependencies git python3-dev libstdc++ g++ gmp-dev cmake make boost-dev gfortran && \
    apk add --no-cache openblas boost-thread boost-filesystem gmp mpfr-dev && \
    apk add --no-cache --virtual test-dependencies bash && \
    apk --no-cache add \
        --repository http://dl-cdn.alpinelinux.org/alpine/edge/testing/ \
        hdf5 py3-numpy py3-h5py py3-pandas && \
    ln -s /usr/include/locale.h /usr/include/xlocale.h && \
    pip3 install --no-cache-dir boto3 && \
    git clone https://github.com/eigenteam/eigen-git-mirror.git && \
    cd /eigen-git-mirror && \
    git checkout tags/3.3.7 && \
    mkdir /eigen-git-mirror/build && \
    cd /eigen-git-mirror/build && \
    cmake -DCMAKE_INSTALL_PREFIX=/usr/local /eigen-git-mirror && \
    make install && \
    cd / && \
    git clone https://github.com/CGAL/cgal.git && \
    cd /cgal && \
    git checkout releases/CGAL-4.14-branch && \
    mkdir /cgal/build && \
    cd /cgal/build && \
    cmake -DCGAL_HEADER_ONLY=ON /cgal && \
    make install && \
    cd / && \
    git clone --recurse-submodules https://github.com/amit112amit/ops-python.git && \
    cd /ops-python && \
    git checkout AWS && \
    mkdir /ops-python/build && \
    cd /ops-python/build && \
    cmake /ops-python/ && \
    make ops && \
    mkdir /simulation && \
    cp /ops-python/build/ops.cpython-37m-x86_64-linux-gnu.so /simulation && \
    cp /ops-python/phasediagramsimulation.py /simulation && \
    cp /ops-python/data/T7.vtk /simulation && \
    cp /ops-python/data/T7_OPS_Asphericity.dat /simulation && \
    cp /ops-python/data/generate_input.py /simulation && \
    cd /simulation && \
    rm -rf /eigen-git-mirror && \
    rm -rf /cgal && \
    rm -rf /ops-python && \
    rm -rf /usr/local/include/CGAL /usr/local/include/eigen3 /usr/local/share/* /usr/local/lib64/cmake && \
    apk del build-dependencies test-dependencies

# Switch to simulation directory
WORKDIR /simulation

# Run the simulation
CMD python3 phasediagramsimulation.py
