# We will install our software on Ubuntu
FROM archlinux/base

# Update the package database
RUN pacman -Syu --noconfirm

# Install compilers and build tools
# Install C++ dependencies
# Install Python and Python libraries
RUN pacman -S --noconfirm \
    make gcc gcc-fortran cmake git openblas boost \
    pybind11 qt5-base cgal eigen vtk \
    python python-numpy python-h5py python-pandas \
    python-boto3 aws-cli

# Clean-up all pacman cache files
RUN pacman -Scc --noconfirm

# Create a non-privileged user and switch to that user
RUN useradd -U -m amit && chown amit /home/amit
USER amit

# Pull our OPS code and build it
RUN cd /home/amit && \
    git clone --recurse-submodules https://github.com/amit112amit/ops-python.git && \
    cd /home/amit/ops-python && \
    mkdir /home/amit/ops-python/build && \
    cd /home/amit/ops-python/build && \
    cmake /home/amit/ops-python/ && \
    make ops

RUN mkdir /home/amit/simulation && \
    cp /home/amit/ops-python/build/ops.cpython-37m-x86_64-linux-gnu.so /home/amit/simulation && \
    cp /home/amit/ops-python/phasediagramsimulation.py /home/amit/simulation && \
    cp /home/amit/ops-python/data/T7.vtk /home/amit/simulation && \
    cp /home/amit/ops-python/data/T7_OPS_Asphericity.dat /home/amit/simulation && \
    cp /home/amit/ops-python/data/generate_input.py /home/amit/simulation && \
    cd /home/amit/simulation && \
    python generate_input.py

# Switch to simulation directory
WORKDIR /home/amit/simulation

# Run the simulation
CMD python phasediagramsimulation.py
