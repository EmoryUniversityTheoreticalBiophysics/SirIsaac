FROM ubuntu:latest
RUN apt-get update && apt-get -y update
FROM python:2
ARG DEBIAN_FRONTEND=noninteractive
RUN mkdir app
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends sudo apt-utils && \
    apt-get install -y --no-install-recommends openssh-server \
        python-dev python-numpy python-pip python-virtualenv python-scipy \
        gcc gfortran libopenmpi-dev openmpi-bin openmpi-common openmpi-doc binutils && \
    apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
WORKDIR /app
RUN useradd --create-home --home-dir /home/docker --shell /bin/bash docker
RUN usermod -a -G sudo docker
RUN echo "docker ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers
RUN wget http://www.open-mpi.org/software/ompi/v4.1/downloads/openmpi-4.1.1.tar.gz
RUN tar xzvf openmpi-4.1.1.tar.gz
WORKDIR openmpi-4.1.1
RUN ./configure --prefix=/usr/local/openmpi4.1.1
RUN make all
RUN make install
RUN MPI_DIR=/usr/local/openmpi4.1.1
ENV LD_LIBRARY_PATH="${MPI_DIR}/lib:${LD_LIBRARY_PATH}"
ENV PATH="/usr/local/openmpi4.1.1/bin:${PATH}"
RUN git clone https://github.com/mpi4py/mpi4py.git ./mpi4py.git
WORKDIR mpi4py.git
RUN python setup.py build --mpicc=/usr/local/openmpi4.1.1/bin/mpicc
RUN python setup.py install
WORKDIR /app
RUN apt-get -y update &&  apt -y install libopenmpi-dev

RUN pip install numpy
RUN pip install jupyter notebook
RUN git clone -b mpi4py-convesion https://github.com/diging/SloppyCell.git
RUN pip install -e SloppyCell/
RUN git clone -b pypar-check  https://github.com/diging/SirIsaac.git
RUN pip install -e SirIsaac/
CMD ["jupyter", "notebook", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root"]
