FROM ubuntu:16.04

RUN mkdir -p /app
WORKDIR /app

RUN apt-get update
RUN apt-get install -y make g++ libboost-all-dev libconfig++-dev python-pip tcsh
RUN apt-get install -y autoconf automake autotools-dev gnuplot libxml2 libz-dev
RUN apt-get install -y libxslt1-dev

COPY requirements.txt .
RUN pip install -r requirements.txt

WORKDIR /app

COPY . /app

RUN mkdir -p /usr/cbs/bio/src
RUN mkdir /scratch
RUN mv netphos-3.1.Linux.tar.Z /usr/cbs/bio/src
WORKDIR /usr/cbs/bio/src
RUN uncompress netphos-3.1.Linux.tar.Z
RUN tar xf netphos-3.1.Linux.tar
RUN sed -i 's/setenv\s\+ARCH.*$/setenv ARCH i386/' /usr/cbs/bio/src/ape-1.0/ape
RUN ln -s /usr/cbs/bio/src/ape-1.0/ape /usr/bin/netphos

WORKDIR /app

RUN ./autogen.sh && ./configure && make && make install


