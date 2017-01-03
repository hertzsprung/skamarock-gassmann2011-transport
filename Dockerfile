FROM ubuntu:16.04
MAINTAINER James Shaw <js102@zepler.net>
RUN apt-get update && apt-get install -y \
    python3-numpy \
    python3-pip \
    git \
    gnuplot
RUN pip3 install pytest==3.0.5
RUN groupadd -r transport && \
    useradd -r -m -d /home/transport -s /sbin/nologin -g transport transport && \
    chown -R transport:transport /home/transport
USER transport
RUN git clone https://github.com/hertzsprung/skamarock-gassmann2011-transport.git /home/transport/high-order-transport
WORKDIR /home/transport/high-order-transport
RUN git checkout high-order-transport-datumedge
RUN make -j all
