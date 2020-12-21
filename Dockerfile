FROM ubuntu:12.04
LABEL MAINTAINER="Pilar Natividad"
LABEL EMAIL="pnatividad@incliva.es"

RUN apt-get update
RUN apt-get install -y   gcc make g++ bowtie libz-dev libncurses5-dev libncursesw5-dev libtbb-dev python
 

RUN mkdir -p /usr/local/src/MapSplice_multi_threads_2.0.1.9/
COPY  MapSplice_multi_threads_2.0.1.9/ /usr/local/src/MapSplice_multi_threads_2.0.1.9/
WORKDIR /usr/local/src/MapSplice_multi_threads_2.0.1.9/
RUN make

ENV PATH="/root/.local/bin:/usr/local/src/MapSplice_multi_threads_2.0.1.9/bin:${PATH}"
RUN echo PATH="/root/.local/bin:/usr/local/src/MapSplice_multi_threads_2.0.1.9/bin":${PATH} >> /etc/profile

VOLUME [ "/media/scratch/" ]
WORKDIR /media/scratch

CMD [ "/bin/bash" ]
