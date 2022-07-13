FROM ubuntu:20.04


WORKDIR /gusApp/project_home


# Installing Software
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y \
    wget \
    openjdk-8-jdk \
    python3-pip \
    git \
    trimmomatic=0.39+dfsg-1 \
 && rm -rf /var/lib/apt/lists/*

RUN pip3 install cython==0.29.30
RUN pip3 install numpy==1.23.1

RUN export HUMANN_GIT_COMMIT_SHA=13e8c7910d9aff4cabb58df19aa652a3c20e101b \
    && git clone https://github.com/wbazant/humann.git \
    && cd humann \
    && git checkout $HUMANN_GIT_COMMIT_SHA \
    && python3 setup.py install

WORKDIR /gusApp/project_home

RUN export KNEADDATA_GIT_COMMIT_SHA=379a7c8a141350d57a988af458d232a21e5bc5f6 \
    && git clone https://github.com/wbazant/kneaddata.git \
    && cd kneaddata \
    && git checkout $KNEADDATA_GIT_COMMIT_SHA \
    && python3 setup.py install
