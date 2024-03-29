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
    libbz2-dev \
    liblzma-dev \
  && rm -rf /var/lib/apt/lists/*

RUN pip3 install cython==0.29.30
RUN pip3 install numpy==1.23.1
RUN pip3 install metaphlan==3.0.14

COPY bin/config.py /gusApp/project_home/

RUN export HUMANN_GIT_COMMIT_SHA=13e8c7910d9aff4cabb58df19aa652a3c20e101b \
    && git clone https://github.com/wbazant/humann.git \
    && cd humann \
    && git checkout $HUMANN_GIT_COMMIT_SHA \
    && cd humann \
    && rm config.py \
    && mv /gusApp/project_home/config.py . \
    && cd .. \
    && python3 setup.py install

WORKDIR /gusApp/project_home

RUN export KNEADDATA_GIT_COMMIT_SHA=969a06fdec01ac142be38d518ed942aafc8535d8 \
    && git clone https://github.com/rdemko2332/kneaddata.git \
    && cd kneaddata \
    && git checkout $KNEADDATA_GIT_COMMIT_SHA \
    && python3 setup.py install

RUN humann_config --update database_folders nucleotide /humann_databases/chocophlan
RUN humann_config --update database_folders protein /humann_databases/uniref
RUN humann_config --update database_folders utility_mapping /humann_databases/utility_mapping

WORKDIR /gusApp/tools

RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-ubuntu64.tar.gz \
    && tar -xvzf sratoolkit.3.0.0-ubuntu64.tar.gz

ENV PATH=$PATH:/gusApp/tools/sratoolkit.3.0.0-ubuntu64/bin

RUN mkdir -p /root/.ncbi
RUN printf '/LIBS/GUID = "%s"\n' `uuidgen` > /root/.ncbi/user-settings.mkfg
RUN printf '/libs/cloud/report_instance_identity = "true"\n' >> /root/.ncbi/user-settings.mkfg
RUN printf '/libs/cloud/accept_aws_charges = "false"\n/libs/cloud/accept_gcp_charges = "false"\n' >> /root/.ncbi/user-settings.mkfg

COPY bin/* /usr/local/bin/

RUN ln -s /usr/bin/python3 /usr/bin/python