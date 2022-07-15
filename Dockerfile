FROM ubuntu:20.04


WORKDIR /gusApp
WORKDIR /gusApp/project_home

RUN apt-get -qq update --fix-missing

# Installing Software
RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y \
    wget \
    openjdk-8-jdk \
    python3-pip \
    git \
    trimmomatic=0.39+dfsg-1

RUN apt-get install -y libbz2-dev liblzma-dev

RUN pip3 install cython==0.29.30
RUN pip3 install numpy==1.23.1
RUN pip3 install metaphlan==3.0.14


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
