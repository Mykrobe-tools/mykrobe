FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
ARG PKGS=" \
  build-essential \
  ca-certificates \
  git \
  gnupg \
  libbz2-dev \
  libcurl4-gnutls-dev \
  libssl-dev \
  liblzma-dev \
  python-is-python3 \
  python3-pip \
  tzdata \
  wget \
  zlib1g-dev \
"

RUN apt update \
    && apt-get install --no-install-recommends -y $PKGS \
    && update-ca-certificates -f

# install mongodb - https://www.mongodb.com/docs/manual/tutorial/install-mongodb-on-ubuntu/

RUN wget -qO - https://www.mongodb.org/static/pgp/server-5.0.asc | apt-key add -
RUN echo "deb [ arch=amd64,arm64 ] https://repo.mongodb.org/apt/ubuntu focal/mongodb-org/5.0 multiverse" | tee /etc/apt/sources.list.d/mongodb-org-5.0.list
RUN apt-get update && apt-get install -y mongodb-org && rm -rf /var/lib/apt/lists/*

ARG PROJECT="mykrobe"
COPY . "/${PROJECT}"

# install mccortex
WORKDIR "/tmp"

RUN git clone --recursive -b geno_kmer_count https://github.com/phelimb/mccortex \
    && cd mccortex \
    && make MAXK=31 \
    && cp bin/mccortex31 /${PROJECT}/src/mykrobe/cortex/ \
    && rm -rf /tmp/mccortex

# install mykrobe
WORKDIR "/${PROJECT}"
RUN python -m pip install requests && python -m pip install . -vv

# download panels
RUN mykrobe panels update_metadata \
    && mykrobe panels update_species all