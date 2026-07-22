# syntax=docker/dockerfile:1.7

ARG UBUNTU_VERSION=22.04

FROM ubuntu:${UBUNTU_VERSION} AS builder

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        autoconf \
        automake \
        build-essential \
        libbamtools-dev \
        libhts-dev \
        time \
        zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /src/telseq
COPY . .

WORKDIR /src/telseq/src

RUN ./autogen.sh && \
    ./configure \
        --with-bamtools=/usr \
        --with-htslib=/usr \
        --prefix=/opt/telseq

RUN make -j2

RUN ../scripts/test_parallel_synthetic.sh \
        ./Telseq/telseq \
        ./Telseq/telseq \
        ./Test/generate_parallel_fixture

RUN ../scripts/test_parallel_scaling.sh \
        ./Telseq/telseq \
        ./Test/generate_parallel_fixture

RUN make install

FROM ubuntu:${UBUNTU_VERSION} AS runtime

ARG DEBIAN_FRONTEND=noninteractive

LABEL org.opencontainers.image.title="TelSeq Parallel" \
      org.opencontainers.image.description="Indexed multithreaded TelSeq with stock-compatible output" \
      org.opencontainers.image.source="https://github.com/michtrofimov/telseq-parallel" \
      org.opencontainers.image.licenses="GPL-3.0-only"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libbamtools2.5.1 \
        libhts3 && \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/telseq/bin/telseq /usr/local/bin/telseq
COPY LICENSE /usr/local/share/licenses/telseq/LICENSE

RUN telseq --version

WORKDIR /data
ENTRYPOINT ["telseq"]
CMD ["--help"]
