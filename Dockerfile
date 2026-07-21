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
        zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /src/telseq
COPY . .

WORKDIR /src/telseq/src

RUN ./autogen.sh && \
    ./configure \
        --with-bamtools=/usr \
        --prefix=/opt/telseq && \
    make -j2 && \
    ./Test/generate_parallel_fixture /tmp/parallel-fixture.bam && \
    ./Telseq/telseq -t 1 /tmp/parallel-fixture.bam \
        > /tmp/telseq-serial.tsv 2> /tmp/telseq-serial.log && \
    ./Telseq/telseq -t 22 /tmp/parallel-fixture.bam \
        > /tmp/telseq-parallel.tsv 2> /tmp/telseq-parallel.log && \
    cmp /tmp/telseq-serial.tsv /tmp/telseq-parallel.tsv && \
    make install

FROM ubuntu:${UBUNTU_VERSION} AS runtime

ARG DEBIAN_FRONTEND=noninteractive

LABEL org.opencontainers.image.title="TelSeq Parallel" \
      org.opencontainers.image.description="Indexed multithreaded TelSeq with stock-compatible output" \
      org.opencontainers.image.source="https://github.com/michtrofimov/telseq-parallel" \
      org.opencontainers.image.licenses="GPL-3.0-only"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libbamtools2.5.1 && \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/telseq/bin/telseq /usr/local/bin/telseq
COPY LICENSE /usr/local/share/licenses/telseq/LICENSE

RUN telseq --version

WORKDIR /data
ENTRYPOINT ["telseq"]
CMD ["--help"]
