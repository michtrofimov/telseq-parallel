# syntax=docker/dockerfile:1.7

ARG UBUNTU_VERSION=22.04

FROM ubuntu:${UBUNTU_VERSION} AS builder

ARG DEBIAN_FRONTEND=noninteractive
ARG BAMTOOLS_VERSION=2.4.0
ARG BAMTOOLS_SHA256=f1fe82b8871719e0fb9ed7be73885f5d0815dd5c7277ee33bd8f67ace961e13e

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        autoconf \
        automake \
        build-essential \
        ca-certificates \
        cmake \
        curl \
        zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /build

RUN curl -fsSL \
        "https://github.com/pezmaster31/bamtools/archive/refs/tags/v${BAMTOOLS_VERSION}.tar.gz" \
        -o bamtools.tar.gz && \
    echo "${BAMTOOLS_SHA256}  bamtools.tar.gz" | sha256sum -c - && \
    tar -xzf bamtools.tar.gz && \
    cmake \
        -S "bamtools-${BAMTOOLS_VERSION}" \
        -B bamtools-build \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/opt/bamtools && \
    cmake --build bamtools-build --parallel "$(nproc)" && \
    cmake --install bamtools-build

WORKDIR /src/telseq
COPY . .

WORKDIR /src/telseq/src

RUN ./autogen.sh && \
    ./configure \
        --with-bamtools=/opt/bamtools \
        --prefix=/opt/telseq && \
    make -j"$(nproc)" && \
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
        libstdc++6 \
        zlib1g && \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/telseq/bin/telseq /usr/local/bin/telseq
COPY --from=builder /opt/bamtools/lib/bamtools/libbamtools.so* /usr/local/lib/
COPY LICENSE /usr/local/share/licenses/telseq/LICENSE

RUN ldconfig && telseq --version

WORKDIR /data
ENTRYPOINT ["telseq"]
CMD ["--help"]
