# ------------------------------------------------------------------------------
# Single-stage Docker container image

FROM python:3.11-slim

ARG PANDAS_VER=1.5.3
ARG MATPLOTLIB_VER=3.7.1
ARG SEABORN_VER=0.12.2
ARG RICH_VER=13.3.2
ARG IPYTHON_VER=8.11.0
ARG DEEPTOOLS_VER=3.5.1

RUN apt-get update \
    && apt-get install --assume-yes --no-install-recommends \
       build-essential python3-dev zlib1g-dev \
    && pip install \
       pandas==${PANDAS_VER} \
       rich==${RICH_VER} \
       matplotlib==${MATPLOTLIB_VER} \
       seaborn==${SEABORN_VER} \
       ipython==${IPYTHON_VER} \
       deepTools==${DEEPTOOLS_VER} \
       --no-cache-dir \
    && apt-get clean \
    && apt-get purge --assume-yes \
    && rm -rf /var/lib/apt/lists*

CMD /usr/local/bin/ipython

