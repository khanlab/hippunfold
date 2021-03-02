FROM khanlab/autotop_deps:v0.2

MAINTAINER alik@robarts.ca

COPY . /src/

RUN apt-get install -y libgraphviz-dev
RUN pip install /src

#pre-download the models here:
ENV HIPPUNFOLD_CACHE_DIR=/opt/hippunfold_cache
RUN hippunfold_download_models

ENTRYPOINT [ "hippunfold" ]

