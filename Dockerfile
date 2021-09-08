FROM khanlab/autotop_deps:v0.4.1

MAINTAINER alik@robarts.ca

COPY . /src/

RUN pip install /src

#pre-download the models here:
ENV HIPPUNFOLD_CACHE_DIR=/opt/hippunfold_cache
RUN hippunfold_download_models

ENTRYPOINT [ "hippunfold" ]

