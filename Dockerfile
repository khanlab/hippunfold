FROM khanlab/autotop_deps:v0.4.2

MAINTAINER alik@robarts.ca

COPY . /src/

RUN cd /src && pip install -r requirements.txt .

#pre-download the models here:
ENV HIPPUNFOLD_CACHE_DIR=/opt/hippunfold_cache
RUN hippunfold_download_models

ENTRYPOINT [ "hippunfold" ]

