FROM khanlab/hippunfold_deps:v0.5.1

MAINTAINER alik@robarts.ca

COPY . /src/

# avoid pre-downloading the models to make for lighter container
# ENV HIPPUNFOLD_CACHE_DIR=/opt/hippunfold_cache
ENV PYTHONNOUSERSITE=1

#install hippunfold and imagemagick (for reports)
RUN pip install --no-cache-dir /src && \
    apt install -y graphviz && \
    wget https://imagemagick.org/archive/binaries/magick && \
    mv magick /usr/bin && chmod a+x /usr/bin/magick 
    

ENTRYPOINT [ "hippunfold" ]

