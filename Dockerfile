FROM khanlab/autotop_deps:v0.1

MAINTAINER alik@robarts.ca


COPY . /src/

RUN pip install /src

ENV AUTOTOP_DIR /src

ENTRYPOINT [ "hippunfold" ]

