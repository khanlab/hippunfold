FROM khanlab/autotop_deps:v0.2

MAINTAINER alik@robarts.ca

COPY . /src/

RUN pip install /src

ENTRYPOINT [ "hippunfold" ]

