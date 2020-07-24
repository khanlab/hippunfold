FROM khanlab/hippocampal_autotop:v0.1

MAINTAINER alik@robarts.ca


RUN pip install pybids

RUN mkdir -p /app
COPY . /app

ENTRYPOINT ["/app/run.py"]

