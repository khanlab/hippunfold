FROM khanlab/hippocampal_autotop:latest

MAINTAINER alik@robarts.ca

# base directory for Hippocampal_AutoTop
RUN mkdir -p /app
COPY . /app

ENTRYPOINT ["/app/run.py"]

