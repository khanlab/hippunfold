FROM khanlab/hippocampal_autotop:reg-fix

MAINTAINER alik@robarts.ca


RUN pip install pybids

RUN mkdir -p /app
COPY . /app

ENTRYPOINT ["/app/run.py"]

