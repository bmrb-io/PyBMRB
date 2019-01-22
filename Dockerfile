FROM jupyter/scipy-notebook

RUN pip install plotly
RUN pip install pynmrstar
RUN pip install pybmrb
