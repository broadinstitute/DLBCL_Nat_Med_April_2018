FROM debian:jessie

RUN apt-get update &&     apt-get install  -y wget tar libgtk2.0-dev libcairo2-dev xvfb xauth xfonts-base 
RUN apt-get install -y r-base-core r-base-dev libxt-dev

WORKDIR /
RUN mkdir /Rscripts/
COPY setup.R /Rscripts/
RUN Rscript /Rscripts/setup.R
COPY *R  /Rscripts/
RUN mkdir src
COPY src /Rscripts/src/
WORKDIR /Rscripts