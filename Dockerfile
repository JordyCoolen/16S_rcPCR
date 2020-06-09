FROM continuumio/miniconda

# add codebase to docker
ADD ./ /workflow

# Update conda
RUN conda update -y conda

# Install miniconda environment
RUN conda env create -f /workflow/env.yml

# activate by bash 16S env
RUN echo "source activate 16S" > ~/.bashrc
ENV PATH /miniconda/envs/16S/bin:$PATH

# NOT TESTED
# install multiqc
#RUN conda install -c bioconda -c conda-forge multiqc


