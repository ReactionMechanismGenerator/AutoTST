# Set base image to anaconda
FROM continuumio/anaconda3
# Install make and gcc
RUN apt-get update && apt-get install -y make gcc g++
# Set up PATH (persistent in final image)
ENV PATH /home/tst-admin/.conda/envs/tst_env/bin:/opt/conda/bin:/home/tst-admin/AutoTST:$PATH
ENV PYTHONPATH /home/tst-admin/AutoTST:$PYTHONPATH
# Set working directory to /tst-admin
WORKDIR /home/tst-admin
# Get environment file from GitHub and create environment
# Do this here to take better advantage of build caching since the environment does not change as often as AutoTST
ADD  https://raw.githubusercontent.com/ReactionMechanismGenerator/AutoTST/master/environment.yml  /tmp/environment.yml
RUN echo "unset SUDO_UID SUDO_GID SUDO_USER" >> /home/tst-admin/.bashrc
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> /home/tst-admin/.bashrc && \
    echo "conda activate tst_env" >> /home/tst-admin/.bashrc
RUN cat /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml
# Clone repos
RUN git clone https://github.com/ReactionMechanismGenerator/AutoTST.git
# Compile RMG
RUN cd AutoTST && make 
