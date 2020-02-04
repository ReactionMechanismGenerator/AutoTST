# Set base image to centos7
FROM centos:centos7
LABEL maintainers.owner="Sai Krishna Sirumalla <sirumalla.s@northeastern.edu>"
SHELL ["/bin/bash", "-c"]
# Install development tools and upgrade 
RUN yum group install "Development Tools" -y
RUN yum install man-pages man-db man -y
RUN yum install centos-release-scl -y 
RUN yum install devtoolset-8-gcc devtoolset-8-gcc-c++ devtoolset-8-gcc-gfortran -y 
RUN scl enable devtoolset-8 -- bash
RUN gcc --version
RUN whereis gcc  
RUN yum install lapack -y
RUN yum install wget -y
#install for rdkit and cclib 
RUN yum install libXrender libXext fontconfig -y 
# Set up a non-root user, `user`, with a group, `group`
ENV HOME=/home/user
RUN mkdir -p $HOME
RUN groupadd -r group && \
    useradd -r -g group -d $HOME -s /sbin/nologin -c "Default user" user
RUN cp /root/.bashrc $HOME

RUN wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh --directory-prefix=$HOME
RUN /bin/bash $HOME/Anaconda3-2019.10-Linux-x86_64.sh -bp /anaconda3
RUN rm $HOME/Anaconda3-2019.10-Linux-x86_64.sh
ENV PATH /anaconda3/bin:$PATH
# Set up PATH (persistent in final image)
ENV PATH $HOME/.conda/envs/tst_env/bin:/opt/conda/bin:$HOME/AutoTST:$PATH
ENV PYTHONPATH $HOME/AutoTST:$PYTHONPATH
# Set working directory to /user
WORKDIR $HOME
# Get environment file from GitHub and create environment
# Do this here to take better advantage of build caching since the environment does not change as often as AutoTST
ADD  https://raw.githubusercontent.com/ReactionMechanismGenerator/AutoTST/master/environment.yml  /tmp/environment.yml 
RUN echo "unset SUDO_UID SUDO_GID SUDO_USER" >> $HOME/.bashrc
RUN echo ". /anaconda3/etc/profile.d/conda.sh" >> $HOME/.bashrc && \
    echo "conda activate tst_env" >> $HOME/.bashrc
RUN cat /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml
# Clone repo
RUN git clone https://github.com/ReactionMechanismGenerator/AutoTST.git
# Run Tests
RUN cat $HOME/.bashrc

RUN . /anaconda3/etc/profile.d/conda.sh && \
    conda activate tst_env && conda info --envs && conda list rdkit
