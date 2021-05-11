FROM gitpod/workspace-full
RUN sudo apt-get install --yes libgsl0-dev

USER gitpod
ENV SHELL /bin/bash
RUN wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
RUN ./bin/micromamba create -n varlociraptor samtools -c bioconda -c conda-forge
RUN yes | ./bin/micromamba shell init -s bash -p /home/gitpod/micromamba
