# Use Miniconda as the base image
FROM continuumio/miniconda3:25.1.1-2

# Set a working directory
WORKDIR /app

# Install system dependencies for R
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    git \
    r-base \
    r-base-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy environment files for each environment
COPY python_environment.yml /tmp/python_environment.yml
COPY r_environment.yml /tmp/r_environment.yml
COPY starsolo_environment.yml /tmp/starsolo_environment.yml
COPY install.R /tmp/install.R
COPY renv.lock /renv.lock
COPY renv /renv

# Create all three Conda environments with distinct configurations
RUN conda env create -f /tmp/python_environment.yml -n myenv_Python && \
    conda env create -f /tmp/r_environment.yml -n myenv_R && \
    conda env create -f /tmp/starsolo_environment.yml -n myenv_StarSolo && \
    conda clean --all --yes && \
    rm /tmp/python_environment.yml /tmp/r_environment.yml /tmp/starsolo_environment.yml

# Install R packages using renv
RUN Rscript /tmp/install.R && rm /tmp/install.R

# Setup Python kernel for each environment
RUN conda run -n myenv_Python python -m pip install ipykernel && \
    conda run -n myenv_Python python -m ipykernel install --user --name=myenv_Python --display-name "Python (myenv_Python)" && \
    conda run -n myenv_R python -m pip install ipykernel && \
    conda run -n myenv_R python -m ipykernel install --user --name=myenv_R --display-name "Python (myenv_R)" && \
    conda run -n myenv_StarSolo python -m pip install ipykernel && \
    conda run -n myenv_StarSolo python -m ipykernel install --user --name=myenv_StarSolo --display-name "Python (myenv_StarSolo)"

# Setup R kernel for R environment
RUN conda run -n myenv_R R -e "install.packages('IRkernel')" && \
    conda run -n myenv_R R -e "IRkernel::installspec(user = FALSE, name = 'myenv_R', displayname = 'R (myenv_R)')"

# Create a script to configure aliases for activating each environment
RUN echo '#!/bin/bash\n\
alias activate_python="conda activate myenv_Python"\n\
alias activate_r="conda activate myenv_R"\n\
alias activate_starsolo="conda activate myenv_StarSolo"\n\
echo "Available environments:"\n\
echo " - myenv_Python (use: activate_python)"\n\
echo " - myenv_R (use: activate_r)"\n\
echo " - myenv_StarSolo (use: activate_starsolo)"\n\
' > /etc/profile.d/conda_env_aliases.sh && \
chmod +x /etc/profile.d/conda_env_aliases.sh

# Make the environment activation script run on container startup
RUN echo "source /etc/profile.d/conda_env_aliases.sh" >> ~/.bashrc

# Expose port for JupyterLab
EXPOSE 8888

# Start JupyterLab by default
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root"]