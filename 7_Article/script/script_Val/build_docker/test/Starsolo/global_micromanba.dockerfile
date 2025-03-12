# Use Micromamba as the base image
FROM mambaorg/micromamba:1.5.6

# Set working directory
WORKDIR /app

USER root

# Install system dependencies for R
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    git \
    r-base \
    r-base-dev \
    && rm -rf /var/lib/apt/lists/*

# Switch back to micromamba user for environment creation
USER micromamba

# Copy environment files
COPY --chown=micromamba:micromamba python_environment.yml /tmp/python_environment.yml
COPY --chown=micromamba:micromamba r_environment.yml /tmp/r_environment.yml
COPY --chown=micromamba:micromamba starsolo_environment.yml /tmp/starsolo_environment.yml
COPY --chown=micromamba:micromamba install.R /tmp/install.R
COPY --chown=micromamba:micromamba renv.lock /renv.lock
COPY --chown=micromamba:micromamba renv /renv

# Create all three environments with distinct configurations
RUN micromamba create -f /tmp/python_environment.yml -n myenv_Python && \
    micromamba create -f /tmp/r_environment.yml -n myenv_R && \
    micromamba create -f /tmp/starsolo_environment.yml -n myenv_StarSolo && \
    micromamba clean --all --yes && \
    rm /tmp/python_environment.yml /tmp/r_environment.yml /tmp/starsolo_environment.yml

# Install R packages using renv
RUN /usr/bin/Rscript /tmp/install.R && rm /tmp/install.R

# Setup Python kernel for each environment
RUN micromamba run -n myenv_Python python -m pip install ipykernel && \
    micromamba run -n myenv_Python python -m ipykernel install --user --name=myenv_Python --display-name "Python (myenv_Python)" && \
    micromamba run -n myenv_R python -m pip install ipykernel && \
    micromamba run -n myenv_R python -m ipykernel install --user --name=myenv_R --display-name "Python (myenv_R)" && \
    micromamba run -n myenv_StarSolo python -m pip install ipykernel && \
    micromamba run -n myenv_StarSolo python -m ipykernel install --user --name=myenv_StarSolo --display-name "Python (myenv_StarSolo)"

# Setup R kernel for R environment
RUN micromamba run -n myenv_R /usr/bin/R -e "install.packages('IRkernel', repos='https://cran.rstudio.com/')" && \
    micromamba run -n myenv_R /usr/bin/R -e "IRkernel::installspec(user = FALSE, name = 'myenv_R', displayname = 'R (myenv_R)')"

# Create a script to configure aliases for activating each environment
USER root
RUN echo '#!/bin/bash\n\
alias activate_python="micromamba activate myenv_Python"\n\
alias activate_r="micromamba activate myenv_R"\n\
alias activate_starsolo="micromamba activate myenv_StarSolo"\n\
echo "Available environments:"\n\
echo " - myenv_Python (use: activate_python)"\n\
echo " - myenv_R (use: activate_r)"\n\
echo " - myenv_StarSolo (use: activate_starsolo)"\n\
' > /etc/profile.d/micromamba_env_aliases.sh && \
chmod +x /etc/profile.d/micromamba_env_aliases.sh

# Make the environment activation script run on container startup
RUN echo "source /etc/profile.d/micromamba_env_aliases.sh" >> /home/micromamba/.bashrc

# Switch back to micromamba user for running the container
USER micromamba

# Expose port for JupyterLab
EXPOSE 8888

# Start JupyterLab by default (need to specify which environment to use)
CMD ["micromamba", "run", "-n", "myenv_Python", "jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root"]