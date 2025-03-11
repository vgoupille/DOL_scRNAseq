# Use Miniconda as the base image
FROM continuumio/miniconda3:25.1.1-2 
    # latest

# Set a working directory
WORKDIR /app

# Copy the environment.yml file into the container
COPY environment.yml /app/environment.yml

# Create and activate the Conda environment
RUN conda env create -f /app/environment.yml && \
    conda clean --all --yes

# Activate the environment by default
SHELL ["conda", "run", "-n", "myenv_StarSolo", "/bin/bash", "-c"]

# Set the default command to start JupyterLab
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root"]