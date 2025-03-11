# Utiliser une version spécifique de R
FROM rocker/r-ver:4.3.3

# Installer des dépendances système pour R et Python
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Installer Miniconda avec une version spécifique
RUN curl -sSL https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh -o /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -f -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    /opt/conda/bin/conda init bash

# Ajouter Conda à PATH
ENV PATH=/opt/conda/bin:$PATH

# Copier et installer les dépendances Conda
COPY myenvR.yml /tmp/myenvR.yml
RUN conda env create -f /tmp/myenvR.yml && rm /tmp/myenvR.yml

# Activer l’environnement Conda par défaut
RUN echo "source activate myenv_R" > ~/.bashrc

# Copier les fichiers pour installer les packages R
COPY install.R /tmp/install.R
COPY renv.lock /renv.lock
COPY renv /renv

# Exécuter le script pour installer les packages R et Bioconductor
RUN Rscript /tmp/install.R && rm /tmp/install.R

# Installer ipykernel pour rendre Python accessible dans Jupyter
RUN /opt/conda/envs/myenv_R/bin/python -m pip install --upgrade pip && \
    /opt/conda/envs/myenv_R/bin/python -m pip install ipykernel && \
    /opt/conda/envs/myenv_R/bin/python -m ipykernel install --user --name=myenv_R --display-name "Python (myenv_R)"

# Ajouter le kernel R à Jupyter
RUN R -e "install.packages('IRkernel')" && \
    R -e "IRkernel::installspec(user = FALSE, name = 'myenv_R', displayname = 'R (myenv_R)')"

# # Ajouter l'environnement Conda à Jupyter
# RUN /opt/conda/envs/myenv_R/bin/python -m ipykernel install --user --name=myenv_R --display-name "Python (myenv_R)"

# # Installer IRkernel et le configurer pour Jupyter Lab
# RUN R -e "install.packages('IRkernel')" && \
#     R -e "IRkernel::installspec(user = FALSE, name = 'myenv_R', displayname = 'R (myenv_R)')"

# Définir la commande par défaut pour Jupyter
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--allow-root"]