FROM rocker/r-ver:latest

# Installer les dépendances système
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Installer le package remotes pour gérer les versions spécifiques
RUN R -e "install.packages('remotes', repos='https://cran.r-project.org')"

# Installer les versions spécifiques des packages R
RUN R -e "remotes::install_version('Seurat', version = '5.2.1', repos = 'https://cran.r-project.org')"
RUN R -e "remotes::install_version('scCustomize', version = '3.0.1', repos = 'https://cran.r-project.org')"
RUN R -e "remotes::install_version('reticulate', version = '1.41.0.1', repos = 'https://cran.r-project.org')"

# Installer la bibliothèque Python anndata
RUN pip3 install anndata==0.11.3   ##########ERRRREUR ICI 

# Définir la commande par défaut
CMD ["R"]