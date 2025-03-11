# Utiliser une version spécifique de l'image de base pour garantir la reproductibilité
FROM rocker/r-ver:latest

# Mettre à jour le système et installer les dépendances de base de manière sécurisée
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    python3 \
    python3-pip \
    python3-dev \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Installer Miniconda (version spécifique) pour la gestion des environnements
RUN curl -sSL https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh -o /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -f -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    /opt/conda/bin/conda init bash

# Ajouter Conda à PATH
ENV PATH=/opt/conda/bin:$PATH

# Créer un fichier environment.yml avec toutes les dépendances
COPY environment.yml /tmp/environment.yml

# Créer l'environnement Conda depuis le fichier environment.yml
RUN conda env create -f /tmp/environment.yml && \
    rm /tmp/environment.yml

# Activer l'environnement Conda 'myenv' par défaut
RUN echo "source activate myenv" > ~/.bashrc

# Installer IRkernel pour rendre R disponible dans Jupyter
RUN R -e "install.packages('IRkernel')"
RUN R -e "IRkernel::installspec(user = FALSE)"  # Enregistre le kernel R pour Jupyter

# Installer ipykernel dans chaque environnement Conda pour que Python soit disponible dans Jupyter
RUN conda install -n myenv ipykernel
RUN conda install -n otherenv ipykernel

# Ajouter les environnements Conda comme kernels dans Jupyter
RUN /opt/conda/envs/myenv/bin/python -m ipykernel install --user --name=myenv --display-name "Python (myenv)"
RUN /opt/conda/envs/otherenv/bin/python -m ipykernel install --user --name=otherenv --display-name "Python (otherenv)"

# Définir la commande par défaut pour lancer Jupyter Notebook
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--allow-root"]