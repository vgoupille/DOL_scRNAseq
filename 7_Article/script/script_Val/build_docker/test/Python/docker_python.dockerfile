# Utiliser Miniconda comme base
FROM continuumio/miniconda3:25.1.1-2

# Définir le dossier de travail
WORKDIR /app

# Copier le fichier d’environnement Conda
COPY environment.yml /tmp/environment.yml

# Créer l’environnement Conda
RUN conda env create -f /tmp/environment.yml && \
    rm /tmp/environment.yml && \
    conda clean --all --yes

# Définir le shell par défaut avec l’environnement activé
SHELL ["conda", "run", "-n", "myenv_Python", "/bin/bash", "-c"]

# Exposer le port pour JupyterLab
EXPOSE 8888

# Définir la commande par défaut pour JupyterLab
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root"]