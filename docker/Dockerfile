# ==============================
# Dockerfile para hormone-analysis
# Basado en un archivo .yml @EXPLICIT
# ==============================

# Imagen base con Miniconda
FROM continuumio/miniconda3

# Variables
ENV ENV_NAME=hormone-analysis

# Copiar el archivo explicit .yml
COPY environment.lock.linux-64.yml /tmp/environment.lock.linux-64.yml

# Instalar dependencias del sistema necesarias
RUN apt-get update && \
    apt-get install -y curl wget sudo git build-essential && \
    rm -rf /var/lib/apt/lists/*

# Crear entorno conda a partir de un spec file @EXPLICIT
RUN conda create --name $ENV_NAME --file /tmp/environment.lock.linux-64.yml && \
    conda clean --all -y

# Configurar Conda y usar el entorno por defecto
SHELL ["conda", "run", "-n", "hormone-analysis", "/bin/bash", "-c"]

# Comando por defecto
CMD ["/bin/bash"]
