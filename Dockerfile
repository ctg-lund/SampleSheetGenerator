FROM python:3.9-slim-bullseye
RUN dasdasdkjnkdfsadsaf

# Upgrade pip before installing other packages.
RUN python -m pip install --upgrade pip

# Numpy must apparently be installed before scikit-bio.
RUN pip install numpy

# Copy only the requirements file and install them all.
COPY requirements.txt /requirements.txt
RUN pip install -r requirements.txt

# Copy the entire application into an "app" folder, and use that as workdir
# going forward.
COPY / /app
WORKDIR /app

# Make Gunicorn the entrypoint program, and then feed it appropriate arguments.
ENTRYPOINT [ "gunicorn" ]
CMD [ "-b", ":8082", "wsgi:application" ]